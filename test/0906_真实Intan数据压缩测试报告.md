# 真实 Intan 16-bit 神经数据压缩测试报告

**日期**：2026-09-06　|　**测试对象**：neuralink-codec 改造版（encode16w / decode16w，去量化 + 符号表 ±511）

## 算法与代码来源

### 基线代码库

- **来源**：GitHub `hxrdxkxvxd/neuralink-codec`（Neuralink 2024 压缩挑战赛的参赛提交，本地已 clone 到 `D:/bcidata/neuralink-codec/`，原版 encode.c / decode.c 完整保留）
- **挑战赛官方成绩**：无损压缩率 **3.4575×**（本地复现一致；评测脚本为逐字节 diff）
- **代码规模**：仅两个 C 文件（encode.c + decode.c，约 330 行），无第三方依赖，gcc -O3 直接编译，无 LICENSE 文件（使用需注意授权风险）

### 算法结构（原版）

```
输入样本 → 8 阶 SSLMS 自适应线性预测（去掉时间相关性，残差逼近白噪声）
        → 4 上下文自适应概率模型
        → 二进制算术/Range 编码器 → 码流
```

- 单通道流式编码（`desc.channels = 1` 写死），128 通道即 128 个独立实例，天然适合 FPGA 分时复用
- 预测器 + 熵编码的组合与学术主流（Buccino 2023 的 WavPack、Dresden DPCM2+熵编码）同属一个范式，差异只在预测器阶数与熵编码选择

### 为何必须改造（原版对 16-bit 数据天生有损）

原版在挑战赛数据上是逐字节无损的，但前提是**官方数据本身是 10-bit 信号放在 ×64 格点上**（encode.su 论坛确认，发现者为 phoboslab）。对 RHD2164 的满 16-bit 原始计数，原版两处代码会直接破坏数据：

| 位置 | 原版代码 | 问题 | 改造版 |
|---|---|---|---|
| `encode.c` L217 | `quant(v) = floor(v/64)` | 编码前丢弃低 6 bit，量化误差 ±37 LSB（比信号噪声 σ≈25 LSB 还大） | 改为恒等映射 `quant(v)=v` |
| `decode.c` L152 | `dequant(q) = round(q×64.0616+31.03)` | 用最小二乘拟合的近似直线还原，仅对格点数据精确可逆 | 改为恒等映射 `dequant(v)=v` |
| `encode.c`/`decode.c` 宏 | `MAX_VAL 256`（符号表 ±127，逃逸符号 255） | 16-bit 残差频繁超限走 escape 路径，合成宽带测试中 CR 崩到 1.1× | `MAX_VAL 1024`（符号表 ±511，逃逸符号 1023） |

改造版源码：`D:/bcidata/neuralink-codec/encode16w.c`、`decode16w.c`（相对原版仅改动上述 4 处，算法主体未动）。原版与"仅去量化"中间版一并保留可对照复现。

## 结论（要点前置）

在**真实在体神经数据**（麻醉大鼠新皮层，Intan RHD-2000 系统采集，128 通道 × 20 kHz × 16-bit）上：

| 指标 | 结果 |
|---|---|
| 无损性 | **128/128 通道逐 bit 无损（MaxErr = 0）**，两个独立时段均验证 |
| 压缩率 | **2.37–2.38×**（6.72–6.74 bit/样本） |
| 折算 128ch × 30 kHz 带宽 | **25.8–25.9 Mbps**（原始 61.44 Mbps）——落在方案文档预测的 17–31 Mbps 无损窗口内，可接 Wi-Fi 6 链路 |
| 原版对照组 | 同数据上 MaxErr = 33–35 LSB（**有损**）；在 Intan 官方 64ch 样本上另有 5/64 通道**码流失步**（数据损坏不可解） |

## 数据来源

### 数据集 1（主证据）：Horváth et al. 2021, Scientific Data 8:180
- **内容**：20 只麻醉大鼠（ketamine/xylazine）新皮层自发宽频记录，109 段，共 7126 个排序单元
- **硬件**：Intan RHD-2000 系统——**一片 64 通道头stage（RHD2164）+ 两片 32 通道头stage = 128 通道**；宽频 0.1–7500 Hz，20 kS/s，16-bit
- **本次使用**：`Rat06_Insertion2_Depth1.nwb`（3.37 GB，15.1 min），提取两段各 60 s（第 300–360 s、700–760 s）全 128 通道
- **数据画像**：每通道 std ≈ 460–546 LSB（0.195 μV/LSB，与 RHD2164 一致）；慢波 LFP 主导 + spike；相邻通道相关 0.998（密集探针共享 LFP，属正常在体特征）
- **许可**：CC BY 4.0；下载自 HuggingFace 镜像 `rokaijano/rat_cortical_128ch`
- **原文**：doi:10.1038/s41597-021-00970-3

### 数据集 2（辅助）：Intan 官方 sampledata.rhd
- 来源：github.com/Intan-Technologies/load-rhd-notebook-python（官方示例文件）
- 64 通道放大器，20 kS/s，1.49 s，RHD 芯片真实输出
- 注意：该文件为评估板演示性质（各通道近共模、互相关 0.989），非在体信号，仅作芯片输出格式验证
- 另测试 Mendeley w767nnk5wh（Orbán，859 MB .rhd）：16 通道互相关 = 1.000（各通道完全相同），判定非有效多通道在体数据，**弃用**

## 测试方法

1. NWB (HDF5/gzip) → h5py 提取 `acquisition/wideband_multichannel_recording/data`（int16 原始计数，conversion = 0.195 μV）
2. 每通道写 16-bit mono wav → 分别运行三版编码器（gcc -O3 编译）
3. 解码后与原始逐样本比对 MaxErr；按文件字节数计 CR

> 注：官方读取器 `load_file` 会把数据乘 0.195 μV 换算成浮点；本测试用 `np.round(x/0.195)` 精确恢复 int16 原始计数后再编码，保证逐 bit 无损判定在 ADC 计数域进行。

## 结果明细

### 大鼠在体数据（Rat06，128ch × 60s × 20kHz × 16bit，两段）

| 版本 | 改动 | CR（中位/总体） | bit/样本 | 折算 128ch@30kHz | MaxErr | 失步 |
|---|---|---|---|---|---|---|
| 原版 | quant(v)=v/64 | 11.82× / 11.61× | 1.38 | 5.3 Mbps | **35 LSB（有损）** | 0/128 |
| 仅去量化 | 恒等 quant | 2.40× / 2.38× | 6.72 | 25.8 Mbps | 0 | 0/128 |
| **去量化+宽表** | 恒等 quant + 符号表 ±127→±511 | **2.40× / 2.38×** | **6.73** | **25.8 Mbps** | **0** | 0/128 |

- 第二时段（700–760 s）：宽表版 MaxErr=0，CR 2.37×，25.9 Mbps —— 与第一段一致
- 宽表版 CR 分布（128 通道）：min 2.18 / 25% 2.26 / 中位 2.40 / 75% 2.51 / max 2.74
- 本数据 LFP 重、残差小（窄表与宽表结果相同）；宽表的意义在宽带/spike 密集场景（合成 AP 带测试中窄表崩到 1.1×，宽表 2.3×），保留宽表配置更稳妥

### Intan 官方 64ch 样本（RHD2164 类芯片输出）

| 版本 | CR | MaxErr | 备注 |
|---|---|---|---|
| 原版 | 10.86×（可解的 59 通道） | 33 LSB | **5/64 通道码流失步，数据损坏** |
| 去量化+宽表 | 2.28× | 0 | 64/64 无损，折算 26.9 Mbps |

## 复现路径（2026-09-06 已全部迁移至 D 盘）

- 改造版源码与 exe：`D:/bcidata/neuralink-codec/`（encode16w.c / decode16w.c，含原版与中间版本）
- 大鼠数据与全部中间产物：`D:/bcidata/`（Rat06_Insertion2_Depth1.nwb、rat128_raw.npy、ratwav/、rres/、ratwav2/、rres2/）
- 官方样本与产物：`D:/bcidata/intan-data/`（sampledata.rhd、ch64wav/、vendor_h5py/）
- 旧合成数据归档：`D:/bcidata/neuralink-codec-archive/`（testdata/、testdata2/、testdata3/、results/）
- markitdown 依赖库：`D:/bcidata/vendor_markitdown/`（docx 转换时 `PYTHONPATH=D:/bcidata/vendor_markitdown`）

## 残留局限（答辩口径）

1. 大鼠新皮层 ≠ 犬类嗅球：嗅球呼吸节律耦合更强、嗅探期 LFP 更大，CR 可能略低于 2.38×；但同为宽带胞外记录，量级可信
2. 麻醉慢波状态下 LFP 相关性高、预测器友好；清醒自由活动数据（伪迹更多）CR 预计再降 10–20%，25.8 Mbps 仍有 Wi-Fi 6 余量（TCP 45–60 Mbps）
3. 30 kHz 与 20 kHz 对算法无本质差异（预测 + 熵编码与采样率无关），折算按 bit/样本等比放大

## 环境说明（2026-09-06 20:00 更新）

C 盘曾满至 0 字节，本次已将全部大文件迁移至 D 盘并校验（详见"复现路径"）；C 盘现剩余约 30 GB。后续大数据一律落 `D:/bcidata/`。
