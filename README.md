# dog-bci-lossless-codec

针对 Intan RHD2164（128 通道 × 30 kS/s × 16 bit）无线 BCI 链路的** 16-bit 无损神经信号压缩器**。基于 Neuralink 2024 压缩挑战赛参赛代码改造而成。

## 代码来源与授权

- **上游**：GitHub [`hxrdxkxvxd/neuralink-codec`](https://github.com/hxrdxkxvxd/neuralink-codec)（Neuralink 压缩挑战赛提交，官方评测无损 3.4575×，本地复现一致）。原版完整保留在 [`upstream/`](upstream/) 目录。
- **⚠️ 上游无 LICENSE 文件**。本仓库为私有研究用途；如需公开，必须先处理上游代码的授权问题。
- 改造版相对原版仅改动 4 处（见下表），算法主体未动。

## 为何必须改造

上游在挑战赛数据上逐字节无损，但官方数据是 **10-bit 信号放在 ×64 格点上**。对满 16-bit 的 RHD2164 原始 ADC 计数，原版是**有损**的：

| 位置 | 原版 | 问题 | 改造版（encode16w / decode16w） |
|---|---|---|---|
| `encode.c` L217 | `quant(v) = floor(v/64)` | 丢弃低 6 bit，误差 ±37 LSB（大于信号噪声 σ≈25 LSB） | 恒等映射 |
| `decode.c` L152 | `dequant(q) = round(q×64.0616+31.03)` | 拟合直线，仅对格点数据可逆 | 恒等映射 |
| 符号表宏 | `MAX_VAL 256`（±127） | 16-bit 残差超限走 escape，CR 崩到 ~1.1× | `MAX_VAL 1024`（±511） |

`encode16.c / decode16.c` 为"仅去量化"中间版本（符号表未扩），保留作对照。

## 真实数据实测结果（2026-09-06）

数据：Horváth et al. 2021（Sci. Data 8:180），麻醉大鼠新皮层在体记录，Intan RHD-2000 采集（RHD2164 类，128ch × 20 kS/s × 16-bit），CC BY 4.0。

| 指标 | 结果 |
|---|---|
| 无损性 | **128/128 通道逐 bit 无损（MaxErr = 0）**，两段独立 60 s 时段验证 |
| 压缩率 | **2.37–2.40×**（6.72–6.74 bit/样本） |
| 折算 128ch × 30 kHz | **25.8–25.9 Mbps**（原始 61.44 Mbps），落在 17–31 Mbps 无损窗口，可接 Wi-Fi 6 |
| 原版对照 | 同数据 MaxErr 33–35 LSB（有损）；Intan 官方 64ch 样本上 5/64 通道码流失步 |
| 速度 | 单线程 ≈6.5 M 样本/s（gcc -O3），128ch×30kHz = 3.84 M 样本/s，单核足够 |

完整报告（含数据来源、方法、局限）：[`docs/0906_真实Intan数据压缩测试报告.md`](docs/0906_真实Intan数据压缩测试报告.md)（附一~附七）

## 仓库结构

| 目录 | 内容 |
|---|---|
| `encode16*.c` / `decode16*.c` | 16-bit 无损压缩器（Neuralink 改造版，见上） |
| [`drice/`](drice/) | C 语言 drice（差分 Rice-Golomb）编解码器——FPGA 移植首选结构，每通道 1 个减法器 |
| [`rtl_demo/`](rtl_demo/) | 首个 RTL 演示核（Verilog）：帧捕获 → 1 阶差分 → 17bit zigzag → 固定 K Rice，含 Python 黄金模型与 Icarus 仿真对拍脚本（黄金模型自检已通过；RTL 仿真待安装 Icarus 后运行 `run_sim.py`） |
| [`bench/`](bench/) | 两套真实数据基准测试脚本（rat128 / dandi001539），含全部实测结果 JSON 与日志 |
| [`docs/`](docs/) | 测试报告（0906/0907）、数据集调研、EDA 图表、ZYNQ 移植规划 |
| [`docs/papers/`](docs/papers/) + [`docs/文献阅读笔记/`](docs/文献阅读笔记/) | 参考文献 PDF 与阅读笔记（首篇：Nicolelis 2014 Nature Methods 无线大规模记录——本项目对标文献，含与压缩路线的对比分析） |
| [`upstream/`](upstream/) | Neuralink 挑战赛原版代码（未改动） |

## 基准测试结论摘要（2026-09-06/07，两套真实在体数据）

| 类别 | 方案 | CR | 备注 |
|---|---|---|---|
| 无损（rat128, 128ch） | 本仓库 encode16w | **2.38×** | 逐 bit 无损，25.9 Mbps @128ch×30kHz |
| 无损（dandi48, 48ch LFP） | T5 双亲LS+FLAC | **2.041×** | 需定点化 |
| 无损（dandi48） | FLAC 立体声配对 | 1.990× | 零自研代码，hub/PC 层适用 |
| 无损（dandi48） | T1+drice | 1.949× | FPGA 亲和（每通道 1 减法器） |
| 近无损（dandi48） | T1 残差+JPEG-LS NEAR=2 | **2.582×** | worst 20µV（残差域量化沿树深放大） |
| 近无损（dandi48） | SZ3 abs=4µV | 2.252× | 误差逐点硬保证 |

**诚实结论**：无重复接线时跨通道预测的真实增益仅 ~4%（dandi48 的 3.3× 是 tetrode 双线冗余红利）；单通道 1 阶差分熵 9.30 bit 的上限下，14 种单通道编码器均未显著超过 FLAC。详见报告附七。

后续 FPGA 移植规划（正点原子领航者 ZYNQ 7020）：[`docs/0907_压缩管线移植领航者ZYNQ规划.md`](docs/0907_压缩管线移植领航者ZYNQ规划.md)

## 算法结构

```
样本 → 8 阶 SSLMS 自适应线性预测（残差白化）→ 4 上下文自适应概率模型 → Range 编码
```

单通道流式（128 通道 = 128 个独立实例），适合 FPGA Hub 分时复用。

## 编译与使用

```bash
gcc -O3 -o encode16w encode16w.c -lm -std=c99
gcc -O3 -o decode16w decode16w.c -lm -std=c99

# 单通道 16-bit mono wav（30 kS/s）进，压缩码流出
./encode16w input.wav output.bw
./decode16w output.bw reconstructed.wav
```

上游原版构建/评测见 [`upstream/`](upstream/)（`build.sh` / `eval.sh`，针对挑战赛数据格式）。
