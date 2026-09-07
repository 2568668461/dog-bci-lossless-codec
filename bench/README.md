# bench/ — 压缩基准测试脚本

针对犬类嗅觉 BCI（Intan RHD2164，128ch × 30 kS/s × 16-bit，61.44 Mbps 无线化）的压缩算法选型实测。

**结论先行**：全部脚本与数据结论的完整叙述见 [`../docs/0906_真实Intan数据压缩测试报告.md`](../docs/0906_真实Intan数据压缩测试报告.md)（附一~附七）。

## 目录

| 子目录 | 数据集 | 内容 |
|---|---|---|
| `rat128/` | Horváth et al. 2021（Sci. Data 8:180）麻醉大鼠 128ch×20kHz 在体记录 | 多轮算法横评（round3 系列）、drice 近无损、熵界、EDA、Neuralink 3× 复现分析 |
| `dandi001539/` | DANDI:001539 清醒大鼠 OB/CA1/PFC LFP 48ch×1.5kHz | 第二基准：图像类/科学压缩器扫描、T1/T5 跨通道预测树、FLAC 立体声、分脑区分解 |

## 环境依赖

- Python 3.12（系统版），`soundfile`、`h5py`、`numpy`、`scipy`、`matplotlib`
- `imagecodecs`（sz3/aec/jpegxl/jpeg2k/htj2k/lerc/zfp 等），Windows 下用 wheel 安装到 `--target vendor` 目录后 `PYTHONPATH` 引入
- ffmpeg / flac（CLI，部分脚本调用）
- C 编译器（drice/Neuralink codec 对照实验）

## 数据获取（不入库，体积大）

- **rat128**：Horváth 2021 数据，Mendeley Data 公开（CC BY 4.0），`chunk_download.sh` / `rat_extract.py` 为下载与提取脚本
- **dandi001539**：DANDI archive `sub-CS39_ses-06.nwb` 等约 2GB，`chunk_download_d.sh` + `extract.py` 为下载与提取脚本

## 路径约定

脚本内为绝对路径 `D:/bcidata/...`（当时的实测工作区布局）。在新机器复现时需要：
1. 重建 `D:/bcidata`（或全局替换脚本内路径）
2. 按上述数据获取方式下载数据
3. `PYTHONPATH` 指向含 imagecodecs 的 vendor 目录

## 脚本索引

### rat128/
| 脚本 | 用途 |
|---|---|
| `round3_test.py` | 第三轮横评：drice/FLAC/WavPack/lzma 等 |
| `round3b_lzma.py` / `round3c_flac.py` | lzma 变体 / FLAC 参数扫描专项 |
| `drice_nl_test.py` | drice 近无损（NEAR 量化）4 方法实测 |
| `near_lossless_test.py` | 近无损方案横评 |
| `headroom_test.py` | 16-bit 余量/饱和统计 |
| `entropy_bound.py` | 1 阶/2 阶差分熵界计算 |
| `why_challenge_3x.py` | Neuralink 挑战赛 3.46× 在 16-bit 真实数据上失效原因分析 |
| `rat_test.py` | rat128 数据基础验证 |
| `eda_full.py` / `eda_supp.py` | rat128 数据集 EDA（输出见 `../docs/eda_figs/`） |

### dandi001539/
| 脚本 | 用途 |
|---|---|
| `extract.py` | 从 NWB 提取 48ch int16（LSB=1µV） |
| `bench.py` | 第一轮：FLAC/WavPack/通用库基线 |
| `bench2.py` | 第二轮：图像类/科学压缩器/航天标准扫描 |
| `bench3.py` | 跨通道预测树（T1/T1d/T4/T5、FLAC 立体声、SZ3/AEC） |
| `bench3b.py` | 分脑区分解：去冗余 vs 真实相关增益 |
| `special_tests.py` | 专项（重复通道、行序效应等） |
| `bench*_results.json` / `bench*_log.txt` | 实测结果与日志（可直接引用数字） |

## ⚠️ 方法学提醒（详见报告附七 §四）

- 跨通道预测必须用**树形结构**（根通道独立编码），成对互预测存在循环依赖不可解码
- NWB 的 electrodes 表是脑区映射唯一权威，不可按通道顺序硬编码
- JPEG-LS NEAR 作用在残差域时误差沿树深度链式放大（worst ≈ NEAR × 链深）
