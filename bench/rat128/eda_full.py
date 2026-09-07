# -*- coding: utf-8 -*-
"""rat128_raw.npy 全面数据集分析：统计 + 质量指标 + 可视化"""
import numpy as np, json, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
rcParams["axes.unicode_minus"] = False

D = "D:/bcidata"
FIG = os.path.join(D, "eda_figs")
os.makedirs(FIG, exist_ok=True)

mm = np.load(os.path.join(D, "rat128_raw.npy"), mmap_mode="r")
print("shape:", mm.shape, "dtype:", mm.dtype)
nch, ns = mm.shape

# 分块读取转 float64 计算（避免整块 300MB int16 -> float64 = 1.2GB 一次性驻留）
x = np.asarray(mm, dtype=np.float64)  # 128*1.2M*8B ~= 1.2GB，可行
# 免转副本方案保险起见监控内存
import gc

# ---------- 1. 逐通道描述统计 ----------
stats = {}
stats["mean"] = x.mean(axis=1)
stats["std"] = x.std(axis=1, ddof=1)
stats["min"] = x.min(axis=1)
stats["max"] = x.max(axis=1)
q = np.quantile(x, [0.01, 0.25, 0.5, 0.75, 0.99], axis=1)
stats["q01"], stats["q25"], stats["q50"], stats["q75"], stats["q99"] = q
# 偏度/峰度（整体样本量足够，用聚合公式）
m1 = stats["mean"]
m2 = ((x - m1[:, None]) ** 2).mean(axis=1)
m3 = ((x - m1[:, None]) ** 3).mean(axis=1)
m4 = ((x - m1[:, None]) ** 4).mean(axis=1)
skew = m3 / m2 ** 1.5
kurt = m4 / m2 ** 2 - 3.0

# ---------- 2. 数据质量指标 ----------
int16_min, int16_max = -32768, 32767
sat_lo = (x <= int16_min + 2).sum(axis=1)
sat_hi = (x >= int16_max - 2).sum(axis=1)
# 重复样本：相邻样本完全相同
dup_adj = (np.diff(x, axis=1) == 0).sum(axis=1)
# 离群点（>5σ，按各通道自身 σ）
sigma = stats["std"]
out5 = (np.abs(x - m1[:, None]) > 5 * sigma[:, None]).sum(axis=1)
out8 = (np.abs(x - m1[:, None]) > 8 * sigma[:, None]).sum(axis=1)
# 缺失值：int16 原始数组不存在 NaN，但检查特殊哨兵值
n_nan_marker = ((x == -32768) | (x == 32767)).sum()

# 全数组重复行（通道间完全相同的行 = 死通道对）
del x
gc.collect()

# 相关性：用子采样（前 300k 样本）计算 128x128 皮尔逊相关
sub = np.asarray(mm[:, :300000], dtype=np.float64)
sub = sub - sub.mean(axis=1, keepdims=True)
sd = sub.std(axis=1, keepdims=True)
subn = sub / np.maximum(sd, 1e-9)
corr = (subn @ subn.T) / sub.shape[1]
np.save(os.path.join(FIG, "corr128.npy"), corr.astype(np.float32))

# 通道间相同行检测（前 300k 内逐样本比较会导致 O(n^2)，改为相关性=1 判断）
ident_pairs = np.argwhere(np.triu(np.abs(corr - 1.0) < 1e-12, 1))
print("完全相同通道对:", len(ident_pairs))

os.makedirs(FIG, exist_ok=True)
summary = {
    "shape": [int(nch), int(ns)],
    "dtype": str(mm.dtype),
    "duration_min": ns / 30000.0 / 60.0,
    "global": {
        "mean": float(m1.mean()), "std_mean": float(sigma.mean()),
        "std_min": float(sigma.min()), "std_max": float(sigma.max()),
        "min": float(stats["min"].min()), "max": float(stats["max"].max()),
        "q50_mean": float(stats["q50"].mean()),
        "skew_mean": float(skew.mean()), "kurt_mean": float(kurt.mean()),
    },
    "sat_total": int(sat_lo.sum() + sat_hi.sum()),
    "dup_adj_total": int(dup_adj.sum()),
    "dup_adj_pct": float(dup_adj.sum() / (nch * (ns - 1)) * 100),
    "out5_total": int(out5.sum()), "out5_pct": float(out5.sum() / (nch * ns) * 100),
    "out8_total": int(out8.sum()), "out8_pct": float(out8.sum() / (nch * ns) * 100),
    "nan_marker_total": int(n_nan_marker),
    "ident_pairs": int(len(ident_pairs)),
    "corr_offdiag_mean": float((corr.sum() - np.trace(corr)) / (nch * nch - nch)),
    "corr_offdiag_max": float(np.max(np.triu(np.abs(corr), 1))),
}
per_ch = {
    "ch": list(range(1, nch + 1)),
    "mean": m1.tolist(), "std": sigma.tolist(),
    "min": stats["min"].tolist(), "max": stats["max"].tolist(),
    "q01": stats["q01"].tolist(), "q25": stats["q25"].tolist(),
    "q50": stats["q50"].tolist(), "q75": stats["q75"].tolist(), "q99": stats["q99"].tolist(),
    "skew": skew.tolist(), "kurt": kurt.tolist(),
    "sat": (sat_lo + sat_hi).tolist(), "dup": dup_adj.tolist(),
    "out5": out5.tolist(),
}
json.dump({"summary": summary, "per_ch": per_ch},
          open(os.path.join(FIG, "eda_stats.json"), "w"), ensure_ascii=False, indent=1)
print(json.dumps(summary, ensure_ascii=False, indent=1))

# ---------- 3. 可视化 ----------
# 3.1 幅值分布直方图（全数据下采样抽样以省内存）
samp = np.asarray(mm[:, ::17], dtype=np.float64).ravel()  # ~9M 点
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
ax[0].hist(samp, bins=201, color="#3b6fb5", edgecolor="none")
ax[0].set_title("全数据幅值分布（线性坐标）"); ax[0].set_xlabel("幅值 (LSB)"); ax[0].set_ylabel("样本数")
ax[1].hist(samp, bins=201, color="#c0504d", edgecolor="none")
ax[1].set_yscale("log"); ax[1].set_title("幅值分布（对数纵轴）"); ax[1].set_xlabel("幅值 (LSB)")
fig.suptitle(f"rat128 全体样本幅值分布（n≈{samp.size/1e6:.0f}M，128 通道合并）", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/01_hist.png", dpi=130); plt.close(fig)

# 3.2 逐通道箱线图（128 通道全画）
fig, ax = plt.subplots(figsize=(14, 5))
bp = ax.boxplot([np.asarray(mm[c, ::29], dtype=np.float64) for c in range(nch)],
                showfliers=False, patch_artist=True, widths=0.7,
                medianprops=dict(color="#d62728"))
for b in bp["boxes"]:
    b.set(facecolor="#9dc3e6", alpha=0.8)
ax.set_xlabel("通道号"); ax.set_ylabel("幅值 (LSB)")
ax.set_title("128 通道幅值箱线图（每通道抽样 ~41k 点，四分位框 + 中位线，不画离群点）")
fig.tight_layout(); fig.savefig(f"{FIG}/02_box.png", dpi=130); plt.close(fig)

# 3.3 通道噪声 σ 条形图（按通道号排列）
fig, ax = plt.subplots(figsize=(14, 4))
colors = ["#c0504d" if s > 2 * np.median(sigma) else "#3b6fb5" for s in sigma]
ax.bar(range(1, nch + 1), sigma, color=colors)
med = np.median(sigma)
ax.axhline(med, color="#333", ls="--", lw=1, label=f"中位数 σ={med:.1f}")
ax.set_xlabel("通道号"); ax.set_ylabel("标准差 σ (LSB)")
ax.set_title("各通道噪声水平（标准差）分布"); ax.legend()
fig.tight_layout(); fig.savefig(f"{FIG}/03_chstd.png", dpi=130); plt.close(fig)

# 3.4 相关性热力图
fig, ax = plt.subplots(figsize=(8.6, 7.4))
im = ax.imshow(corr, cmap="RdBu_r", vmin=-0.5, vmax=0.5)
ax.set_title("128×128 通道间皮尔逊相关系数热力图（前 30 万样本）")
ax.set_xlabel("通道号"); ax.set_ylabel("通道号")
fig.colorbar(im, ax=ax, shrink=0.85, label="相关系数 r")
fig.tight_layout(); fig.savefig(f"{FIG}/04_corr.png", dpi=130); plt.close(fig)

# 3.5 代表性通道时域波形（高噪/中位/低噪各 1，200ms 窗口）
order = np.argsort(sigma)
picks = [int(order[0]), int(order[nch // 2]), int(order[-1])]
t = np.arange(6000) / 30000.0 * 1000.0
fig, axes = plt.subplots(3, 1, figsize=(12, 7), sharex=True)
for axi, c in zip(axes, picks):
    w = np.asarray(mm[c, 100000:106000], dtype=np.float64) * 0.195
    axi.plot(t, w, lw=0.6, color="#3b6fb5")
    axi.set_ylabel("µV"); axi.set_title(f"通道 {c+1}（σ={sigma[c]:.1f} LSB = {sigma[c]*0.195:.1f} µV）", fontsize=10)
axes[-1].set_xlabel("时间 (ms)")
fig.suptitle("代表性通道时域波形（200 ms，低噪 / 中位 / 高噪通道）", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/05_wave.png", dpi=130); plt.close(fig)

# 3.6 幅值分位数剖面（每通道 q01/q25/q50/q75/q99 散点）
fig, ax = plt.subplots(figsize=(12, 4.5))
xs = np.arange(1, nch + 1)
for name, col in [("q01", "#c0504d"), ("q25", "#e8a33d"), ("q50", "#333"),
                  ("q75", "#e8a33d"), ("q99", "#c0504d")]:
    ax.plot(xs, stats[name], "o", ms=3, color=col, label=name)
ax.set_xlabel("通道号"); ax.set_ylabel("幅值 (LSB)"); ax.legend(ncol=5)
ax.set_title("各通道分位数剖面（q01/q25/q50/q75/q99）")
fig.tight_layout(); fig.savefig(f"{FIG}/06_quant.png", dpi=130); plt.close(fig)

print("figures done:", os.listdir(FIG))
