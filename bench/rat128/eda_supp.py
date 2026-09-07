# -*- coding: utf-8 -*-
"""补充分析：PCA 共模分量、相邻差分、相关性分布"""
import numpy as np, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
rcParams["axes.unicode_minus"] = False

D = "D:/bcidata"; FIG = f"{D}/eda_figs"
mm = np.load(f"{D}/rat128_raw.npy", mmap_mode="r")

sub = np.asarray(mm[:, :300000], dtype=np.float64)
sub = sub - sub.mean(axis=1, keepdims=True)

# PCA：特征分解相关矩阵
corr = np.load(f"{FIG}/corr128.npy").astype(np.float64)
w, v = np.linalg.eigh(corr)
w = w[::-1]
evr = w / w.sum()
cum = np.cumsum(evr)
print("PC1 explained: %.3f, PC2: %.3f, PC1-5 cum: %.3f, PC1-10 cum: %.3f" %
      (evr[0], evr[1], cum[4], cum[9]))

# 相邻差分统计（时域相关性 -> 差分熵/压缩相关）
diff = np.diff(sub, axis=1)
dst = diff.std(axis=1)
sig = sub.std(axis=1)
print("sigma mean %.1f, |diff| sigma mean %.1f, lag1 corr mean %.3f" %
      (sig.mean(), dst.mean(), ((sub[:, :-1]*sub[:, 1:]).mean(axis=1)/sig**2).mean()))

# 相关性分布（排除对角）
off = corr[np.triu_indices(128, 1)]
print("corr percentiles:", np.percentile(off, [1, 25, 50, 75, 99]).round(3))

# 高通（去共模）：减去通道逐样本中位数 = 共模估计
common = np.median(sub, axis=0)
resid = sub - common
rs = resid.std(axis=1)
print("去共模后 sigma mean %.1f (原 %.1f), 共模占比 %.1f%%" %
      (rs.mean(), sig.mean(), 100*(1-(rs**2).mean()/ (sig**2).mean())))

# 图 7：相关性分布直方图 + PCA 累计方差
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
ax[0].hist(off, bins=81, color="#3b6fb5", edgecolor="none")
ax[0].set_title("通道对相关系数分布（8128 对）"); ax[0].set_xlabel("r"); ax[0].set_ylabel("通道对数")
ax[1].plot(np.arange(1, 129), cum, color="#c0504d", lw=1.8)
ax[1].axhline(0.9, ls="--", color="#333", lw=1)
ax[1].set_xlabel("主成分序号"); ax[1].set_ylabel("累计方差占比")
ax[1].set_title("PCA 累计方差（PC1 = %.1f%%）" % (evr[0]*100))
fig.suptitle("通道间冗余结构：相关分布与主成分", fontsize=12)
fig.tight_layout(); fig.savefig(f"{FIG}/07_pca_corr.png", dpi=130); plt.close(fig)

# 图 8：共模去除前后对比（一个通道）
c = 64
fig, ax = plt.subplots(2, 1, figsize=(12, 5.4), sharex=True)
t = np.arange(6000)/30.0
ax[0].plot(t, sub[c, 100000:106000]*0.195, lw=0.6, color="#3b6fb5")
ax[0].set_title("原始通道 65（含共模分量，σ=%.1f µV）" % (sig[c]*0.195), fontsize=10)
ax[0].set_ylabel("µV")
ax[1].plot(t, resid[c, 100000:106000]*0.195, lw=0.6, color="#2e8b57")
ax[1].set_title("去除共模（通道中位数）后（σ=%.1f µV）" % (rs[c]*0.195), fontsize=10)
ax[1].set_xlabel("时间 (ms)"); ax[1].set_ylabel("µV")
fig.tight_layout(); fig.savefig(f"{FIG}/08_common.png", dpi=130); plt.close(fig)

json.dump({"pc1": float(evr[0]), "pc2": float(evr[1]), "cum5": float(cum[4]),
           "cum10": float(cum[9]), "corr_med": float(np.median(off)),
           "corr_p99": float(np.percentile(off, 99)),
           "sig_mean": float(sig.mean()), "diff_sig_mean": float(dst.mean()),
           "lag1_mean": float(((sub[:, :-1]*sub[:, 1:]).mean(axis=1)/sig**2).mean()),
           "resid_sig_mean": float(rs.mean()),
           "common_pct": float(100*(1-(rs**2).mean()/(sig**2).mean()))},
          open(f"{FIG}/eda_supp.json", "w"), indent=1)
print("done")
