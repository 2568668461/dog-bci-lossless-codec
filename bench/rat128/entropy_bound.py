# -*- coding: utf-8 -*-
"""验证：为什么本数据集无损压缩率卡在 ~2.4x —— 计算残差熵的理论上限"""
import numpy as np
import math

D = np.load('D:/bcidata/rat128_raw.npy')  # (128, NS) int16
NCH, NS = D.shape
print(f'data: {NCH}ch x {NS} samples, {D.dtype}, LSB=0.195uV')

def hist_entropy(a):
    """离散整数残差的经验熵 (bits/sample)"""
    mn, mx = a.min(), a.max()
    cnt = np.bincount(a - mn, minlength=mx - mn + 1).astype(np.float64)
    p = cnt / cnt.sum()
    p = p[p > 0]
    return float(-(p * np.log2(p)).sum())

# ---------- 1. 差分残差熵（一阶预测下界近似） ----------
d = np.diff(D.astype(np.int32), axis=1)
H_delta = np.mean([hist_entropy(d[c]) for c in range(NCH)])
sigma_delta = d.std()
print(f'\n[一阶差分]  sigma={sigma_delta:.1f} LSB ({sigma_delta*0.195:.1f} uV)')
print(f'  经验熵={H_delta:.2f} bits  ->  熵编码上限 CR = {16/H_delta:.2f}x')

# ---------- 2. 最优一阶线性预测 (x_t = a*x_{t-1}) ----------
Hs, sigs = [], []
for c in range(NCH):
    x = D[c].astype(np.float64)
    a = np.dot(x[1:], x[:-1]) / np.dot(x[:-1], x[:-1])
    r = x[1:] - a * x[:-1]
    Hs.append(hist_entropy(np.round(r).astype(np.int64)))
    sigs.append(r.std())
H_lpc1 = float(np.mean(Hs))
print(f'\n[最优一阶预测]  a={np.mean([ (np.dot(D[c].astype(np.float64)[1:],D[c].astype(np.float64)[:-1])/np.dot(D[c].astype(np.float64)[:-1],D[c].astype(np.float64)[:-1])) for c in range(5)]):.3f}..  sigma={np.mean(sigs):.1f} LSB ({np.mean(sigs)*0.195:.1f} uV)')
print(f'  经验熵={H_lpc1:.2f} bits  ->  熵编码上限 CR = {16/H_lpc1:.2f}x')

# ---------- 3. 二阶预测 x_t = 2x_{t-1} - x_{t-2} ----------
r2 = D[:, 2:].astype(np.int32) - 2*D[:, 1:-1].astype(np.int32) + D[:, :-2].astype(np.int32)
H_2nd = np.mean([hist_entropy(r2[c]) for c in range(NCH)])
print(f'\n[二阶预测]  经验熵={H_2nd:.2f} bits  ->  CR = {16/H_2nd:.2f}x')

# ---------- 4. 跨通道实验：用相邻通道预测能省多少 ----------
# 检查通道间相关系数（去自身时间趋势后）
H_nc, gain = [], []
for c in range(NCH - 1):
    x = D[c].astype(np.float64); y = D[c+1].astype(np.float64)
    b = np.dot(x, y) / np.dot(x, x)
    r = y - b * x   # 纯跨通道预测（不用时间信息）
    H_nc.append(hist_entropy(np.round(r).astype(np.int64)))
    gain.append(H_nc[-1])
H_cross = float(np.mean(gain))
rho = np.corrcoef(D.astype(np.float64))[np.triu_indices(NCH, 1)]
print(f'\n[纯跨通道预测(邻道)]  通道对相关系数均值={rho.mean():.3f} (最大={rho.max():.3f})')
print(f'  经验熵={H_cross:.2f} bits  ->  CR = {16/H_cross:.2f}x   (vs 时间预测 {H_lpc1:.2f} bits)')

# ---------- 5. 关键实验：丢掉 k 个低位后再算（文献"3.x"的来源） ----------
print('\n[丢低位后的熵与CR]  (误差界 = (2^k - 1) * 0.195 uV)')
for k in range(0, 5):
    Dq = (D.astype(np.int32) >> k)
    dq = np.diff(Dq, axis=1)
    Hq = np.mean([hist_entropy(dq[c]) for c in range(NCH)])
    err_uv = (2**k - 1) * 0.195
    print(f'  丢{k}位: 残差熵={Hq:.2f} bits  ->  CR={16/Hq:.2f}x   最大误差 {err_uv:.2f} uV (噪声~{sigma_delta*0.195/1.414:.1f} uV 的 {err_uv/(sigma_delta*0.195/1.414)*100:.0f}%)')

# ---------- 6. 模拟 12-bit 系统（spike-band ASIC 常用） ----------
# 假设信号本身只占 12 bit（即噪声以 12bit 计数表示）
print(f'\n[12-bit 等效系统]  同样物理噪声折算到 12bit: 熵~{H_lpc1-4:.2f} bits/12bit样本 -> CR={12/(H_lpc1-4):.2f}x')

# ---------- 7. 对照：我们的实测 ----------
print('\n[实测对照] encode16w=2.38x  FLAC=2.34x  drice=2.29x  JPEG-LS=2.27x')
print(f'理论最优一阶预测熵上限: {16/H_lpc1:.2f}x  <- 我们离上限差 {(16/H_lpc1-2.38)/(16/H_lpc1)*100:.1f}%')
