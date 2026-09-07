import numpy as np, glob, soundfile as sf

print("=== 一、Neuralink 挑战赛数据检验：低6bit是否恒为0 ===")
for f in sorted(glob.glob('C:/Users/MI/WorkBuddy/2026-09-04-20-43-57/neuralink-codec/testdata/ch*.wav'))[:6]:
    x, sr = sf.read(f, dtype='int16')
    x = x.astype(np.int64)
    n64 = np.sum(x % 64 != 0)
    n32 = np.sum(x % 32 != 0)
    # 有效位数
    bits = int(np.max(np.abs(x))).bit_length()
    print(f"{f.split('/')[-1]}: sr={sr} n={len(x)} 非64倍数样本={n64} 非32倍数={n32} 峰值位数={bits}bit")

print()
print("=== 二、无损剩余余量测试（rat128 真实数据） ===")
D = np.load('D:/bcidata/rat128_raw.npy').astype(np.float64)  # (128, NS)
NCH, NS = D.shape

def ent_int(r):
    """量化到整数的经验熵（bit/样本）"""
    q = np.rint(r).astype(np.int64)
    v, c = np.unique(q, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())

# 1) LS 最优 p 阶线性预测器（每通道独立）
def resid_ar(x, p):
    """最小二乘 AR(p) 残差"""
    X = np.lib.stride_tricks.sliding_window_view(x, p)[:-1]  # (n-p, p)
    y = x[p:]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ coef

for p in [1, 2, 4, 8]:
    hs = [ent_int(resid_ar(D[c], p)) for c in range(NCH)]
    h = np.mean(hs)
    print(f"AR({p}) LS最优预测:  H={h:.2f} bits -> CR上限={16/h:.3f}x")

# 2) 时空联合预测: x_t ~ x_{t-1} + 邻道_t + 邻道_{t-1}
hs = []
for c in range(NCH):
    x = D[c]; nb = D[(c+1) % NCH]
    X = np.column_stack([x[:-1], nb[1:], nb[:-1]])
    y = x[1:]
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    hs.append(ent_int(y - X @ coef))
h = np.mean(hs)
print(f"时空联合预测(自滞后1+邻道当前+邻道滞后1): H={h:.2f} bits -> CR上限={16/h:.3f}x")

# 3) 上下文建模余量: 按前一残差绝对值分桶后的条件熵（AR(1)残差）
hs_c = []
for c in range(NCH):
    r = resid_ar(D[c], 1)
    ctx = np.digitize(np.abs(r[:-1]), [2, 8, 32])  # 4个上下文桶
    total = 0; n = len(r) - 1
    for k in range(4):
        m = ctx == k
        if m.sum() == 0: continue
        total += m.sum() / n * ent_int(r[1:][m])
    hs_c.append(total)
h = np.mean(hs_c)
print(f"AR(1)+4上下文条件熵: H={h:.2f} bits -> CR上限={16/h:.3f}x")

# 4) AR(8)+上下文 双重叠加
hs_c = []
for c in range(NCH):
    r = resid_ar(D[c], 8)
    ctx = np.digitize(np.abs(r[:-1]), [2, 8, 32])
    total = 0; n = len(r) - 1
    for k in range(4):
        m = ctx == k
        if m.sum() == 0: continue
        total += m.sum() / n * ent_int(r[1:][m])
    hs_c.append(total)
h = np.mean(hs_c)
print(f"AR(8)+4上下文条件熵: H={h:.2f} bits -> CR上限={16/h:.3f}x")
