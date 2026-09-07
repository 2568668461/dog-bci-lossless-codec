import numpy as np, soundfile as sf

def ent_int(r):
    q = np.rint(r).astype(np.int64)
    v, c = np.unique(q, return_counts=True)
    p = c / c.sum()
    return float(-(p * np.log2(p)).sum())

def analyze(name, x):
    x = x.astype(np.float64)
    d = np.diff(x)
    print(f"[{name}] n={len(x)} 峰值={int(np.abs(x).max())} ({int(np.abs(x).max()).bit_length()}bit)")
    print(f"  差分残差 sigma={d.std():.1f} LSB   差分熵={ent_int(d):.2f} bits -> CR上限={16/ent_int(d):.2f}x")
    # AR(1)
    X = x[:-1]; y = x[1:]
    a = np.dot(y, X)/np.dot(X, X)
    r = y - a*X
    print(f"  AR(1) 熵={ent_int(r):.2f} bits -> CR上限={16/ent_int(r):.2f}x  (a={a:.4f})")
    # 频谱：各频带噪声功率（预测残差近似白噪声，看每 Hz 密度）
    f = np.fft.rfftfreq(len(d), 1/30000)
    P = np.abs(np.fft.rfft(d))**2
    bands = [(0,300),(300,1000),(1000,5000),(5000,10000),(10000,15000)]
    tot = P.mean()
    print("  残差功率谱占比:", "  ".join(f"{lo}-{lo2}Hz:{P[(f>=lo)&(f<lo2)].mean()/tot*100:.0f}%" for lo,lo2 in bands))
    print()

print("===== 挑战赛风格数据（archive testdata, 30kHz, 峰值13bit）=====")
for ch in ['ch000','ch063','ch127']:
    x, sr = sf.read(f'D:/bcidata/neuralink-codec-archive/testdata/{ch}.wav', dtype='int16')
    analyze(ch, x)

print("===== 我们的真实 Intan 大鼠数据（rat128, 20kHz）=====")
D = np.load('D:/bcidata/rat128_raw.npy')
for c in [0, 63, 127]:
    analyze(f'rat_ch{c:03d}', D[c])
