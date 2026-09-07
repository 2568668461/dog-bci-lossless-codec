# -*- coding: utf-8 -*-
"""近无损策略实测 v2
P0 族(锚点): 全局 NEAR = 0/2/4/8  -> 误差上界 = NEAR LSB 恒定
P1 族:       每通道 NEAR_c ~ sigma_c / {4, 2} (噪声比例, 误差随通道噪声缩放)
P2 族:       spike 感知分块: 活动块 NEAR=1(保 spike, ~0.2uV) + 背景块 NEAR~sigma/4 或 sigma/2
数据: rat128_raw.npy (128ch x 1.2M, int16, 0.195uV/LSB)
"""
import numpy as np, imagecodecs, time, json

t0 = time.time()
D = np.load('D:/bcidata/rat128_raw.npy')
NCH, NS = D.shape
L = 3000
NB = NS // L
print(f'data {NCH}x{NS}, blocks {NB}x{L}', flush=True)

dall = np.abs(np.diff(D.astype(np.int32), axis=1))
mad_d = np.median(dall, axis=1)
sig_c = 1.4826 * mad_d / np.sqrt(2)
thr_c = 4.0 * mad_d

def u16(ch):
    return (ch.astype(np.int32) + 32768).astype(np.uint16)

rows = []
def report(name, nbytes, note=''):
    cr = D.nbytes / nbytes
    rows.append({'policy': name, 'CR': round(cr, 3), 'Mbps': round(61.44 / cr, 2),
                 'bits/sample': round(16 / cr, 3), 'note': note})
    print(f'{name:34s} CR={cr:6.3f}x {61.44/cr:6.2f} Mbps bits={16/cr:.3f}  {note}', flush=True)

def enc_level_all(level):
    n = 0
    for c in range(NCH):
        n += len(imagecodecs.jpegls_encode(u16(np.ascontiguousarray(D[c])).reshape(NB, L), level=level))
    return n

# ---------- 锚点：全局固定 NEAR ----------
for lv in [0, 1, 2, 4, 8]:
    report(f'P0 全局NEAR={lv}', enc_level_all(lv), f'err<={lv}LSB')

# ---------- P1：噪声比例每通道 ----------
for denom, tag in [(4, 'sig/4'), (2, 'sig/2')]:
    ne = np.clip(np.round(sig_c / denom), 1, 16).astype(int)
    n = 0
    for c in range(NCH):
        n += len(imagecodecs.jpegls_encode(u16(np.ascontiguousarray(D[c])).reshape(NB, L), level=int(ne[c])))
    report(f'P1 噪声自适应NEAR~{tag}', n, f'NEAR_c {ne.min()}/{int(np.median(ne))}/{ne.max()}')

# ---------- P2：spike 感知分块（阈值灵敏度扫描） ----------
ne4 = np.clip(np.round(sig_c / 4), 1, 16).astype(int)   # 背景: err<=sig/4
ne2 = np.clip(np.round(sig_c / 2), 1, 16).astype(int)   # 背景: err<=sig/2

def p2_run(neq, f):
    n = 0; act_tot = 0
    for c in range(NCH):
        r = D[c][:NB * L].reshape(NB, L)
        md = np.abs(np.diff(r.astype(np.int32), axis=1)).max(axis=1)
        act = md > f * mad_d[c]
        act_tot += int(act.sum())
        u = u16(np.ascontiguousarray(D[c])).reshape(NB, L)
        for bi in range(NB):
            lv = 1 if act[bi] else int(neq[c])
            n += len(imagecodecs.jpegls_encode(u[bi:bi + 1], level=lv))
    return n, act_tot / (NCH * NB)

for f in [5, 6, 7, 8]:
    n, frac = p2_run(ne2, f)
    print(f'[thr={f}xMAD] active blocks {frac*100:.2f}%  (activity scan)', flush=True)

for tag, neq in [('sig/4', ne4), ('sig/2', ne2)]:
    n, frac = p2_run(neq, 7.0)
    report(f'P2 spike感知(thr=7xMAD,活动NEAR=1,背景~{tag})', n, f'活动块 {frac*100:.2f}%')

print('elapsed %.0fs' % (time.time() - t0), flush=True)

# ---------- 校验：抽验误差上界 ----------
rng = np.random.default_rng(1)
def maxerr(u2d, lv):
    r = imagecodecs.jpegls_decode(imagecodecs.jpegls_encode(u2d, level=int(lv)))
    if isinstance(r, tuple):
        r = r[0]
    return int(np.abs(r.astype(np.int32) - u2d.astype(np.int32)).max())

worst = 0
for c in rng.choice(NCH, 4, replace=False):
    u = u16(np.ascontiguousarray(D[c])).reshape(NB, L)
    for lv in [0, 2, 8, int(ne2[c]), int(ne4[c])]:
        m = maxerr(u, lv); worst = max(worst, m - lv)
print('verify uniform/adaptive passes OK, worst exceed =', worst, flush=True)
# P2 抽 80 块
for _ in range(80):
    c = int(rng.integers(NCH)); bi = int(rng.integers(NB))
    r = D[c][:NB * L].reshape(NB, L)
    md = np.abs(np.diff(r.astype(np.int32), axis=1)).max(axis=1)
    act = md > 7.0 * mad_d[c]
    lv = 1 if act[bi] else int(ne2[c])
    u = u16(np.ascontiguousarray(D[c])).reshape(NB, L)
    m = maxerr(u[bi:bi + 1], lv); worst = max(worst, m - lv)
print('verify P2 blocks OK, worst exceed =', worst, flush=True)

np.savez('D:/bcidata/nearloss_stats.npz', sig_c=sig_c, ne4=ne4, ne2=ne2, thr_c=thr_c)
json.dump(rows, open('D:/bcidata/nearloss_rows.json', 'w'), ensure_ascii=False, indent=1)
print('DONE', flush=True)
