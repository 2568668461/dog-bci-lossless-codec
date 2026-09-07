# -*- coding: utf-8 -*-
"""bench3b —— 分解分析: 跨通道收益多少来自重复通道冗余, 多少来自真实相关?
1) T1d 端到端可逆性校验
2) 分脑区: FLAC单声道 vs FLAC立体声配对 vs T1残差FLAC 的 CR 对比
"""
import os
import numpy as np
import imagecodecs
import soundfile as sf

D = np.load('D:/bcidata/dandi001539/ob48_lfp.npy')
NCH, NS = D.shape
FS = 1500.0
OUT = 'D:/bcidata/dandi001539/flac_out3'

X = D.astype(np.float64)
Xn = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
C = Xn @ Xn.T / NS
import h5py
_f = h5py.File('D:/bcidata/dandi001539/sub-CS39_ses-06.nwb', 'r')
_el = _f['general/extracellular_ephys/electrodes']
_loc = np.array([x.decode() if isinstance(x, bytes) else str(x)
                 for x in np.atleast_1d(_el['location'][()])])
_f.close()
REG = [(r, [int(i) for i in np.where(_loc == r)[0]]) for r in ['OB', 'CA1', 'PFC']]
print('脑区映射 (NWB electrodes):', REG)

# ---- 1. T1d 可逆性 ----
# 重建 T1 残差 (与 bench3 相同的 Prim 树)
def build_tree(chs):
    order = [chs[0]]; parent = {}; resolved = set([chs[0]])
    while len(resolved) < len(chs):
        best = None
        for i in chs:
            if i in resolved: continue
            for j in resolved:
                if best is None or C[i, j] > best[0]: best = (C[i, j], i, j)
        _, i, j = best
        parent[i] = j; resolved.add(i); order.append(i)
    return order, parent

topo, parent, roots = [], {}, []
for rn, chs in REG:
    o, p = build_tree(chs); roots.append(o[0]); topo += o; parent.update(p)
resid = np.zeros((NCH, NS), dtype=np.int32)
for i, j in parent.items(): resid[i] = D[i].astype(np.int32) - D[j].astype(np.int32)
chmap = [i for i in topo if i not in roots]
enc = sorted(chmap)
r16 = np.ascontiguousarray(resid[enc]).astype(np.int16)
rd = np.empty_like(r16)
rd[:, 0] = r16[:, 0]
rd[:, 1:] = (r16[:, 1:].astype(np.int32) - r16[:, :-1].astype(np.int32)).astype(np.int16)
c = imagecodecs.zstd_encode(rd, level=19)
back = np.frombuffer(imagecodecs.zstd_decode(c), dtype=np.int16).reshape(len(enc), NS)
x = np.cumsum(back.astype(np.int64), axis=1)
r_rec = ((x + 32768) % 65536 - 32768).astype(np.int16)
print('T1d zstd 解压+反差分 == 原残差:', np.array_equal(r_rec, r16))

# ---- 2. 分脑区 CR ----
def flac_bytes(arr2d, tag):
    total = 0
    for k in range(arr2d.shape[0]):
        fn = '%s/%s_%02d.flac' % (OUT, tag, k)
        sf.write(fn, arr2d[k], int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
    return total

def flac_stereo_bytes(pairs, tag):
    total = 0; ok = True
    for a, b in pairs:
        fn = '%s/%s_%d_%d.flac' % (OUT, tag, a, b)
        sf.write(fn, np.stack([D[a], D[b]]).T, int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
        bb, _ = sf.read(fn, dtype='int16')
        ok = ok and np.array_equal(bb.T, np.stack([D[a], D[b]]))
    return total, ok

print('\n%-6s %10s %12s %12s %12s' % ('脑区', 'FLAC单声道', 'FLAC立体声配对', 'T1残差FLAC', '最高相关'))
for rn, chs in REG:
    sub = D[chs]
    n_mono = flac_bytes(sub, 'm_' + rn)
    # 配对
    pool = list(chs); pairs = []
    while len(pool) >= 2:
        best = None
        for x1 in range(len(pool)):
            for x2 in range(x1 + 1, len(pool)):
                cc = C[pool[x1], pool[x2]]
                if best is None or cc > best[0]: best = (cc, pool[x1], pool[x2])
        _, a, b = best; pairs.append((a, b)); pool.remove(a); pool.remove(b)
    n_st, ok = flac_stereo_bytes(pairs, 's_' + rn)
    # T1 残差 (本脑区内树)
    r_ch = [i for i in chs if i in parent]
    resid_sub = np.ascontiguousarray(resid[r_ch]).astype(np.int16)
    n_res = flac_bytes(resid_sub, 'r_' + rn) + 0  # 根通道成本已在单声道里, 近似用根FLAC单通道
    root = [i for i in chs if i not in parent][0]
    n_root = flac_bytes(np.ascontiguousarray(D[root][None, :]), 'rt_' + rn)
    n_res_total = n_res + n_root
    raw = sub.nbytes
    cmax = C[np.ix_(chs, chs)][np.triu_indices(len(chs), 1)].max()
    print('%-6s %8.3fx %12.3fx %14.3fx %12.3f' % (rn, raw / n_mono, raw / n_st, raw / n_res_total, cmax))

# ---- 3. 去重后的全局 FLAC 立体声 (剔除硬件冗余后的真实增益) ----
uniq = []
for rn, chs in REG:
    pool = list(chs)
    while len(pool) >= 2:
        best = None
        for x1 in range(len(pool)):
            for x2 in range(x1 + 1, len(pool)):
                if C[pool[x1], pool[x2]] > 0.9999:
                    best = (2.0, pool[x1], pool[x2]); break
            if best: break
        if best:
            _, a, b = best  # 保留一个, 丢弃重复的
            pool.remove(b)
        else:
            break
    uniq += pool
print('\n去重后通道数: %d (原 48)' % len(uniq))
sub = D[uniq]
n_mono = flac_bytes(sub, 'um')
raw = sub.nbytes
pool = list(uniq); pairs = []
while len(pool) >= 2:
    best = None
    for x1 in range(len(pool)):
        for x2 in range(x1 + 1, len(pool)):
            cc = C[pool[x1], pool[x2]]
            if best is None or cc > best[0]: best = (cc, pool[x1], pool[x2])
    _, a, b = best; pairs.append((a, b)); pool.remove(a); pool.remove(b)
if pool:
    fn = '%s/u_single_%d.flac' % (OUT, pool[0])
    sf.write(fn, D[pool[0]], int(FS), format='FLAC', subtype='PCM_16')
    n_st = flac_stereo_bytes(pairs, 'us')[0] + os.path.getsize(fn)
else:
    n_st = flac_stereo_bytes(pairs, 'us')[0]
print('去重后: FLAC单声道 %.3fx, FLAC立体声配对 %.3fx' % (raw / n_mono, raw / n_st))
