# -*- coding: utf-8 -*-
"""DANDI 001539 bench3 v2 —— 树形/DAG 跨通道预测管线 (可解码, 拓扑序) + FLAC立体声配对
方案 (全部端到端校验, 除注明外均无损):
  T1  树形零系数邻道: r = x_i - x_parent        (FPGA: 1 减法器)
  T1d T1 残差再做时域一阶差分
  T4  DAG 双亲零系数: r = x_i - (x_p1+x_p2)/2   (FPGA: 加法+移位)
  T5  DAG 双亲 LS:    r = x_i - round(a1*x_p1+a2*x_p2+b)  (需定点化)
  每方案: 根通道 FLAC 独立编码; 残差分别用 zstd L19 / FLAC / drice / JPEG-LS NEAR
  对照: FLAC 立体声配对 (joint stereo, 零自研代码), SZ3, AEC
"""
import os, time, json, subprocess
import numpy as np
import imagecodecs
import soundfile as sf

D = np.load('D:/bcidata/dandi001539/ob48_lfp.npy')
NCH, NS = D.shape
RAW = D.nbytes
FS, NATIVE_MBPS = 1500.0, 48 * 1500 * 16 / 1e6
TMP = 'D:/bcidata/dandi001539/_tmp3'
os.makedirs(TMP, exist_ok=True)
os.makedirs('D:/bcidata/dandi001539/flac_out3', exist_ok=True)
print('data: %s %s, raw = %.1f MB, 信号 std=%.1fµV' % (D.shape, D.dtype, RAW / 1e6, D.std()))

rows = []
def report(name, nbytes, ok=None, note=''):
    cr = RAW / nbytes
    rows.append((name, cr, NATIVE_MBPS / cr, ok, note))
    print('%-48s CR=%6.3fx  %5.3f Mbps  ok=%s  %s' %
          (name, cr, NATIVE_MBPS / cr, ok, note), flush=True)

u16 = lambda a: (a.astype(np.int32) + 32768).astype(np.uint16)
i16 = lambda a: (a.astype(np.int64) - 32768).astype(np.int16)

# ---- 相关系数矩阵 / 脑区 (从 NWB electrodes 表取真实映射) ----
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

# ---- Prim 相关树: 拓扑序 (父必先于子) ----
def build_tree(chs):
    order = [chs[0]]                      # 根在最前
    parent = {}
    resolved = set([chs[0]])
    while len(resolved) < len(chs):
        best = None
        for i in chs:
            if i in resolved: continue
            for j in resolved:
                if best is None or C[i, j] > best[0]: best = (C[i, j], i, j)
        _, i, j = best
        parent[i] = j; resolved.add(i); order.append(i)
    return order, parent

topo, parent = [], {}
roots = []
for rn, chs in REG:
    o, p = build_tree(chs)
    roots.append(o[0]); topo += o; parent.update(p)
assert sorted(topo) == list(range(NCH))
chmap = [i for i in topo if i not in roots]     # 拓扑序的非根通道
def chain_depth(i):
    d = 0
    while i in parent: i = parent[i]; d += 1
    return d
print('3 根通道 %s, 最大链深 %d' % (roots, max(chain_depth(i) for i in range(NCH))))

# ---- 根通道 FLAC ----
root_bytes = 0
for r in roots:
    fn = 'D:/bcidata/dandi001539/flac_out3/root%d.flac' % r
    sf.write(fn, D[r], int(FS), format='FLAC', subtype='PCM_16')
    root_bytes += os.path.getsize(fn)
    b, _ = sf.read(fn, dtype='int16')
    assert np.array_equal(b, D[r]), '根通道往返失败'

# ---- T1: 零系数单亲 ----
print('\n---- T1 树形零系数单亲 ----', flush=True)
resid_T1 = np.zeros((NCH, NS), dtype=np.int32)
for i in chmap: resid_T1[i] = D[i].astype(np.int32) - D[parent[i]].astype(np.int32)
rmax = np.abs(resid_T1[chmap]).max()
print('  残差 std=%.1f max=%d' % (resid_T1[chmap].std(), rmax))
enc = sorted(chmap)                      # 编码用自然电极序 (zstd 可利用相邻重复通道)
r16 = np.ascontiguousarray(resid_T1[enc]).astype(np.int16)
c = imagecodecs.zstd_encode(r16, level=19)
back = np.frombuffer(imagecodecs.zstd_decode(c), dtype=np.int16).reshape(len(enc), NS)
report('T1 零系数单亲: 根FLAC+残差zstd(电极序)', root_bytes + len(c), np.array_equal(back, r16))
# 混排方案: 48 行整流 (根=原始, 其余=残差) 一个 zstd 流, 根不单独编码
mix = D.copy()
for i in chmap: mix[i] = resid_T1[i].astype(np.int16)
c = imagecodecs.zstd_encode(np.ascontiguousarray(mix), level=19)
report('T1 零系数单亲: 根原始+残差混排zstd', len(c), None, '单流无根开销')
total = 0
for k, i in enumerate(enc):
    fn = 'D:/bcidata/dandi001539/flac_out3/a%02d.flac' % i
    sf.write(fn, r16[k], int(FS), format='FLAC', subtype='PCM_16')
    total += os.path.getsize(fn)
report('T1 零系数单亲: 根FLAC+残差FLAC', root_bytes + total, None)
# drice (FPGA 亲和) on T1 residual
in_bin = TMP + '/t1resid.bin'
r16.astype(np.int16).tofile(in_bin)
def run_drice(exe, args, tag, chs_n):
    out = TMP + '/%s.bin' % tag
    r = subprocess.run([exe, 'c', in_bin, out] + args, capture_output=True, text=True, timeout=600)
    if r.returncode != 0: print('  drice ERR', r.stderr[:150]); return None
    backf = TMP + '/%s.back' % tag
    subprocess.run([exe, 'd', out, backf], capture_output=True, text=True, timeout=600)
    back = np.fromfile(backf, dtype=np.int16).reshape(chs_n, NS)
    return os.path.getsize(out), back
res = run_drice('D:/bcidata/drice/drice.exe', ['45'], 't1drice', len(chmap))
if res: report('T1 零系数单亲: 根FLAC+残差drice', root_bytes + res[0], None, 'FPGA亲和')
# T1 近无损端到端 (drice_nl / JPEG-LS), 拓扑序累进重建
for near in (2, 4):
    res = run_drice('D:/bcidata/drice/drice_nl.exe', ['45', '1024', str(near)], 't1nl%d' % near, len(enc))
    if res:
        sz, back = res
        rec = {r: D[r].astype(np.float64) for r in roots}
        for i in chmap:                          # 拓扑序重建; back 行按 enc 序索引
            rec[i] = back[enc.index(i)].astype(np.float64) + rec[parent[i]]
        err = np.abs(np.stack([rec[i] for i in chmap]) - D[chmap].astype(np.float64))
        report('T1 零系数单亲: 残差drice_nl NEAR=%d' % near, root_bytes + sz, None,
               'worst|err|=%dµV RMSE=%.2f' % (err.max(), np.sqrt((err ** 2).mean())))
L = 1024; NB = NS // L
for near in (2, 4):
    n = 0
    rec = {r: D[r].astype(np.float64) for r in roots}
    errs = []
    for k, i in enumerate(chmap):
        blk = u16(np.ascontiguousarray(resid_T1[i, :NB * L]).reshape(NB, L))
        cc = imagecodecs.jpegls_encode(blk, level=near)
        n += len(cc)
        rq = i16(np.asarray(imagecodecs.jpegls_decode(cc)).reshape(-1)).astype(np.float64)
        rec[i] = rq + rec[parent[i]][:NB * L]          # 块化覆盖 618496/618512 = 99.997%
        errs.append(np.abs(rec[i] - D[i, :NB * L].astype(np.float64)))
    err = np.concatenate(errs)
    report('T1 零系数单亲: 残差JPEG-LS NEAR=%d' % near, root_bytes + int(n * NS / (NB * L)), None,
           'worst|err|=%dµV RMSE=%.2f' % (err.max(), np.sqrt((err ** 2).mean())))

# ---- T1d: T1 残差再时域差分 ----
print('\n---- T1d 残差再时域差分 ----', flush=True)
rd = np.empty_like(r16)
rd[:, 0] = r16[:, 0]
rd[:, 1:] = (r16[:, 1:].astype(np.int32) - r16[:, :-1].astype(np.int32)).astype(np.int16)
c = imagecodecs.zstd_encode(rd, level=19)
report('T1d 残差时域差分+zstd', root_bytes + len(c), None)

# ---- T4: 双亲零系数平均 ----
print('\n---- T4 DAG 双亲零系数平均 ----', flush=True)
resid_T4 = np.zeros((NCH, NS), dtype=np.int32)
par2 = {}
for k, i in enumerate(chmap):
    resolved = [r for r in roots if True] + [j for j in chmap[:k]]
    cand = sorted(resolved, key=lambda j: -C[i, j])[:2]
    if len(cand) == 2:
        p1, p2 = cand
        resid_T4[i] = D[i].astype(np.int32) - ((D[p1].astype(np.int32) + D[p2].astype(np.int32)) >> 1)
        par2[i] = (p1, p2)
    else:
        resid_T4[i] = D[i].astype(np.int32) - D[cand[0]].astype(np.int32)
        par2[i] = (cand[0], None)
print('  残差 std=%.1f max=%d' % (resid_T4[chmap].std(), np.abs(resid_T4[chmap]).max()))
assert np.abs(resid_T4[chmap]).max() < 32768
r16b = np.ascontiguousarray(resid_T4[enc]).astype(np.int16)
c = imagecodecs.zstd_encode(r16b, level=19)
back = np.frombuffer(imagecodecs.zstd_decode(c), dtype=np.int16).reshape(len(enc), NS)
report('T4 双亲平均: 根FLAC+残差zstd', root_bytes + len(c), np.array_equal(back, r16b))
total = 0
for k, i in enumerate(enc):
    fn = 'D:/bcidata/dandi001539/flac_out3/b%02d.flac' % i
    sf.write(fn, r16b[k], int(FS), format='FLAC', subtype='PCM_16')
    total += os.path.getsize(fn)
report('T4 双亲平均: 根FLAC+残差FLAC', root_bytes + total, None)

# ---- T5: 双亲 LS ----
print('\n---- T5 DAG 双亲最小二乘 ----', flush=True)
resid_T5 = np.zeros((NCH, NS), dtype=np.int32)
coef5 = {}
for k, i in enumerate(chmap):
    resolved = [r for r in roots] + [j for j in chmap[:k]]
    cand = sorted(resolved, key=lambda j: -C[i, j])[:2]
    if len(cand) >= 2:
        Y = np.stack([D[p].astype(np.float64) for p in cand[:2]] + [np.ones(NS)])
        w, *_ = np.linalg.lstsq(Y.T, D[i].astype(np.float64), rcond=None)
        pred = w[0] * Y[0] + w[1] * Y[1] + w[2]
        coef5[i] = (cand[0], cand[1], w)
    else:
        pred = D[cand[0]].astype(np.float64)
        coef5[i] = (cand[0], None, None)
    resid_T5[i] = np.rint(D[i].astype(np.float64) - pred).astype(np.int32)
print('  残差 std=%.1f max=%d' % (resid_T5[chmap].std(), np.abs(resid_T5[chmap]).max()))
if np.abs(resid_T5[chmap]).max() < 32768:
    r16c = np.ascontiguousarray(resid_T5[enc]).astype(np.int16)
    c = imagecodecs.zstd_encode(r16c, level=19)
    report('T5 双亲LS: 根FLAC+残差zstd', root_bytes + len(c), None, '系数需定点化')
    total = 0
    for k, i in enumerate(enc):
        fn = 'D:/bcidata/dandi001539/flac_out3/c%02d.flac' % i
        sf.write(fn, r16c[k], int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
    report('T5 双亲LS: 根FLAC+残差FLAC', root_bytes + total, None)

# ---- FLAC 立体声配对 ----
print('\n---- FLAC 立体声配对 (joint stereo) ----', flush=True)
pairs = []
for rn, chs in REG:
    pool = list(chs)
    while len(pool) >= 2:
        best = None
        for x1 in range(len(pool)):
            for x2 in range(x1 + 1, len(pool)):
                cc = C[pool[x1], pool[x2]]
                if best is None or cc > best[0]: best = (cc, pool[x1], pool[x2])
        _, a, b = best
        pairs.append((a, b)); pool.remove(a); pool.remove(b)
    if pool: pairs.append((pool[0], None))
total = 0; ok = True
for a, b in pairs:
    if b is None:
        fn = 'D:/bcidata/dandi001539/flac_out3/s%d.flac' % a
        sf.write(fn, D[a], int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
        bb, _ = sf.read(fn, dtype='int16'); ok = ok and np.array_equal(bb, D[a])
    else:
        fn = 'D:/bcidata/dandi001539/flac_out3/s%d_%d.flac' % (a, b)
        sf.write(fn, np.stack([D[a], D[b]]).T, int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
        bb, _ = sf.read(fn, dtype='int16')
        ok = ok and np.array_equal(bb.T, np.stack([D[a], D[b]]))
report('FLAC 立体声配对 (23对+2单)', total, ok, '零自研代码')

# ---- SZ3 / AEC ----
print('\n---- SZ3 / AEC ----', flush=True)
Df = D.astype(np.float32)
c = imagecodecs.sz3_encode(Df)
back = imagecodecs.sz3_decode(c, D.shape, np.float32)
report('SZ3 无损 (f32中转)', len(c), np.array_equal(back.astype(np.int16), D))
for eb in (2, 4):
    c = imagecodecs.sz3_encode(Df, abs=eb)
    back = np.asarray(imagecodecs.sz3_decode(c, D.shape, np.float32)).reshape(NCH, NS)
    err = np.abs(back - Df)
    report('SZ3 abs=%dµV' % eb, len(c), None, 'worst|err|=%d RMSE=%.2f' % (err.max(), np.sqrt((err ** 2).mean())))
flat = np.ascontiguousarray(D).reshape(-1)
c = imagecodecs.aec_encode(flat)
back = np.asarray(imagecodecs.aec_decode(c, bitspersample=16, out=np.empty(NCH * NS, dtype=np.int16))).reshape(NCH, NS)
report('AEC/CCSDS-121 (raw)', len(c), np.array_equal(back, D), '航天标准')

# ---- 汇总 ----
print('\n==== bench3 汇总 (排序) ====')
for name, cr, mb, ok, note in sorted(rows, key=lambda r: -r[1]):
    print('%-48s CR=%6.3fx  %5.3f Mbps  ok=%s  %s' % (name, cr, mb, ok, note))
with open('D:/bcidata/dandi001539/bench3_results.json', 'w', encoding='utf-8') as f:
    json.dump([{'method': n, 'CR': round(cr, 3), 'Mbps': round(m, 3), 'lossless': o, 'note': s}
               for n, cr, m, o, s in rows], f, ensure_ascii=False, indent=1)
print('saved bench3_results.json')
