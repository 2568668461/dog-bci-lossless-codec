# -*- coding: utf-8 -*-
"""001539 专项: 重复通道 / 分脑区压缩 / 跨通道预测"""
import numpy as np, soundfile as sf, os, h5py, imagecodecs

os.chdir('D:/bcidata/dandi001539')
D = np.load('ob48_lfp.npy')
f = h5py.File('sub-CS39_ses-06.nwb', 'r')
el = f['general/extracellular_ephys/electrodes']
loc = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in np.atleast_1d(el['location'][()])])

print('== 重复通道检测 ==')
dup = []
for i in range(48):
    for j in range(i + 1, 48):
        if np.array_equal(D[i], D[j]):
            dup.append((i, j, loc[i]))
print('完全相同通道对:', dup if dup else '无')

print('\n== 分脑区压缩 ==')
L = 1024
NB = D.shape[1] // L
def u16(a):
    return (a.astype(np.int32) + 32768).astype(np.uint16)
for name in ['OB', 'CA1', 'PFC']:
    idx = np.where(loc == name)[0]
    sub = D[idx]
    total = 0
    for ch_i in idx:
        fn = 'flac_out/rg%d.flac' % ch_i
        sf.write(fn, D[ch_i], 1500, format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
    cr_flac = sub.nbytes / total
    n0 = n2 = 0
    for ch_i in idx:
        b = u16(np.ascontiguousarray(D[ch_i, :NB * L])).reshape(NB, L)
        n0 += len(imagecodecs.jpegls_encode(b, level=0))
        n2 += len(imagecodecs.jpegls_encode(b, level=2))
    ns = sub.shape[1]
    print('%s (%dch): sigma中位 %4.0f uV | FLAC CR %.3f | JLS-0 CR %.3f | JLS-2 CR %.3f' %
          (name, len(idx), np.median(sub.astype(float).std(axis=1)), cr_flac,
           sub.nbytes / n0 * ns / (NB * L), sub.nbytes / n2 * ns / (NB * L)))

total = 0
for ch in range(48):
    fn = 'flac_out/l8_%d.flac' % ch
    sf.write(fn, D[ch], 1500, format='FLAC', subtype='PCM_16', compression_level=1.0)
    total += os.path.getsize(fn)
print('\nFLAC 最高压缩级: CR %.3f' % (D.nbytes / total))

corr = np.corrcoef(D.astype(np.float64))
ob = np.where(loc == 'OB')[0]
co = corr[np.ix_(ob, ob)]
pairs = [(ob[i], ob[j], co[i, j]) for i in range(len(ob)) for j in range(i + 1, len(ob)) if co[i, j] > 0.95]
print('OB 内 r>0.95 通道对:', [(int(a), int(b), round(r, 4)) for a, b, r in pairs])

if pairs:
    a, b, r0 = pairs[0]
    xa = D[a].astype(np.float64)
    xb = D[b].astype(np.float64)
    k = np.dot(xa, xb) / np.dot(xa, xa)
    res = np.round(xb - k * xa)
    res = np.clip(res, -32768, 32767).astype(np.int16)
    sf.write('flac_out/res.flac', res, 1500, format='FLAC', subtype='PCM_16')
    sz_res = os.path.getsize('flac_out/res.flac')
    sf.write('flac_out/orig.flac', D[b], 1500, format='FLAC', subtype='PCM_16')
    sz_orig = os.path.getsize('flac_out/orig.flac')
    print('\n跨通道预测: ch%d=a, ch%d=b, r=%.4f, k=%.3f' % (a, b, r0, k))
    print('  ch%d 单独 FLAC %.1f KB; 残差 FLAC %.1f KB (省 %.1f%%)' %
          (b, sz_orig / 1e3, sz_res / 1e3, 100 * (1 - sz_res / sz_orig)))
    err = np.abs((xb - k * xa) - res).max()
    print('  量化引入 worst err: %.2f uV (LSB=1uV)' % err)
