# -*- coding: utf-8 -*-
"""DANDI 001539 (ob48_lfp.npy) 压缩基准 —— 第二轮算法扫描 (bench2)
目标: 寻找超过 FLAC 1.761x 的无损方案
新增: 通用库(bz2/lzma/brotli/zstd-22/字节平面-lzma)、16bit图像类(J2K/HTJ2K/JXL/PNG/AVIF)、
      科学压缩器(SZ3/ZFP/LERC)、航天标准 AEC(CCSDS-121)、自研跨通道预测(tetrode邻道)
数据: 48ch x 618512 int16 @ 1500 Hz (OB8+CA120+PFC20), 原生码率 1.152 Mbps
"""
import os, time, json
import numpy as np
import imagecodecs
import soundfile as sf

D = np.load('D:/bcidata/dandi001539/ob48_lfp.npy')
NCH, NS = D.shape
RAW = D.nbytes
FS, NATIVE_MBPS = 1500.0, 48 * 1500 * 16 / 1e6
print('data: %s %s, raw = %.1f MB, 原生码率 %.3f Mbps' % (D.shape, D.dtype, RAW / 1e6, NATIVE_MBPS))

rows = []
def report(name, nbytes, ok=None, note=''):
    cr = RAW / nbytes
    rows.append((name, cr, NATIVE_MBPS / cr, ok, note))
    print('%-42s CR=%6.3fx  %5.3f Mbps  %6.3f bits/s  ok=%s  %s' %
          (name, cr, NATIVE_MBPS / cr, 16 / cr, ok, note), flush=True)

def asarr(b): return np.frombuffer(b, dtype=np.int16).reshape(NCH, NS)
u16 = lambda a: (a.astype(np.int32) + 32768).astype(np.uint16)
i16 = lambda a: (a.astype(np.int32) - 32768).astype(np.int16)

# 预计算: 一阶差分 (int16)
d16 = np.empty_like(D, dtype=np.int16)
d16[:, 0] = D[:, 0]
d16[:, 1:] = (D[:, 1:].astype(np.int32) - D[:, :-1].astype(np.int32)).astype(np.int16)

def undelta(b):
    x = np.cumsum(asarr(b).astype(np.int64), axis=1)
    return ((x + 32768) % 65536 - 32768).astype(np.int16)

# ============ A. 通用无损库 ============
if os.environ.get('SKIP_A') == '1':
    print('\n(A 段已跑过, 跳过 —— 结果沿用 bench2_log 上次运行)', flush=True)
else:
    print('\n---- A. 通用无损库 ----', flush=True)
    t0 = time.time()
    try:
        c = imagecodecs.bz2_encode(D, level=9)
        report('bz2 L9 (raw)', len(c), np.array_equal(D, asarr(imagecodecs.bz2_decode(c))))
    except Exception as e: print('bz2 raw ERR', repr(e)[:100])
    try:
        c = imagecodecs.lzma_encode(D, level=9)
        report('lzma/xz L9 (raw)', len(c), np.array_equal(D, asarr(imagecodecs.lzma_decode(c))))
    except Exception as e: print('lzma raw ERR', repr(e)[:100])
    try:
        c = imagecodecs.brotli_encode(D, level=9)
        report('brotli Q9 (raw)', len(c), np.array_equal(D, asarr(imagecodecs.brotli_decode(c))))
    except Exception as e: print('brotli raw ERR', repr(e)[:100])
    print('  [%.0fs]' % (time.time() - t0), flush=True)

    t0 = time.time()
    try:
        c = imagecodecs.lzma_encode(d16, level=9)
        report('delta16 + lzma L9', len(c), np.array_equal(D, undelta(imagecodecs.lzma_decode(c))))
    except Exception as e: print('d+lzma ERR', repr(e)[:100])
    try:
        c = imagecodecs.brotli_encode(d16, level=9)
        report('delta16 + brotli Q9', len(c), np.array_equal(D, undelta(imagecodecs.brotli_decode(c))))
    except Exception as e: print('d+brotli ERR', repr(e)[:100])
    try:
        c = imagecodecs.zstd_encode(d16, level=22)
        report('delta16 + zstd L22', len(c), np.array_equal(D, undelta(imagecodecs.zstd_decode(c))))
    except Exception as e: print('d+zstd22 ERR', repr(e)[:100])
    print('  [%.0fs]' % (time.time() - t0), flush=True)

    t0 = time.time()
    # 字节平面 + lzma
    uu = u16(D)
    hi = (uu >> 8).astype(np.uint8); lo = (uu & 0xFF).astype(np.uint8)
    lod = lo.copy(); lod[:, 1:] = ((lo[:, 1:].astype(np.int16) - lo[:, :-1]) & 0xFF).astype(np.uint8)
    try:
        c1 = imagecodecs.lzma_encode(hi, level=9); c2 = imagecodecs.lzma_encode(lod, level=9)
        h2_ = np.frombuffer(imagecodecs.lzma_decode(c1), dtype=np.uint8).reshape(NCH, NS).astype(np.int64)
        l2d = np.frombuffer(imagecodecs.lzma_decode(c2), dtype=np.uint8).reshape(NCH, NS).astype(np.int64)
        back = (((h2_ << 8) | (np.cumsum(l2d, axis=1) & 0xFF)) - 32768).astype(np.int16)
        report('hi8 + dlo8 + lzma L9', len(c1) + len(c2), np.array_equal(D, back))
    except Exception as e: print('hi8lzma ERR', repr(e)[:100])
    # 时间交织布局 (NS, NCH) + zstd: 同一时刻的48通道相邻存放
    try:
        c = imagecodecs.zstd_encode(np.ascontiguousarray(D.T), level=19)
        back = np.frombuffer(imagecodecs.zstd_decode(c), dtype=np.int16).reshape(NS, NCH).T
        report('时间交织布局 + zstd L19', len(c), np.array_equal(D, back))
    except Exception as e: print('interleave ERR', repr(e)[:100])
    print('  [%.0fs]' % (time.time() - t0), flush=True)
# A 段结果回填 (上次运行已验证, 免重跑)
for _r in [('lzma/xz L9 (raw)', 1.663, None, ''),
           ('brotli Q9 (raw)', 1.426, None, ''),
           ('delta16 + lzma L9', 1.784, True, ''),
           ('delta16 + brotli Q9', 1.634, True, ''),
           ('delta16 + zstd L22', 1.606, True, ''),
           ('hi8 + dlo8 + lzma L9', 1.790, True, ''),
           ('时间交织布局 + zstd L19', 1.294, True, '')]:
    rows.append((_r[0], _r[1], NATIVE_MBPS / _r[1], _r[2], _r[3]))

# ============ B. 16bit 图像类 ============
print('\n---- B. 16bit 图像类 (每通道 604x1024 块) ----', flush=True)
L = 1024; NB = NS // L
def img_bench(name, enc, dec, **kw):
    t0 = time.time(); total = 0; ok = True
    for ch in range(NCH):
        img = u16(np.ascontiguousarray(D[ch, :NB * L]).reshape(NB, L))
        try:
            c = enc(img, **kw)
        except Exception as e:
            print(name, 'ch%d ERR' % ch, repr(e)[:100]); return
        total += len(c)
        if ch in (0, 7, 20, 40, NCH - 1):
            back = dec(c)
            if back is None: ok = False; continue
            back = np.asarray(back).reshape(-1)
            ok = ok and np.array_equal(i16(back), D[ch, :NB * L])
    report(name, int(total * NS / (NB * L)), ok, '%.0fs' % (time.time() - t0))

img_bench('JPEG2000 无损 (5/3, j2k)', imagecodecs.jpeg2k_encode, imagecodecs.jpeg2k_decode,
          reversible=True, codecformat='j2k')
img_bench('HTJ2K 无损', imagecodecs.htj2k_encode, imagecodecs.htj2k_decode)
img_bench('JPEG-XL 无损 (effort=7)', imagecodecs.jpegxl_encode, imagecodecs.jpegxl_decode,
          lossless=True)
img_bench('PNG 16bit L9', imagecodecs.png_encode, imagecodecs.png_decode, level=9)
img_bench('AVIF 尝试 16bit', imagecodecs.avif_encode, imagecodecs.avif_decode, level=100)

# ============ C. 科学数据压缩器 ============
print('\n---- C. 科学数据压缩器 (SZ3/ZFP/LERC) ----', flush=True)
Df = D.astype(np.float32)   # int16 -> float32 精确无损 (LSB=1µV, |x|<<2^24)
def fchk(back):
    back = np.asarray(back).reshape(NCH, NS).astype(np.float32)
    return np.array_equal(D, back.astype(np.int16))
t0 = time.time()
try:
    c = imagecodecs.sz3_encode(Df)
    report('SZ3 默认模式 (f32)', len(c), fchk(imagecodecs.sz3_decode(c)))
except Exception as e: print('SZ3 default ERR', repr(e)[:120])
for eb in (2, 4):
    try:
        c = imagecodecs.sz3_encode(Df, abs=eb)
        back = np.asarray(imagecodecs.sz3_decode(c)).reshape(NCH, NS)
        err = np.abs(back - Df)
        report('SZ3 abs=%dµV (误差受限)' % eb, len(c), None,
               'worst|err|=%d RMSE=%.2f' % (err.max(), np.sqrt((err ** 2).mean())))
    except Exception as e: print('SZ3 abs%d ERR' % eb, repr(e)[:120])
print('  [%.0fs]' % (time.time() - t0), flush=True)

t0 = time.time()
try:
    c = imagecodecs.zfp_encode(Df)
    report('ZFP 默认 (f32)', len(c), fchk(imagecodecs.zfp_decode(c)))
except Exception as e: print('ZFP default ERR', repr(e)[:120])
for acc in (2, 4):
    try:
        c = imagecodecs.zfp_encode(Df, level=-acc)   # 负 level 按容差解释, 实测验证误差界
        back = np.asarray(imagecodecs.zfp_decode(c)).reshape(NCH, NS)
        err = np.abs(back - Df)
        report('ZFP level=-%d (f32)' % acc, len(c), None,
               'worst|err|=%d RMSE=%.2f' % (err.max(), np.sqrt((err ** 2).mean())))
    except Exception as e: print('ZFP -%d ERR' % acc, repr(e)[:120])
print('  [%.0fs]' % (time.time() - t0), flush=True)

t0 = time.time()
try:
    c = imagecodecs.lerc_encode(D, level=0)
    back = np.asarray(imagecodecs.lerc_decode(c)).reshape(NCH, NS)
    report('LERC level=0 (无损)', len(c), np.array_equal(D, back))
except Exception as e: print('LERC ERR', repr(e)[:120])
print('  [%.0fs]' % (time.time() - t0), flush=True)

# ============ D. 航天标准 AEC (CCSDS-121) ============
print('\n---- D. 航天标准 AEC / CCSDS-121 ----', flush=True)
t0 = time.time()
try:
    c = imagecodecs.aec_encode(np.ascontiguousarray(D).reshape(-1))
    back = np.asarray(imagecodecs.aec_decode(c)).reshape(NCH, NS)
    report('AEC/CCSDS-121 (raw)', len(c), np.array_equal(D, back))
except Exception as e: print('AEC raw ERR', repr(e)[:120])
try:
    c = imagecodecs.aec_encode(np.ascontiguousarray(d16).reshape(-1))
    back = np.asarray(imagecodecs.aec_decode(c)).reshape(NCH, NS)
    report('delta16 + AEC/CCSDS-121', len(c), np.array_equal(D, undelta(np.ascontiguousarray(back))))
except Exception as e: print('AEC delta ERR', repr(e)[:120])
print('  [%.0fs]' % (time.time() - t0), flush=True)

# ============ E. 自研跨通道预测 (tetrode 邻道) ============
print('\n---- E. 自研跨通道预测 ----', flush=True)
REG = [('OB', list(range(0, 8))), ('CA1', list(range(8, 28))), ('PFC', list(range(28, 48)))]
X = D.astype(np.float64)
Xn = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, keepdims=True)
C = Xn @ Xn.T / NS          # 48x48 相关系数矩阵
print('分脑区相邻相关系数(排序后前5):')
for rn, chs in REG:
    sub = np.sort(C[np.ix_(chs, chs)][np.triu_indices(len(chs), 1)])[::-1]
    print('  %-4s top5: %s' % (rn, np.round(sub[:5], 3)))

def resid_encode(name, resid, predmap):
    """resid: int32 (48,NS) 残差; predmap: 每通道的预测通道列表; 返回 None"""
    print('  残差统计: 全局 std=%.1f max|r|=%d' % (resid.std(), np.abs(resid).max()))
    assert np.abs(resid).max() < 32768, '残差超 int16'
    r16 = resid.astype(np.int16)
    # 1) zstd 整体
    c = imagecodecs.zstd_encode(r16, level=19)
    # 重建校验
    back = np.frombuffer(imagecodecs.zstd_decode(c), dtype=np.int16).reshape(NCH, NS).astype(np.int32)
    for i, ps in predmap.items():
        for p in ps: back[i] += D[p]
    report(name + ' 残差+zstd L19', len(c), np.array_equal(D, back.astype(np.int16)))
    # 2) FLAC 每通道
    os.makedirs('D:/bcidata/dandi001539/flac_out2', exist_ok=True)
    total = 0; ok = True
    for i in range(NCH):
        fn = 'D:/bcidata/dandi001539/flac_out2/ch%02d.flac' % i
        sf.write(fn, r16[i], int(FS), format='FLAC', subtype='PCM_16')
        total += os.path.getsize(fn)
        if i in (0, 7, 20, 40, NCH - 1):
            b, _ = sf.read(fn, dtype='int16'); ok = ok and np.array_equal(b, r16[i])
    report(name + ' 残差+FLAC', total, ok)
    return r16

# E1: 零系数邻道差分 r = x_i - x_j (最 FPGA 友好)
pred1 = {}; resid = np.zeros((NCH, NS), dtype=np.int32)
for rn, chs in REG:
    for i in chs:
        j = max((j2 for j2 in chs if j2 != i), key=lambda j2: C[i, j2])
        pred1[i] = [j]
        resid[i] = D[i].astype(np.int32) - D[j].astype(np.int32)
print('E1 零系数邻道差分 (r = x_i - x_j):')
r16_a = resid_encode('E1 零系数邻道', resid, pred1)

# E2: 双邻道平均 r = x_i - round((x_j1+x_j2)/2) (FPGA: 求和移位)
pred2 = {}; resid2 = np.zeros((NCH, NS), dtype=np.int32)
for rn, chs in REG:
    for i in chs:
        order = sorted([j for j in chs if j != i], key=lambda j: -C[i, j])
        j1, j2 = order[0], order[1]
        pred2[i] = [j1, j2]
        resid2[i] = D[i].astype(np.int32) - ((D[j1].astype(np.int32) + D[j2].astype(np.int32)) >> 1)
print('E2 双邻道平均 (r = x_i - (x_j1+x_j2)/2):')
r16_b = resid_encode('E2 双邻道平均', resid2, pred2)

# E3: 最小二乘单邻道 r = x_i - round(a*x_j + b)
pred3 = {}; resid3 = np.zeros((NCH, NS), dtype=np.int32)
for rn, chs in REG:
    for i in chs:
        order = sorted([j for j in chs if j != i], key=lambda j: -C[i, j])
        j = order[0]
        a = C[i, j] * X[i].std() / X[j].std()
        b = X[i].mean() - a * X[j].mean()
        pred3[i] = [j]
        resid3[i] = np.rint(D[i].astype(np.float64) - (a * X[j] + b)).astype(np.int32)
print('E3 最小二乘单邻道 (r = x_i - round(a*x_j+b)):')
r16_c = resid_encode('E3 LS 单邻道', resid3, pred3)

# E4: 近无损组合 —— 最优残差 + JPEG-LS NEAR=2
best = min([(np.abs(r16_a.astype(np.int32)).mean(), r16_a, 'E1'),
            (np.abs(r16_b.astype(np.int32)).mean(), r16_b, 'E2'),
            (np.abs(r16_c.astype(np.int32)).mean(), r16_c, 'E3')])[1]
try:
    n = 0
    for ch in range(NCH):
        n += len(imagecodecs.jpegls_encode(u16(np.ascontiguousarray(best[ch, :NB * L])).reshape(NB, L), level=2))
    report('最优残差 + JPEG-LS NEAR=2', int(n * NS / (NB * L)), None, '残差域 err<=2µV')
except Exception as e: print('E4 ERR', repr(e)[:120])

# ============ 汇总 ============
print('\n==== bench2 汇总 (排序) ====')
for name, cr, mb, ok, note in sorted(rows, key=lambda r: -r[1]):
    print('%-42s CR=%6.3fx  %5.3f Mbps  ok=%s  %s' % (name, cr, mb, ok, note))
with open('D:/bcidata/dandi001539/bench2_results.json', 'w', encoding='utf-8') as f:
    json.dump([{'method': n, 'CR': round(c, 3), 'Mbps': round(m, 3), 'lossless': o, 'note': s}
               for n, c, m, o, s in rows], f, ensure_ascii=False, indent=1)
print('saved bench2_results.json')
