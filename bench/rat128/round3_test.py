# -*- coding: utf-8 -*-
"""第三轮：通用压缩库 + FLAC 在真实 Intan 128ch 数据上的无损测试。
数据: rat128_raw.npy (128, 1200000) int16, 大鼠新皮层 20kHz x 60s
Mbps 换算: 数据 60s@20kHz; 目标系统 30kS/s。压缩后 Mbps = comp_bytes*8/(128*1.2e6)*128*30000/1e6
         = comp_bytes / 5000   （与 encode16w CR2.38->25.8Mbps 口径一致）
"""
import os, time
import numpy as np
import imagecodecs

D = np.load('D:/bcidata/rat128_raw.npy')
NCH, NS = D.shape
raw_bytes = D.nbytes
print(f'data: {D.shape} {D.dtype}, raw = {raw_bytes/1e6:.1f} MB')

results = []
def mbps(cb): return cb / 5e6
def report(name, cb, ok):
    cr = raw_bytes / cb
    results.append((name, cr, mbps(cb), ok))
    print(f'{name:36s} CR={cr:6.3f}x  {mbps(cb):6.1f} Mbps  lossless={ok}')

def asarr(b): return np.frombuffer(b, dtype=np.int16).reshape(NCH, NS)

t0 = time.time()
# ---------- 1. zlib / DEFLATE (对应 WangXuan95/FPGA-GZIP-compressor 的算法) ----------
for lvl in (6, 9):
    c = imagecodecs.zlib_encode(D, level=lvl)
    d = imagecodecs.zlib_decode(c)
    report(f'DEFLATE/zlib L{lvl} (raw)', len(c), np.array_equal(D, asarr(d)))
print(f'  [{time.time()-t0:.0f}s]')

# ---------- 2. delta 预处理 + DEFLATE（delta 层在 FPGA 上几乎免费） ----------
d16 = np.empty_like(D, dtype=np.int32)
d16[:, 0] = D[:, 0]; d16[:, 1:] = D[:, 1:].astype(np.int32) - D[:, :-1].astype(np.int32)
d16 = d16.astype(np.int16)
c = imagecodecs.zlib_encode(d16, level=9)
d = imagecodecs.zlib_decode(c)
back = np.cumsum(asarr(d).astype(np.int32), axis=1)
back = ((back + 32768) % 65536 - 32768).astype(np.int16)
report('delta16 + DEFLATE/zlib L9', len(c), np.array_equal(D, back))

# 8-bit 字节域拆分: 高字节平面 + 低字节差分平面（8bit FPGA 友好）
u = (D.astype(np.int32) + 32768).astype(np.uint16)
hi = (u >> 8).astype(np.uint8)
lo = (u & 0xFF).astype(np.uint8)
lod = lo.copy()
lod[:, 1:] = ((lo[:, 1:].astype(np.int16) - lo[:, :-1]) & 0xFF).astype(np.uint8)
c1 = imagecodecs.zlib_encode(hi, level=9); c2 = imagecodecs.zlib_encode(lod, level=9)
h2 = np.frombuffer(imagecodecs.zlib_decode(c1), dtype=np.uint8).reshape(NCH, NS)
l2 = np.frombuffer(imagecodecs.zlib_decode(c2), dtype=np.uint8).reshape(NCH, NS).copy()
l2[:, 1:] = ((l2[:, 1:].astype(np.int16) + l2[:, :-1]) & 0xFF).astype(np.uint8)
back = (((h2.astype(np.int32) << 8) | l2.astype(np.int32)) - 32768).astype(np.int16)
report('hi8 + dlo8 + DEFLATE L9 (split)', len(c1) + len(c2), np.array_equal(D, back))

# ---------- 3. zstd（软件上界参考，FPGA 极难） ----------
for lvl in (3, 19):
    t1 = time.time()
    c = imagecodecs.zstd_encode(D, level=lvl)
    d = imagecodecs.zstd_decode(c)
    report(f'zstd L{lvl}', len(c), np.array_equal(D, asarr(d)))
    print(f'    ({time.time()-t1:.0f}s)')

# ---------- 4. LZ4（结构最简单，FPGA 有实现） ----------
c = imagecodecs.lz4_encode(D)
d = imagecodecs.lz4_decode(c)
report('LZ4', len(c), np.array_equal(D, asarr(d)))

# ---------- 5. FLAC（LPC+Rice；FLAC 单流最多 8ch，逐通道 mono 编码） ----------
import soundfile as sf
os.makedirs('D:/bcidata/flac_out', exist_ok=True)
total = 0; ok = True
t1 = time.time()
for ch in range(NCH):
    fn = f'D:/bcidata/flac_out/ch{ch}.flac'
    sf.write(fn, D[ch].astype(np.float64), 20000, format='FLAC', subtype='PCM_16')
    total += os.path.getsize(fn)
    if ch in (0, 1, NCH - 1):
        b, _ = sf.read(fn, dtype='int16')
        ok = ok and np.array_equal(b, D[ch])
report('FLAC (per-ch mono, LPC+Rice)', total, ok)
print(f'    ({time.time()-t1:.0f}s)')

print('\n===== 汇总 =====')
for name, cr, mb, okk in sorted(results, key=lambda r: -r[1]):
    print(f'{name:36s} CR={cr:6.3f}x  {mb:6.1f} Mbps  {"OK" if okk else "FAIL"}')
