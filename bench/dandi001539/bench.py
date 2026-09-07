# -*- coding: utf-8 -*-
"""DANDI 001539 (ob48_lfp.npy) 压缩基准 —— 与 rat128 同口径对比
数据: 48ch x 618512 int16 @ 1500 Hz (清醒大鼠 OB/CA1/PFC LFP, LSB=1µV, 412 s)
原生码率: 48 * 1500 * 16 = 1.152 Mbps
"""
import os, time, subprocess, json
import numpy as np
import imagecodecs
import soundfile as sf

D = np.load('D:/bcidata/dandi001539/ob48_lfp.npy')
NCH, NS = D.shape
RAW = D.nbytes
FS, LSB = 1500.0, 1.0
NATIVE_MBPS = NCH * FS * 16 / 1e6
print('data: %s %s, raw = %.1f MB, 原生码率 %.3f Mbps' % (D.shape, D.dtype, RAW/1e6, NATIVE_MBPS))

rows = []
def report(name, nbytes, ok=None, note=''):
    cr = RAW / nbytes
    rows.append((name, cr, NATIVE_MBPS/cr, ok, note))
    print('%-38s CR=%6.3fx  %5.2f Mbps  %6.3f bits/s  ok=%s  %s' %
          (name, cr, NATIVE_MBPS/cr, 16/cr, ok, note), flush=True)

def asarr(b): return np.frombuffer(b, dtype=np.int16).reshape(NCH, NS)

# ---------- 0. 熵锚点 ----------
t0 = time.time()
v = D.astype(np.int32).ravel()
u = (v + 32768).astype(np.uint16)
h = np.bincount(u, minlength=65536).astype(float); h /= h.sum(); h = h[h>0]
e0 = -(h*np.log2(h)).sum()
d1 = np.diff(D.astype(np.int32), axis=1).ravel()
h1 = np.bincount((d1+2048).clip(0,4095), minlength=4096).astype(float); h1 /= h1.sum(); h1 = h1[h1>0]
e1 = -(h1*np.log2(h1)).sum()
d2 = np.diff(D.astype(np.int32), axis=1, n=2).ravel()
h2 = np.bincount((d2+4096).clip(0,8191), minlength=8192).astype(float); h2 /= h2.sum(); h2 = h2[h2>0]
e2 = -(h2*np.log2(h2)).sum()
report('熵锚点: 0阶熵', int(RAW*e0/16), note='%.2f bits/样本' % e0)
report('熵锚点: 1阶差分熵', int(RAW*e1/16), note='%.2f bits/样本' % e1)
report('熵锚点: 2阶差分熵', int(RAW*e2/16), note='%.2f bits/样本' % e2)
print('  [%.0fs]' % (time.time()-t0))

# ---------- 1. 通用无损库 ----------
t0 = time.time()
c = imagecodecs.zlib_encode(D, level=9)
report('DEFLATE/zlib L9 (raw)', len(c), np.array_equal(D, asarr(imagecodecs.zlib_decode(c))))
d16 = np.empty_like(D, dtype=np.int32)
d16[:, 0] = D[:, 0]; d16[:, 1:] = D[:, 1:].astype(np.int32) - D[:, :-1].astype(np.int32)
d16 = d16.astype(np.int16)
c = imagecodecs.zlib_encode(d16, level=9)
back = np.cumsum(asarr(imagecodecs.zlib_decode(c)).astype(np.int32), axis=1)
back = ((back + 32768) % 65536 - 32768).astype(np.int16)
report('delta16 + zlib L9', len(c), np.array_equal(D, back))
# 字节平面拆分
uu = (D.astype(np.int32) + 32768).astype(np.uint16)
hi = (uu >> 8).astype(np.uint8); lo = (uu & 0xFF).astype(np.uint8)
lod = lo.copy(); lod[:, 1:] = ((lo[:, 1:].astype(np.int16) - lo[:, :-1]) & 0xFF).astype(np.uint8)
c1 = imagecodecs.zlib_encode(hi, level=9); c2 = imagecodecs.zlib_encode(lod, level=9)
h2_ = np.frombuffer(imagecodecs.zlib_decode(c1), dtype=np.uint8).reshape(NCH, NS)
l2d = np.frombuffer(imagecodecs.zlib_decode(c2), dtype=np.uint8).reshape(NCH, NS).astype(np.int64)
lo_rec = np.cumsum(l2d, axis=1) & 0xFF   # 累积和还原 delta 低字节平面
back = (((h2_.astype(np.int32) << 8) | lo_rec.astype(np.int32)) - 32768).astype(np.int16)
report('hi8 + dlo8 + zlib L9', len(c1)+len(c2), np.array_equal(D, back))
print('  [%.0fs]' % (time.time()-t0))

t0 = time.time()
for lvl in (3, 19):
    c = imagecodecs.zstd_encode(D, level=lvl)
    report('zstd L%d' % lvl, len(c), np.array_equal(D, asarr(imagecodecs.zstd_decode(c))))
c = imagecodecs.zstd_encode(d16, level=19)
back = np.cumsum(asarr(imagecodecs.zstd_decode(c)).astype(np.int32), axis=1)
back = ((back + 32768) % 65536 - 32768).astype(np.int16)
report('delta16 + zstd L19', len(c), np.array_equal(D, back))
try:
    c = imagecodecs.lz4_encode(D)
    report('LZ4', len(c), np.array_equal(D, asarr(imagecodecs.lz4_decode(c))))
except Exception as e:
    print('LZ4 skipped:', repr(e)[:80])
print('  [%.0fs]' % (time.time()-t0))

# ---------- 2. FLAC ----------
t0 = time.time()
os.makedirs('D:/bcidata/dandi001539/flac_out', exist_ok=True)
total = 0; ok = True
for ch in range(NCH):
    fn = 'D:/bcidata/dandi001539/flac_out/ch%d.flac' % ch
    sf.write(fn, D[ch], int(FS), format='FLAC', subtype='PCM_16')   # int16 直写
    total += os.path.getsize(fn)
    if ch in (0, 20, NCH-1):
        b, _ = sf.read(fn, dtype='int16'); ok = ok and np.array_equal(b, D[ch])
report('FLAC (48ch mono, L5)', total, ok)
print('  [%.0fs]' % (time.time()-t0))

# ---------- 3. JPEG-LS (imagecodecs/CharLS) ----------
t0 = time.time()
L = 1024; NB = NS // L   # 604 块/通道, 覆盖 99.94%
def u16(a): return (a.astype(np.int32) + 32768).astype(np.uint16)
for lv in [0, 1, 2, 4]:
    n = 0
    for c in range(NCH):
        n += len(imagecodecs.jpegls_encode(u16(np.ascontiguousarray(D[c, :NB*L])).reshape(NB, L), level=lv))
    if lv == 0:
        # 无损校验
        blk = u16(np.ascontiguousarray(D[7, :NB*L])).reshape(NB, L)
        back = imagecodecs.jpegls_decode(imagecodecs.jpegls_encode(blk, level=0))
        ok = np.array_equal(back.reshape(-1).astype(np.int32) - 32768, D[7][:NB*L])
        report('JPEG-LS NEAR=0 (无损)', int(n * NS / (NB*L)), ok)
    else:
        report('JPEG-LS NEAR=%d (近无损)' % lv, int(n * NS / (NB*L)), None, 'err<=%dµV' % lv)
print('  [%.0fs]' % (time.time()-t0))

# ---------- 4. WavPack ----------
t0 = time.time()
wv = 'D:/bcidata/wavpack_bin/wavpack.exe'
wav_in = 'D:/bcidata/dandi001539/full48.wav'
wv_out = 'D:/bcidata/dandi001539/full48.wv'
sf.write(wav_in, D.T, int(FS), format='WAV', subtype='PCM_16')   # int16 直写
r = subprocess.run([wv, '-hhx', '-q', wav_in, '-o', wv_out], capture_output=True, text=True, timeout=600)
osz = os.path.getsize(wv_out) if os.path.exists(wv_out) else -1
# 校验
r2 = subprocess.run(['D:/bcidata/wavpack_bin/wvunpack.exe', '-q', wv_out, 'D:/bcidata/dandi001539/full48_back.wav'],
                    capture_output=True, text=True, timeout=600)
b, sr = sf.read('D:/bcidata/dandi001539/full48_back.wav', dtype='int16')
ok = np.array_equal(b.T, D)
report('WavPack 5.9 -hhx', osz, ok)
print('  [%.0fs]' % (time.time()-t0))

# ---------- 5. 自研 drice / drice_nl ----------
t0 = time.time()
TMP = 'D:/bcidata/dandi001539/_tmp'
os.makedirs(TMP, exist_ok=True)
in_bin = TMP + '/raw.bin'
D.astype(np.int16).tofile(in_bin)
def run_drice(exe, args, tag):
    out = TMP + '/%s.bin' % tag
    cmd = [exe, 'c', in_bin, out] + args
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print('ERR', tag, r.stderr[:200]); return None, None
    osz = os.path.getsize(out)
    back = TMP + '/%s.back' % tag
    r2 = subprocess.run([exe, 'd', out, back], capture_output=True, text=True, timeout=600)
    ok = None
    if r2.returncode == 0:
        ok = bool(np.array_equal(np.fromfile(back, dtype=np.int16).reshape(NCH, NS), D))
    return osz, ok

osz, ok = run_drice('D:/bcidata/drice/drice.exe', ['48'], 'drice')
report('drice (自研, 无损)', osz, ok)
# drice_nl NEAR=2 / 4
for near in [2, 4]:
    osz, ok = run_drice('D:/bcidata/drice/drice_nl.exe', ['48', '1024', str(near)], 'dnl%d' % near)
    # 近无损: 校验误差上界
    out = TMP + '/dnl%d.bin' % near
    back = TMP + '/dnl%d.back' % near
    if os.path.exists(back):
        err = np.abs(np.fromfile(back, dtype=np.int16).reshape(NCH, NS).astype(np.int32) - D.astype(np.int32)).max()
        report('drice_nl NEAR=%d' % near, osz, None, 'worst|err|=%d' % err)
    else:
        report('drice_nl NEAR=%d' % near, osz, None, '')
print('  [%.0fs]' % (time.time()-t0))

# ---------- 汇总 ----------
print('\n==== 汇总 (排序) ====')
for name, cr, mb, ok, note in sorted(rows, key=lambda r: -r[1]):
    print('%-38s CR=%6.3fx  %5.2f Mbps  %s' % (name, cr, mb, note))
with open('D:/bcidata/dandi001539/bench_results.json', 'w', encoding='utf-8') as f:
    json.dump([{'method': n, 'CR': round(c,3), 'Mbps': round(m,3), 'lossless': o, 'note': s} for n,c,m,o,s in rows],
              f, ensure_ascii=False, indent=1)
print('saved bench_results.json')
