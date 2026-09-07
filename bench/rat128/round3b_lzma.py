# -*- coding: utf-8 -*-
"""LZMA 补测（对应 WangXuan95/FPGA-LZMA-compressor 的算法族）"""
import time, lzma
import numpy as np

D = np.load('D:/bcidata/rat128_raw.npy')
raw = D.tobytes()
print(f'data: {D.shape}, raw = {len(raw)/1e6:.1f} MB')

t0 = time.time()
# preset 1 ≈ 7zip fast（FPGA 版对标 fast）；preset 6 = 默认
for preset, tag in ((1, 'fast(对标FPGA版)'), (6, 'default')):
    c = lzma.compress(raw, preset=preset)
    d = lzma.decompress(c)
    ok = np.array_equal(np.frombuffer(d, dtype=np.int16), D.ravel())
    print(f'LZMA {tag:16s} CR={len(raw)/len(c):6.3f}x  {len(c)/5e6:6.1f} Mbps  lossless={ok}  [{time.time()-t0:.0f}s]')

# delta + LZMA
d16 = np.empty_like(D, dtype=np.int32)
d16[:, 0] = D[:, 0]; d16[:, 1:] = D[:, 1:].astype(np.int32) - D[:, :-1].astype(np.int32)
c = lzma.compress(d16.astype(np.int16).tobytes(), preset=1)
print(f'LZMA fast + delta16     CR={len(raw)/len(c):6.3f}x  {len(c)/5e6:6.1f} Mbps')
