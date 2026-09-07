# -*- coding: utf-8 -*-
"""FLAC 修正版：直接写 int16（不走 float 缩放），逐通道全部校验"""
import os
import numpy as np
import soundfile as sf

D = np.load('D:/bcidata/rat128_raw.npy')
NCH, NS = D.shape
raw_bytes = D.nbytes
os.makedirs('D:/bcidata/flac_out2', exist_ok=True)

total = 0
bad = []
for ch in range(NCH):
    fn = f'D:/bcidata/flac_out2/ch{ch}.flac'
    sf.write(fn, D[ch], 20000, format='FLAC', subtype='PCM_16')  # int16 直写
    total += os.path.getsize(fn)
    b, _ = sf.read(fn, dtype='int16')
    if not np.array_equal(b, D[ch]):
        bad.append(ch)

cr = raw_bytes / total
print(f'FLAC(int16直写) CR={cr:.3f}x  {total/5e6:.1f} Mbps  全128通道校验: '
      f'{"全部无损" if not bad else f"FAIL ch={bad[:5]}"}')

# 也测 8 通道一组封装（FLAC 上限 8ch/流，但组内只做声道去相关，128ch 全组不必要）
# 上面 per-ch mono 已是最保守口径
