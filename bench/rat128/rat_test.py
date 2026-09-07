# -*- coding: utf-8 -*-
"""大鼠皮层在体数据 (Horvath 2021, Intan RHD-2000, 128ch, 20kHz, 16bit) 压缩测试
在 D:/bcidata 运行: PYTHONPATH=C:/Users/MI/WorkBuddy/2026-09-04-20-43-57/intan-data/vendor_h5py python rat_test.py
"""
import sys, os
import numpy as np

sys.path.insert(0, 'C:/Users/MI/WorkBuddy/2026-09-04-20-43-57/intan-data/vendor_h5py')
import h5py

f = h5py.File('D:/bcidata/Rat06_Insertion2_Depth1.nwb', 'r')

def walk(g, prefix=''):
    out = []
    for k in g:
        item = g[k]
        if isinstance(item, h5py.Dataset):
            out.append(f'{prefix}/{k}  shape={item.shape} dtype={item.dtype} comp={item.compression}')
        else:
            out.extend(walk(item, prefix + '/' + k))
    return out

print('== 文件结构 ==')
for line in walk(f):
    print(' ', line)

# 找宽频连续数据
path = None
def find(g, prefix=''):
    global path
    for k in g:
        item = g[k]
        if isinstance(item, h5py.Dataset):
            if item.ndim == 2 and ('wideband' in k or 'recording' in k):
                path = prefix + '/' + k
        else:
            find(item, prefix + '/' + k)
find(f)
print('数据集:', path)
d = f[path]
print('shape:', d.shape, 'dtype:', d.dtype, 'compression:', d.compression)
print('attrs:', {k: (v if not hasattr(v, 'shape') else '...') for k, v in d.attrs.items()})
print('rate attr:', d.attrs.get('rate'), d.attrs.get('starting_time'))

# 提取中间 60 秒（跳过开头稳定段），128 通道全量
FS = 20000
START_SEC, DUR_SEC = 300, 60
if d.shape[0] < d.shape[1]:
    n_ch, n_t = d.shape
    sl = (slice(None), slice(START_SEC*FS, (START_SEC+DUR_SEC)*FS))
else:
    n_t, n_ch = d.shape
    sl = (slice(START_SEC*FS, (START_SEC+DUR_SEC)*FS), slice(None))
print(f'总时长 {n_t/FS/60:.1f} min, 提取 {START_SEC}s 起的 {DUR_SEC}s x {n_ch}ch')
x = d[sl]
print('提取块:', x.shape, x.dtype, 'range:', x.min(), x.max())

if x.dtype.kind == 'f':
    conv = float(d.attrs.get('conversion', 1.0) or 1.0)
    off = float(d.attrs.get('offset', 0.0) or 0.0)
    print('conversion:', conv, 'offset:', off)
    # 若为伏特, 0.195 uV/LSB
    v = x * conv + off
    raw = np.round(v / 0.195e-6).astype(np.int16)
else:
    raw = x.astype(np.int16)
np.save('D:/bcidata/rat128_raw.npy', raw)
print('每通道 std (LSB):', np.round(raw.std(axis=1)[:12], 1))
print('全局 range:', raw.min(), raw.max())
c = np.corrcoef(raw[:8].astype(float))
print('ch0-7 互相关 |r| 均值: %.3f' % np.abs(c[~np.eye(8,dtype=bool)]).mean())
print('SAVED rat128_raw.npy', raw.shape)
