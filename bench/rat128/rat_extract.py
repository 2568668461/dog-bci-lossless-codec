# -*- coding: utf-8 -*-
"""大鼠皮层在体数据 (Horvath 2021, Intan RHD-2000, 128ch, 20kHz, 16bit) 压缩测试
用法: PYTHONPATH=vendor_h5py python rat_test.py [秒数] [起始秒]
"""
import sys, os
import numpy as np

DUR = int(sys.argv[1]) if len(sys.argv) > 1 else 60      # 提取时长（秒）
START = int(sys.argv[2]) if len(sys.argv) > 2 else 60     # 起始偏移（秒，跳过开头的稳定段）

import h5py
f = h5py.File('Rat06_Insertion2_Depth1.nwb', 'r')

def walk(g, prefix=''):
    for k in g:
        item = g[k]
        if isinstance(item, h5py.Dataset):
            print(f'{prefix}/{k}  shape={item.shape} dtype={item.dtype} compression={item.compression}')
        else:
            walk(item, prefix + '/' + k)

print('== 文件结构 ==')
walk(f)
print()

# 定位宽频连续数据
path = None
def find(g, prefix=''):
    global path
    for k in g:
        item = g[k]
        if isinstance(item, h5py.Dataset):
            if 'wideband' in k or ('recording' in k and item.ndim == 2):
                path = prefix + '/' + k
        else:
            find(item, prefix + '/' + k)
find(f)
print('数据集路径:', path)
d = f[path]
print('shape:', d.shape, 'dtype:', d.dtype)
print('attrs:', dict(d.attrs))

fs = None
for cand in ['rate', 'sampling_rate', 'starting_time']:
    if cand in d.attrs:
        print('attr', cand, '=', d.attrs[cand])
# 尝试从 general 找采样率
try:
    for k, v in f['general'].attrs.items():
        print('general attr:', k, v)
except Exception as e:
    print('general attrs err', e)

n_ch = d.shape[0] if d.shape[0] < d.shape[1] else d.shape[1]
print('通道数:', n_ch)

# 数据布局: (ch, time) 或 (time, ch)
if d.shape[0] < d.shape[1]:
    sl = (slice(None), slice(START*20000, (START+DUR)*20000))
else:
    sl = (slice(START*20000, (START+DUR)*20000), slice(None))
x = d[sl]
print('提取块 shape:', x.shape, x.dtype, 'range:', x.min(), x.max())

if x.dtype.kind == 'f':
    # float 存储：找换算系数并恢复 int16
    conv = 1.0
    for k in ['conversion']:
        if k in d.attrs:
            conv = float(d.attrs[k])
    print('conversion:', conv)
    # volts -> uV -> LSB
    raw = np.round(x / (0.195e-6 / conv if conv else 0.195e-6)).astype(np.int16) if False else None
    # 简化: 假设 float 已是原始计数*conv; 若 conv 使数值处于 V 量级则除以 0.195uV
    scale = abs(x).max()
    if scale < 1.0:   # 伏特量级
        raw = np.round(x / 0.195e-6).astype(np.int32).astype(np.int16)
    else:             # 已经是计数
        raw = np.round(x).astype(np.int16)
else:
    raw = x.astype(np.int16)
np.save('rat128_raw.npy', raw)
print('已保存 rat128_raw.npy', raw.shape)
print('每通道 std (LSB):', np.round(raw.std(axis=1)[:10], 1), '...')
