# -*- coding: utf-8 -*-
"""DANDI 001539 提取: 48ch LFP (OB/CA1/PFC) x 618512 samples @ 1500 Hz, µV float64 -> int16
1) 检查浮点 µV 的量化粒度(推断原始 ADC LSB)
2) 按推断 LSB 量化为 int16 存 npy
"""
import h5py, numpy as np

f = h5py.File('D:/bcidata/dandi001539/sub-CS39_ses-06.nwb', 'r')
d = f['acquisition/ElectricalSeries/data']   # (618512, 48) float64 µV
el = f['general/extracellular_ephys/electrodes']
loc = [x.decode() if isinstance(x, bytes) else str(x) for x in np.atleast_1d(el['location'][()])]
regions = np.array(loc)
print('regions:', {r: int((regions == r).sum()) for r in set(loc)})

# 全量读入（618512*48*8 = 237 MB，可承受）
X = d[:, :]            # (ns, nch) µV
X = X.T                # (nch, ns)
print('shape:', X.shape, 'range: %.1f ~ %.1f µV' % (X.min(), X.max()))

# --- 量化粒度检查：取每通道一段，看 unique 值间距 ---
sub = X[:, 100000:160000]
# 找最小非零差
for name in ['OB', 'CA1', 'PFC']:
    ch = np.where(regions == name)[0][0]
    v = np.sort(np.unique(sub[ch].round(9)))
    dv = np.diff(v)
    dv = dv[dv > 1e-9]
    if len(dv):
        print('%s ch%d: 最小步长 %.6f µV, 步长中位 %.6f µV, unique 值 %d / 60000' %
              (name, ch, dv.min(), np.median(dv), len(v)))
    else:
        print('%s ch%d: 无重复, 全连续' % (name, ch))

# --- 用最小公共步长量化 ---
# 试常见 LSB: Intan 0.195µV; 也试数据自身的最小步长
sub_all = np.sort(np.unique(X[:, ::7].round(9)))   # 抽样
dv = np.diff(sub_all); dv = dv[dv > 1e-9]
lsb_data = dv.min()
print('\n全体抽样最小步长: %.6f µV' % lsb_data)

for lsb in [lsb_data, 0.195, 0.05, 0.01]:
    q = np.round(X / lsb)
    err = np.abs(X - q * lsb).max()
    iq = q.astype(np.int64)
    print('LSB=%.4f µV: int 值域 [%d, %d], 最大量化误差 %.6f µV, int16 可容纳: %s' %
          (lsb, iq.min(), iq.max(), err, iq.min() >= -32768 and iq.max() <= 32767))
