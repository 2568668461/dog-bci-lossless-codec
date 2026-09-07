"""drice_nl_test.py — 在 rat128 真实数据上对比 4 种近无损方法
   方法1:Hub 端截 2 LSB + 原 drice（极简硬件：右移 2 位 + drice）
   方法2:drice_nl NEAR=2（闭环死区量化，等价 NEAR=2 但有界保证）
   方法3:drice_nl NEAR=8（更激进，11.6 Mbps 档）
   方法4:drice_nl 通道噪声自适应（每通道 NEAR=σ̂_c/2 LSB）
"""
import numpy as np
import subprocess, os, json, time, struct

PY = r"C:/Users/MI/AppData/Local/Programs/Python/Python312/python.exe"
DRICE = r"D:/bcidata/drice/drice.exe"
DRICE_NL = r"D:/bcidata/drice/drice_nl.exe"
DATA_NPY = r"D:/bcidata/rat128_raw.npy"
TMP = r"D:/bcidata/_drice_nl_tmp"

os.makedirs(TMP, exist_ok=True)

print("加载数据...", flush=True)
D = np.load(DATA_NPY)
NCH, NS = D.shape
TOTAL = D.nbytes
print(f"  数据: {NCH}ch × {NS} samples, total {TOTAL/1e6:.1f} MB, 时长 {NS/30000:.1f}s @30kHz")

# 估算每通道噪声 σ̂（用一阶差分的 1.4826×MAD 抗差估计）
print("估算各通道噪声 σ̂...", flush=True)
Dint = D.astype(np.int32)
d1 = np.diff(Dint, axis=1)
mad = np.median(np.abs(d1 - np.median(d1, axis=1, keepdims=True)), axis=1)
sig_lsb = (1.4826 * mad / np.sqrt(2)).astype(np.int32)
print(f"  σ̂ LSB: min={sig_lsb.min()} med={int(np.median(sig_lsb))} max={sig_lsb.max()}")
print(f"  σ̂ µV:  min={sig_lsb.min()*0.195:.2f} med={int(np.median(sig_lsb))*0.195:.2f} max={sig_lsb.max()*0.195:.2f}")


def run_drice(bin_in, near_table=None):
    """调用 drice 或 drice_nl，返回编码字节数、解码数据"""
    if near_table is None:
        exe = DRICE
        cmd = [exe, "c", bin_in, bin_in + ".out", str(NCH), "1024"]
    else:
        # per-channel NEAR 表
        nt_path = bin_in + ".near"
        with open(nt_path, "wb") as f:
            for v in near_table: f.write(struct.pack("<H", int(v)))
        exe = DRICE_NL
        cmd = [exe, "c", bin_in, bin_in + ".out", str(NCH), "1024", "--per-near=" + nt_path]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print("STDERR:", r.stderr); raise RuntimeError("encode failed")
    # 从 stderr 取 CR（"CR=..."）
    cr = None
    for line in r.stderr.splitlines():
        if "CR=" in line:
            cr = float(line.split("CR=")[1].split()[0]); break
    out_path = bin_in + ".out"
    osz = os.path.getsize(out_path)
    # 解码
    r2 = subprocess.run([exe, "d", out_path, bin_in + ".back"], capture_output=True, text=True, timeout=600)
    back = np.fromfile(bin_in + ".back", dtype=np.int16).reshape(NCH, NS).copy()
    return osz, cr, back


def report(name, osz, cr, back, near_ref=None):
    err = (back.astype(np.int32) - D.astype(np.int32))
    worst = int(np.abs(err).max())
    p99_9 = float(np.percentile(np.abs(err), 99.9))
    p99 = float(np.percentile(np.abs(err), 99))
    cr_t = TOTAL / osz
    mbps = osz / NS * 30 / 1e6 * 1000  # 编码字节 / 样本数 × 30k × 1e-6  = Mbps  ≈ osz*30e3/NS/1e6 = osz/NS*30/1e3
    mbps = osz * 30 / NS / 1000.0  # ← 该公式与 TOTAL/TOTAL * 30 = 30 Mbps 不匹配，重写
    # Mbps = osz(字节) × 8(bit) × 30kHz / NS(样本) / 1e6
    mbps = osz * 8.0 * 30000 / NS / 1e6
    cr_v = TOTAL / osz
    near_str = f" NEAR_ref={near_ref}" if near_ref is not None else ""
    print(f"\n[{name}] CR={cr_v:.3f}x  {mbps:.1f} Mbps  worst|err|={worst} LSB ({worst*0.195:.2f} µV)  p99={p99:.0f} p99.9={p99_9:.0f}{near_str}")
    return worst


results = {}

# ===== 方法 1: 原始 drice 无损（基线） =====
print("\n========== 方法 0: 原 drice 无损（基线） ==========", flush=True)
in_bin = os.path.join(TMP, "raw.bin")
D.astype(np.int16).tofile(in_bin)
t = time.time()
osz, cr, back = run_drice(in_bin)
report("M0 原 drice 无损", osz, cr, back)
results["M0_原drice无损"] = (osz, cr, 0, back)
print(f"  用时 {time.time()-t:.1f}s")

# ===== 方法 1: Hub 端截 2 LSB + drice =====
print("\n========== 方法 1: 截 2 LSB + drice ==========", flush=True)
in_bin = os.path.join(TMP, "trunc2.bin")
(D & ~3).astype(np.int16).tofile(in_bin)
t = time.time()
osz, cr, back = run_drice(in_bin)
worst = report("M1 截2LSB+drice", osz, cr, back)
# 注：误差是 drice 解码 vs 原始，理论最大 2 LSB（实际因 drice 内部一致，应该完全≤2）
results["M1_截2LSB"] = (osz, cr, worst, back)
print(f"  用时 {time.time()-t:.1f}s")

# ===== 方法 2: drice_nl NEAR=2 =====
print("\n========== 方法 2: drice_nl NEAR=2 ==========", flush=True)
in_bin2 = os.path.join(TMP, "raw.bin")  # 用原始
t = time.time()
# drice_nl 但走全局 NEAR
exe = DRICE_NL
cmd = [exe, "c", in_bin2, in_bin2 + ".out", str(NCH), "1024", "2"]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
cr = float([l for l in r.stderr.splitlines() if "CR=" in l][0].split("CR=")[1].split()[0])
out_path = in_bin2 + ".out"; osz = os.path.getsize(out_path)
subprocess.run([exe, "d", out_path, in_bin2 + ".back"], capture_output=True, text=True, timeout=600)
back = np.fromfile(in_bin2 + ".back", dtype=np.int16).reshape(NCH, NS).copy()
worst = report("M2 drice_nl NEAR=2", osz, cr, back, near_ref=2)
results["M2_drice_nl_NEAR2"] = (osz, cr, worst, back)
print(f"  用时 {time.time()-t:.1f}s")

# ===== 方法 3: drice_nl NEAR=8 =====
print("\n========== 方法 3: drice_nl NEAR=8 ==========", flush=True)
t = time.time()
cmd = [exe, "c", in_bin2, in_bin2 + ".out", str(NCH), "1024", "8"]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
cr = float([l for l in r.stderr.splitlines() if "CR=" in l][0].split("CR=")[1].split()[0])
out_path = in_bin2 + ".out"; osz = os.path.getsize(out_path)
subprocess.run([exe, "d", out_path, in_bin2 + ".back"], capture_output=True, text=True, timeout=600)
back = np.fromfile(in_bin2 + ".back", dtype=np.int16).reshape(NCH, NS).copy()
worst = report("M3 drice_nl NEAR=8", osz, cr, back, near_ref=8)
results["M3_drice_nl_NEAR8"] = (osz, cr, worst, back)
print(f"  用时 {time.time()-t:.1f}s")

# ===== 方法 4: drice_nl 通道噪声自适应（每通道 NEAR=σ̂/2） =====
print("\n========== 方法 4: drice_nl 通道噪声自适应 NEAR≈σ̂/2 ==========", flush=True)
t = time.time()
near_table = (sig_lsb / 2).clip(1, 31).astype(np.int32)  # NEAR≥1，且≤31（编码器限制）
# 写入文件
nt_path = os.path.join(TMP, "per_ch_near.bin")
with open(nt_path, "wb") as f:
    for v in near_table: f.write(struct.pack("<H", int(v)))
cmd = [exe, "c", in_bin2, in_bin2 + ".out", str(NCH), "1024", "0", "--per-near=" + nt_path]
r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
if r.returncode != 0:
    print("M4 ENCODE FAILED rc=", r.returncode)
    print("STDERR:", r.stderr); print("STDOUT:", r.stdout)
    raise RuntimeError("M4 encode failed")
cr = float([l for l in r.stderr.splitlines() if "CR=" in l][0].split("CR=")[1].split()[0])
out_path = in_bin2 + ".out"; osz = os.path.getsize(out_path)
print(f"  M4 编码输出文件大小: {osz/1e6:.1f} MB (per-ch NEAR min/max={near_table.min()}/{near_table.max()})")
subprocess.run([exe, "d", out_path, in_bin2 + ".back"], capture_output=True, text=True, timeout=600)
back = np.fromfile(in_bin2 + ".back", dtype=np.int16).reshape(NCH, NS).copy()
worst = report("M4 噪声自适应 NEAR≈σ̂/2", osz, cr, back, near_ref="per-ch σ̂/2")
results["M4_drice_nl_噪声自适应"] = (osz, cr, worst, back)
print(f"  NEAR 表: min={near_table.min()} med={int(np.median(near_table))} max={near_table.max()}")
print(f"  对应 µV: min={near_table.min()*0.195:.2f} med={int(np.median(near_table))*0.195:.2f} max={near_table.max()*0.195:.2f}")
print(f"  用时 {time.time()-t:.1f}s")

# ===== 总表 =====
print("\n" + "="*90)
print(" 总表（rat128 真实数据，128ch×60s×30kHz, 30 Mbps 原始带宽）")
print("="*90)
print(f"{'方法':<28} {'CR':>6} {'Mbps':>7} {'worst|err| LSB':>13} {'p99 LSB':>9} {'p99.9 LSB':>10}")
print("-"*90)
for name, (osz, cr, worst, back) in results.items():
    err = (back.astype(np.int32) - D.astype(np.int32))
    mbps = osz * 8.0 * 30000 / NS / 1e6
    p99 = float(np.percentile(np.abs(err), 99))
    p999 = float(np.percentile(np.abs(err), 99.9))
    cr_v = TOTAL / osz
    print(f"{name:<28} {cr_v:>6.3f} {mbps:>7.1f} {worst:>13} {p99:>9.0f} {p999:>10.0f}")

# 保存
with open(os.path.join(TMP, "results.json"), "w") as f:
    json.dump({k: {"osz": v[0], "worst_err": v[2]} for k, v in results.items()}, f, indent=2)
print(f"\n中间产物: {TMP}")