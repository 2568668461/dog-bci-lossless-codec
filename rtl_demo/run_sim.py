"""Run Icarus simulation and compare emitted bits with the golden model."""

from pathlib import Path
import shutil
import subprocess

from golden_model import encode_bits

ROOT = Path(__file__).resolve().parent
iverilog = shutil.which("iverilog")
vvp = shutil.which("vvp")
if not iverilog or not vvp:
    raise SystemExit(
        "Icarus Verilog is not installed or not on PATH. "
        "Install it, then rerun: python rtl_demo/run_sim.py"
    )

sim = ROOT / "sim.out"
subprocess.run(
    [iverilog, "-g2012", "-o", str(sim), str(ROOT / "rice_demo_encoder.v"), str(ROOT / "tb_rice_demo_encoder.v")],
    check=True,
    cwd=ROOT,
)
subprocess.run([vvp, str(sim)], check=True, cwd=ROOT)

actual = [int(line) for line in (ROOT / "rtl_bits.txt").read_text().splitlines() if line.strip()]
expected = encode_bits([100, 101, 99, 99, 104, 100, 100, 101], k=2)
if actual != expected:
    raise SystemExit(f"FAIL: RTL/golden mismatch: actual={actual} expected={expected}")
print(f"PASS: {len(actual)} bits match the Python golden model")
