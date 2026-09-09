# Zynq-7020 RTL compression demo

This is the first hardware-learning milestone for the project. It is deliberately
small and deterministic:

```text
signed 16-bit samples -> first-order delta -> 17-bit zigzag -> fixed-K Rice bits
```

`rice_demo_encoder.v` captures one frame and emits one compressed bit per clock.
The first sample is sent as 16 bits, MSB first. Each later delta is encoded as
`0...01` followed by the low `K` bits. The default frame has 16 samples and
`K=2`.

The Python reference is `golden_model.py`. Run the self-check with:

```powershell
python -m pytest rtl_demo/test_golden.py
```

The RTL has no Vivado-specific primitives, so it can later be simulated by
Vivado XSim, Icarus, or Verilator. The next milestone is a testbench that feeds
the fixed vector and packs `out_bit` into bytes for comparison with the Python
model.

The testbench is included. Once Icarus Verilog is installed, run:

```powershell
python rtl_demo/run_sim.py
```

It compiles the RTL, runs the testbench, and compares every emitted bit with
the Python golden model.
