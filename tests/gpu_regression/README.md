# GPU regression harness

Run the same T=0 input through the GPU-enabled executable twice: once with
`do_gpu N` (Fortran reference path) and once with `do_gpu Y`. The harness copies
each fixture to separate temporary directories, compares final restart moments and
all rows of `totenergy.*.out`, and prints a compact pass/fail table.

```bash
python3 tests/gpu_regression/run.py --binary ./bin/sd.f95.cuda
```

Useful variants:

```bash
python3 tests/gpu_regression/run.py --list
python3 tests/gpu_regression/run.py --binary ./bin/sd.f95.cuda --case bcc_tail_block_1000 --keep
```

The manifest covers scalar exchange, DM plus uniaxial anisotropy, tensor exchange,
reduced Hamiltonians, four ensembles, a 1000-site tail block, an odd 3·N·M count,
the convolution backend, GPU-side energy measurement for scalar and convolution
paths, GPU Monte Carlo, and a single-spin easy-axis uniaxial MC invariant. Cases
normally force `do_gpu_measurements N`; the two `*_gpu_measurements` cases override
only the GPU run to `Y` and compare its device-generated energy output to the
Fortran reference. It must be executed on a host with a
visible CUDA/HIP device; compile-only CI intentionally does not run it.
