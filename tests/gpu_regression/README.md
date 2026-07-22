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

## Convolution benchmark

`bench.py` compares GPU real-space and device convolution Hamiltonians. Its
simple-cubic RKKY sweeps cover 16³, 40³, and 100³ cells with 5,000 steps each
at fixed `rmax=5` (514 neighbours per atom), plus an `rmax=3, 5, 8` sweep at
100³ cells (122, 514, and 2,108 neighbours per atom). The `rmax=5`, 100³ case
is shared by the two sweeps. It also includes a BCC matrix, a two-atom-basis
FeCo case, and the 131,072-atom dhcp-Nd long-range fixture (four-site basis,
1,340 neighbours/atom). The dhcp-Nd case uses 1,000 steps by default to keep
routine benchmark time reasonable. Every configuration has one
warm-up plus at least three timed runs of each backend; the convolution run is
rejected if the executable did not print its activation line (preventing an
unnoticed sparse fallback).  The harness forces `do_reduced Y`: convolution
operates on the unit-cell Hamiltonian (`nHam == NA`), while the FFT supplies the
full periodic lattice field.  It prints the convolution activation line once
per benchmark; if setup falls back, the reported error includes UppASD's
reason. It enables `do_gpu_timings Y` and reports both the GPU measurement-phase
total and the GPU Hamiltonian time. Steps/s and both convolution speedups are
derived from those GPU-event times; the process and Fortran measurement-phase
wall times are retained as supplementary data.

The simple-cubic cases are based on `tests/SC_conv`. For every isolated run,
the harness invokes that fixture's `rkky_gen.py` to regenerate `jij` with the
requested cutoff before measuring the sparse/convolution crossover. Select one
explicitly, for example `--case sc_rkky_100_r8`.

```bash
python3 tests/gpu_regression/bench.py --list
python3 tests/gpu_regression/bench.py --binary ./bin/sd.f95.cuda \
  --output /tmp/uppasd-convolution-bench.json --csv /tmp/uppasd-convolution-bench.csv
```

Use `--case sc_rkky_40 --steps 5000` to narrow or scale a run.  The terminal
report gives median wall time, steps/s, and convolution speedup; JSON retains
the full environment and summary, while CSV contains every timed sample.
