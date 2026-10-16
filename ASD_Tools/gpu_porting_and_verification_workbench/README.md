<!--
SPDX-FileCopyrightText: 2026 CSC – IT Center for Science

SPDX-License-Identifier: Apache-2.0
-->

# GPU porting and verification workbench for UppASD

A workbench for porting and optimizing the [UppASD](https://github.com/UppASD/UppASD) multiscale module for GPUs.

The goal is to simplify benchmarking, verification, and development of new GPU-ported functions by providing a minimal, isolated environment with few moving parts before copying the solutions to UppASD.

## Limitations

- Currently verification is done by a common-sense comparison of floating point numbers and accepting relative/absolute errors below certain units in the last place (ULPs). The code should be rechecked that it is sufficiently robust also for production runs.
- The initial condition is inspired by UppASD setup and may not reflect actual values used for production runs. Should as well be validated before using for production.
- Currently the Fortran model solution is copied over from UppASD. Ideally one would interface UppASD and use the native function directly to eliminate translation errors in this aspect.


## Building

Module setup on LUMI:

```bash
module load PrgEnv-gnu craype-accel-amd-gfx90a
module load LUMI/25.03 partition/G rocm/6.3.4-extras
```

```bash
cmake -S . -B build
cmake --build build --parallel
```

## Running

```bash
build/src/test-host   [--atoms N] [--ensembles N] [--no-verify]
build/src/test-device [--atoms N] [--ensembles N] [--no-verify]
```

## File layout

| File | Role |
| --- | --- |
| `src/model.f90` | Fortran reference implementation |
| `src/host-implementation.c` | Serial C implementation |
| `src/device-implementation.cu` | Parallel GPU implementation |
| `src/testcase.c` | Initial condition setup |
| `src/testcase-device.c` | Device allocations and setup |
| `src/testcase-verify.c` | Verification |
| `src/test-host.c`, `src/test-device.c` | Main functions for `test-host` and `test-device` programs |


## Profiling on LUMI

Run on a GPU node, e.g. via `srun --account=<project> --partition=dev-g --gpus-per-node=1 --ntasks=1 --time=00:10:00`.

Modern `rocprof` (>= 6.3)

```bash
# Generating profiles (compute and memory unit utilization, occupancy, etc)
rocprof-compute profile --name test -- build/src/test-device --atoms 100000 --ensembles 16 --no-verify
rocprof-compute analyze --path workloads/test/MI200
```

```bash
# Generating traces
rocprof-sys-run -- build/src/test-device --atoms 100000 --ensembles 16 --no-verify
```

Older `rocprof`.

```bash
# Kernel statistics (*_kernel_stats.csv)
rocprofv3 --stats --kernel-trace --truncate-kernels --output-directory profile -- build/src/test-device --no-verify

# Timeline (*.pftrace). Open in https://ui.perfetto.dev
rocprofv3 --hip-trace --kernel-trace --memory-copy-trace --output-format pftrace --output-directory timeline -- build/src/test-device --no-verify
```

## License

Code derived from UppASD is GPL-3.0-or-later. Other code is Apache-2.0 or GPL-3.0-or-later, as stated in each file header. License texts are in `LICENSES/`. The headers follow [REUSE](https://reuse.software).
