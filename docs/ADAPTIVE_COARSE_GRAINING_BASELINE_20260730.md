# Adaptive coarse-graining pre-remediation baseline

Date: 2026-07-30  
Commit: `f8e27f1eb6b95ed114d8a0bc8801686f6c529203`  
Commit subject: `Complete adaptive CG DMI and anisotropy production support`

## Reproducibility state

Before this task, `git status --short` reported one modified tracked file and
97 untracked paths (98 changed paths total). The modified file was
`docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`. Existing build trees and local
fixtures were not used as source evidence.

The fresh CPU command was:

```sh
cmake -S . -B /tmp/uppasd-cg-baseline-cpu-20260730 \
  -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF \
  -DUPPASD_PRECISION=DOUBLE -DBUILD_TESTING=ON -DUSE_OPENMP=ON \
  -DUSE_MKL=OFF -DUSE_FFTW=OFF
cmake --build /tmp/uppasd-cg-baseline-cpu-20260730 --parallel 2
```

Configuration succeeded with GNU Fortran 13.3.0 (`/bin/f95`), GNU C 12.4.0,
OpenMP, generic BLAS/LAPACK, backend `OFF`, and precision `DOUBLE`. The build
failed at 70% while compiling `source/MonteCarlo/montecarlo.f90:58`:

```text
Fatal Error: Cannot open module file ‘adaptivetimestepping.mod’ for reading at (1): No such file or directory
gmake[2]: *** [CMakeFiles/asdlib.dir/build.make:1232: CMakeFiles/asdlib.dir/source/MonteCarlo/montecarlo.f90.o] Error 1
gmake[1]: *** [CMakeFiles/Makefile2:536: CMakeFiles/asdlib.dir/all] Error 2
gmake: *** [Makefile:146: all] Error 2
```

The complete raw build log is retained at
`/tmp/uppasd-cg-baseline-cpu-20260730/build.log`; CTest discovery is at
`/tmp/uppasd-cg-baseline-cpu-20260730/ctest-list.log`. This is a failed
pre-fix baseline, not a release or parity run. `HEAD` contains the live import
at `source/MonteCarlo/montecarlo.f90:58`, but no tracked
`source/CoarseGraining/adaptivetimestepping.f90` and no CMake source entry for
that module.

## Adaptive-CG build and test inventory

| Target or test | Executable/role | Backend and precision | CTest label or benchmark status |
| --- | --- | --- | --- |
| `asdlib` + `sd` | Production library and ordinary executable; adaptive production sources include `adaptivecgproduction.f90` and `adaptivehybridsolver.f90` | CPU `OFF`/DOUBLE; CUDA or HIP; GPU DOUBLE or SINGLE | `adaptive-cg-production-e2e`: `cg13-cpu`, `cg13-cuda`, or `cg13-hip` |
| `block_topology_tests`, `stiffness_material_tests`, `coarse_tensor_operator_tests`, `multichannel_coarse_tensor_operator_tests`, `smooth_projected_operator_tests`, `static_hybrid_operator_tests`, `block_selector_tests`, `adaptive_hybrid_solver_tests` | Eight focused reference executables | CPU reference; retained in GPU configurations | `coarse-graining;cg13;cg13-cpu;reference`; GPU configurations append `cg13-cuda` or `cg13-hip` |
| `adaptive-cg-setup-rejection-matrix` | `run_setup_rejection_matrix.py` against ordinary `sd` | CPU or selected GPU-configured ordinary executable | `coarse-graining;cg13;cg13-cpu;rejection` or matching backend label |
| `gpu_adaptive_runtime_tests` | Backend runtime parity/ownership executable | CUDA or HIP; DOUBLE by default, SINGLE when configured | `coarse-graining;cg13;cg13-cuda;reference` or `cg13-hip;reference` |
| `gpu_adaptive_runtime_benchmark` | Hardware adaptive-runtime benchmark | CUDA or HIP; DOUBLE or SINGLE | Not a CTest; reports feature-off overhead and active-DOF crossover |
| `adaptive-cg-fixture-dependencies` | Tracked-source/package audit | Backend-independent | `coarse-graining;cg13;packaging` |

Preset coverage is CPU `cpu`/`cpu-debug` (DOUBLE, backend OFF), CUDA
`cuda`, `cuda-debug`, `cuda-relwithdebinfo`, `cuda-fastcopy`, and `cuda-ptds`
(DOUBLE), CUDA `cuda-single` (SINGLE), and HIP `hip-cdna`, `hip-rdna`, and
`lumi` (DOUBLE defaults). `MIXED` is rejected by CMake as not implemented.
The CPU CTest discovery contained the eight reference tests, the two ordinary
executable harnesses, and no GPU runtime target.

The tracked legacy paths are `source/Makefile.gfortran`,
`source/Makefile.legacy`, and `source/Makefile.legacy.gfortran`. Their source
lists do not include the adaptive-CG production package or repair the dangling
`AdaptiveTimeStepping` import, so they are inventoried as unsupported for this
pre-remediation adaptive-CG baseline. The CMake CPU/CUDA/HIP paths are the
paths represented by the target and label matrix above.

## Executable fixture inventory

The production harnesses reference these 32 e2e directories:

```text
adaptive_mixed, feature_off, gpu_adaptive_mixed, gpu_fft_static_mixed,
gpu_static_mixed, initial_phase_g, initial_phase_heat_bath, initial_phase_mc,
initial_phase_mc_gpu, initial_phase_q, initial_phase_q_gpu, initial_phase_sd_cpu,
initial_phase_sd_gpu, initial_phase_y, initial_phase_z, initmag_spin_spiral,
invalid_mask, invalid_partial_block, missing_alat, parity_adaptive_cpu,
parity_adaptive_gpu, parity_static_cpu, parity_static_gpu, static_all_coarse,
static_all_fine, static_mixed, unsupported_initial_phase_x,
unsupported_temperature, anisotropy_uniform_coarse, anisotropy_uniform_fine,
dmi_anisotropy_mixed, dmi_anisotropy_mixed_gpu
```

The production harness also launches
`examples/AdaptiveCoarseGraining/{static_mixed,adaptive,initial_phase_texture}`;
the audit therefore checks 35 fixture directories in total.
`audit_fixture_dependencies.py` resolves every harness case and every input
file named by its `inpsd.dat` against `git ls-files`; generated restart/output
files are not source dependencies.

The audit result on this baseline is:

```text
adaptive-CG fixture dependency audit: FAIL
checked 35 fixture directories and 57 input paths
```

It reports the five local-only e2e directories `missing_alat`,
`anisotropy_uniform_fine`, `anisotropy_uniform_coarse`,
`dmi_anisotropy_mixed`, and `dmi_anisotropy_mixed_gpu`, plus their missing
`kfile_cg`, `kfile_cg_x`, and `dmfile_cg` dependencies. The audit is a
machine-checkable packaging failure; RCG-00 does not add the fixtures.

## Stale incremental artifacts

The existing incremental trees contain 21 locations with
`adaptivetimestepping.f90.o` and its CMake stamp, and 19 retained
`adaptivetimestepping.mod` files. For example:

```text
build_cpu/mod/adaptivetimestepping.mod
build_cpu/CMakeFiles/asdlib.dir/source/CoarseGraining/adaptivetimestepping.f90.o
```

These artifacts conceal the defect because an incremental build's dependency
scanner sees a valid module in its configured module search path and an
up-to-date object/stamp, even though the tracked source and CMake graph no
longer provide that module. The stale module can satisfy `use
AdaptiveTimeStepping` and the stale object can remain linkable, so only a new
out-of-tree build removes the false dependency.

## Claim freeze

The existing CUDA RTX A4000 fp64 `gpu_adaptive_runtime_benchmark` claim is a
synthetic adaptive-runtime measurement: reported active-DOF crossover
`1.3016x` at ratio `0.813232` after kernel optimization. It is not a
comparison with UppASD's production atomistic GPU path.

For remediation purposes the following claims are **SUSPENDED**:

- DMI production support and chirality claims;
- CPU/GPU end-to-end parity claims, including DMI/anisotropy parity;
- all GPU speedup and crossover claims.

They remain suspended until the clean build, tracked-fixture, accepted-DMI,
moving-state parity, and production-baseline benchmark gates close. No
physics or kernel optimization was performed in RCG-00.
