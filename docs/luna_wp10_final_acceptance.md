# WP10 final acceptance and closure review

Date: 2026-07-28  
Reviewer: Luna  
Evaluated revision: `07d14cc` (`Complete WP10.8 OPEN_FFT performance and fp32 gate`)

## Decision

**GO — WP10.9 is complete for the accepted CUDA scope.**

The documentation audit is consistent. The device-independent gates and the
serial CUDA rerun pass for fp64 and the explicitly accepted CUDA fp32 scope.
HIP is absent and is not used as a gate, per the author update. No
implementation repair was made.

## Reviewed authority and commit chain

The review read the WP10 blueprint, WP10 acceptance/performance reports,
`FFT-dipole_implementation_plan.md`, `GPU_FFT_DIPOLE_DESIGN.md`, the final
Fortran input parser and startup/validation path, and
`tests/dipole_validation/README.md`.

Relevant implementation/acceptance commits are:

| Commit | Scope |
|---|---|
| `9a5ec4f` | independent OPEN_FFT host kernel |
| `1224d60` | CUDA OPEN_FFT production mode |
| `0b2ac56` | full-run GPU memory accounting |
| `0830897` | basis-resolved block-one derivative stabilization |
| `8d101f0` | finite uniform coarse projection |
| `52a5446` | coarse-block production integration/acceptance |
| `07d14cc` | WP10.8 performance and CUDA fp32 decision |

The final acceptance entry is checked in the blueprint. HIP parity, CPU repair,
and 2D periodicity remain separately gated work.

## Current commands and results

Build and host-side checks:

```text
cmake --build build_gpu_wp10_cuda -j2                                      PASS
cmake --build build_gpu_fp32 -j2                                           PASS
python3 tests/dipole_validation/test_open_fft_oracle.py -v                 PASS, 9/9
python3 tests/dipole_validation/test_periodic_ewald_reference.py -v        PASS
python3 tests/dipole_validation/test_open_kernel_goldens.py \
  --driver build_gpu_wp10_cuda/bin/dipole_open_kernel_host_driver             PASS
ctest --test-dir build_gpu_wp10_cuda --output-on-failure \
  -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-ewald-host-builder'  PASS, 4/4
ctest --test-dir build_gpu_fp32 --output-on-failure \
  -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-ewald-host-builder'  PASS, 4/4
ctest --test-dir build_gpu_wp10_cuda --output-on-failure \
  -R 'regression-test|asd-tests'                                           PASS, 2/2
python3 tests/dipole_validation/run_wp5_e2e.py \
  --binary build_gpu_wp10_cuda/bin/sd.f95.cuda                              PASS;
  max_field=4.3719150199989277e-40 T, max_energy=3.8376175788916486e-43
python3 tests/dipole_validation/run_wp10_open_coarse_e2e.py \
  --binary build_gpu_wp10_cuda/bin/sd.f95.cuda                              PASS;
  max_field=4.999999675228202e-39 T, max_energy=0
python3 tests/dipole_validation/run_wp5_e2e.py \
  --binary build_gpu_fp32/bin/sd.f95.cuda --fp32                            PASS;
  max_field=3.0000000013187062e-36 T, max_energy=2.8999999998223618e-39
python3 tests/dipole_validation/run_wp10_open_coarse_e2e.py \
  --binary build_gpu_fp32/bin/sd.f95.cuda --fp32                            PASS;
  max_field=4.5500000015329382e-36 T, max_energy=1.0999999992719863e-41
```

The current isolated CUDA runtime probes were:

```text
ctest --test-dir build_gpu_wp10_cuda --output-on-failure                 PASS, 11/11
ctest --test-dir build_gpu_fp32 --output-on-failure \
  -R 'dipole-gpu-fft-convolution|dipole-open-fft-layout'                 PASS, 2/2
compute-sanitizer --tool memcheck --error-exitcode=99 \
  build_gpu_wp10_cuda/bin/dipole_open_fft_layout_tests                    PASS, ERROR SUMMARY: 0 errors
compute-sanitizer --tool memcheck --error-exitcode=99 \
  build_gpu_fp32/bin/dipole_open_fft_layout_tests                         PASS, ERROR SUMMARY: 0 errors
```

Hardware was two NVIDIA RTX A4000 GPUs, driver/KMD `610.43.02`, CUDA UMD
`13.3`, with GPU 0 showing 2540 MiB in use and GPU 1 showing 15 MiB in use at
the start of the rerun. The toolkit is `/usr/local/cuda-13.3/bin/nvcc`
(`13.3.73`); Compute Sanitizer is `2026.2.1.0`. `hipcc`, `rocminfo`,
`/dev/kfd`, and an AMD device are absent.

## Numerical maxima

The following current-run maxima are reproducible:

| Check | Maximum |
|---|---:|
| Frozen OPEN_FFT host goldens, field | `5.3290705182007514e-15` |
| Host direct finite convolution | `7.815970093361102e-14` |
| Host dimensionless energy identity | `3.5527136788005009e-15` |
| Host reciprocity | `7.1054273576010019e-15` |
| Host projected coarse tensor, skew case | `9.770e-15` |
| Host block-one nonuniform convergence field | `6.661e-16` |

| CUDA check | Maximum |
|---|---:|
| CUDA fp64 production E2E field | `4.3719150199989277e-40 T` |
| CUDA fp64 production E2E energy | `3.8376175788916486e-43 mRy/atom` |
| CUDA fp64 OPEN_FFT coarse E2E field | `4.999999675228202e-39 T` |
| CUDA fp64 OPEN_FFT coarse E2E energy | `0` |
| CUDA layout impulse field | `2.2737367544323206e-13` |
| CUDA layout impulse energy | `2.7284841053187847e-12` |
| CUDA coarse production field | `9.9475983006414026e-14` |
| CUDA coarse production energy | `3.4106051316484809e-13` |
| CUDA fp32 independent field | `3.20507819822069e-05` |
| CUDA fp32 independent energy | `3.8058345853642095e-04` absolute; `1.442401007589701e-07` relative |
| CUDA fp32 coarse production field | `2.641373406220282e-06` |
| CUDA fp32 coarse production energy | `2.4842694870130799e-06` absolute; `6.0712226246373291e-08` relative |

For traceability, these are consistent with the prior accepted CUDA reports:

| Accepted CUDA scope | Field maximum | Energy maximum |
|---|---:|---:|
| OPEN_FFT production matrix, fp64 formatted output | `4.602701e-39 T` | `7.192777e-41 mRy/atom` |
| Basis-resolved block one, NA=2, internal | `4.951753e-38 T` | `8.407791e-45` |
| Coarse production projection | `9.9475983006414026e-14` | `3.4106051316484809e-13` |
| CUDA layout production seam | `2.2737367544323206e-13` | `2.7570656868647347e-10` |
| CUDA fp32 independent seam | `3.20507819822069e-05` | `3.8058345853642095e-04` absolute; `1.442401007589701e-07` relative |

The accepted/current fp32 derivative maximum was `1.5934170747300413e-05`, within
the documented `5e-5` mixed budget. It remains accepted for CUDA only; HIP
fp32 remains explicitly rejected pending its own oracle, E2E, and sanitizer
matrix.

## Capability and rejection matrix

| Dimension | Accepted/implemented scope | Rejected or separately gated |
|---|---|---|
| Physics selector | `OPEN_FFT`: finite/open, zero-padded convolution | `EWALD3D_FFT`: separate 3D-periodic tin-foil solver; OPEN_FFT is not its fallback |
| Atomic baseline | Only `block=(1,1,1)`, `NA=1` | `NA>1`, block one is basis-resolved block one, not atomic baseline |
| Coarse representation | Existing basis-resolved macrocell channels; positive uniform blocks that divide every `N` axis; exact finite projection | Partial edges, unequal populations, irregular/deleted/occupancy layouts |
| Boundary conditions | OPEN_FFT requires `BC 0 0 0`; EWALD3D_FFT requires `BC P P P` | Mixed BC and the other mode’s BC contract |
| Dynamics | GPU spin dynamics/SD only; legacy `do_dip=0` | MC and GPU MC; no double application with legacy dipoles |
| Precision | CUDA fp64; CUDA fp32 under the explicit WP10.8 mixed budget | HIP fp32; HIP runtime parity is unavailable and non-gating here |
| FFT layout | Uniform active grid embedded in legal `P`, with `P_i >= 2G_i-1` and cleared padding | Partial/irregular layouts and unproven atom ordering |

The parser accepts `gpu_dipole_mode`, surface, tolerance, alpha/rcut, and
mesh keys. Defaults are `OFF`, `TINFOIL`, zero alpha/rcut, zero mesh, and
`gpu_dipole_tol=1e-10`. Startup validation rejects unsupported mode names,
`PME3D`, wrong BC, nonzero `do_dip`, MC, non-default surface/overrides, and
non-divisible blocks before the GPU dipole allocation path. Startup messages
distinguish `OPEN_FFT` geometry/allocation from `EWALD3D_FFT` geometry/memory.

## Memory and performance evidence

The accepted CUDA accounting report used the same full-run scope for the
preflight and live tracker:

| Input | Full-run projected B | OPEN_FFT phase peak B | Tracker peak B | Release B |
|---|---:|---:|---:|---:|
| `2x1x1`, M=1 | 25,056 | 1,040 | 25,056 | 0 |
| `2x2x1`, M=4 | 34,400 | 6,080 | 34,400 | 0 |
| `4x3x1`, M=1 | 37,336 | 9,440 | 37,336 | 0 |

The serial CUDA device sweep (steady-state microseconds per evaluation) was:

| Shape | FFT grid P | Builder ms | Setup ms | Total us | Persistent B | Construction B | Workspace B | Tracker peak B |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| thin `128x64x1` | `255x127x1` | 2.507 | 7.300 | 252.262 | 5,456,072 | 2,331,720 | 5,182,544 | 13,592,984 |
| thin `128x64x1` | `256x128x2` | 3.676 | 9.533 | 251.095 | 11,071,640 | 4,718,592 | 0 | 16,412,880 |
| thin `128x64x1` | `256x128x1` | 2.509 | 5.811 | 148.479 | 5,535,896 | 2,359,296 | 0 | 8,517,840 |
| racetrack `256x16x1` | `511x31x1` | 1.298 | 5.292 | 165.251 | 2,665,160 | 1,140,552 | 2,537,040 | 6,654,104 |
| racetrack `256x16x1` | `512x32x2` | 1.843 | 6.014 | 138.752 | 5,520,536 | 2,359,296 | 0 | 8,191,184 |
| racetrack `256x16x1` | `512x32x1` | 1.228 | 4.010 | 85.391 | 2,760,344 | 1,179,648 | 0 | 4,251,344 |
| 3D `32x24x16` | `63x47x31` | 7.043 | 12.608 | 441.969 | 15,595,880 | 6,608,952 | 0 | 23,138,776 |
| 3D `32x24x16` | `64x48x32` | 7.344 | 12.987 | 307.364 | 16,883,864 | 7,077,888 | 0 | 24,895,696 |
| 3D `32x24x16` | `64x48x36` | 7.655 | 14.223 | 347.556 | 18,994,328 | 7,962,624 | 0 | 27,890,896 |

The phase timings are recorded in the command output and the complete
benchmark definition remains in `docs/terra_wp10_8_open_performance.md`. No
padding policy, cache, singleton-axis specialization, or reduced construction
allocation was enabled. CUDA memcheck reported `ERROR SUMMARY: 0 errors` for
the fp64 and fp32 seams.

## Representation and deferred-work audit

The required documentation is consistent on the central representation:
the existing basis-resolved macrocell map is the sole runtime source
representation. `NA=1` with unit blocks is the atomic baseline; `NA>1` at
unit blocks is explicitly basis-resolved; coarse blocks are the projected
macrocell Hamiltonian. No required document claims a second atomistic GPU
backend or a hidden atom-array FFT path.

WP11 remains independently open: it covers repair and separate acceptance of
legacy CPU `do_dip=1/2/3`, including known non-cubic FFTW defects. WP12 remains
independently open: it requires a true 2D-periodic Green function/layer
correction, zero-mode and surface derivation, and its own CUDA/HIP oracle and
acceptance. The proposed `EWALD2D_FFT` name is not enabled by WP10.

## Closure

No required WP10 work remains for the accepted CUDA scope. The CUDA fp64 and
fp32 matrices, production E2E checks, rejection paths, memory accounting,
performance sweep, and sanitizer checks pass. HIP numerical/fp32/sanitizer
parity remains a documented non-gating limitation. WP11 CPU repair and WP12
2D periodicity remain independently open and are not implied complete.
