# WP10.8 Terra performance report — OPEN_FFT

Date: 2026-07-28  
Revision: working tree based on `52a5446`

## Scope and method

`OPEN_FFT` remains the accepted finite, zero-padded point/projection operator:
the active macro grid is embedded at the origin of a cleared FFT grid and the
complete finite kernel is embedded for the *selected* legal grid `P`.  This
report does not compare it to periodic Ewald physics.

The benchmark target now has an OPEN_FFT path.  It constructs the production
host fp64 kernel, accepts an explicit `--fft-grid P1 P2 P3`, uploads it through
`uploadOpenKernel()`, then measures the full steady-state path:

```text
pack/clear -> forward FFT -> spectral contraction -> inverse FFT
           -> atom scatter -> active-macro energy reduction
```

It prints a total CUDA/HIP event time per evaluation, event-backed per-phase
times from the existing production stopwatch, host builder time, plan/upload
setup time, persistent/construction/workspace estimates, and the live tensor
allocation tracker current/peak.  Timers are disabled during warmup, so every
printed phase is divided only by the measured iteration count.  The harness
uses the same basis-resolved active macro layout as production and does not
substitute a periodic or truncated kernel.

## Required device sweep

The following commands use `Pmin=2G-1`, `Peven=2G`, and an additional legal
candidate that leaves singleton axes unpadded where that is legal.  They are
deliberately separate rows: no conclusion may assume that a nearby FFT extent
is faster.

```bash
cmake --build build_gpu_wp10_cuda --target dipole_gpu_fft_benchmark -j2

# Thin film: G=128x64x1
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 128 64 1 --fft-grid 255 127 1 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 128 64 1 --fft-grid 256 128 2 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 128 64 1 --fft-grid 256 128 1 --warmup 10 --iterations 100

# Long racetrack-like sample: G=256x16x1
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 256 16 1 --fft-grid 511 31 1 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 256 16 1 --fft-grid 512 32 2 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 256 16 1 --fft-grid 512 32 1 --warmup 10 --iterations 100

# Fully 3D sample: G=32x24x16
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 32 24 16 --fft-grid 63 47 31 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 32 24 16 --fft-grid 64 48 32 --warmup 10 --iterations 100
build_gpu_wp10_cuda/bin/dipole_gpu_fft_benchmark --open --grid 32 24 16 --fft-grid 64 48 36 --warmup 10 --iterations 100

ctest --test-dir build_gpu_wp10_cuda --output-on-failure \
  -R 'dipole-open-fft-layout|dipole-open-fft-coarse-e2e'
compute-sanitizer --tool memcheck --error-exitcode=99 \
  build_gpu_wp10_cuda/bin/dipole_open_fft_layout_tests
```

Run the corresponding fp32 build only as an error sweep; its production
OPEN_FFT selector must remain rejected unless the independent-oracle,
field/energy, and sanitizer gates are explicitly accepted.  Run the same
matrix on HIP with the same independent oracle before enabling an optimization
that relies on a CUDA/HIP FFT-plan property.

## CUDA results

The sweep was run on 2026-07-28 using CUDA 13.3.73 and NVIDIA driver
610.43.02.  Two RTX A4000 GPUs were visible; at the pre-run snapshot GPU 0
also held a 2.398 GiB Python allocation, so these are one-sample engineering
measurements, not a noise-characterized throughput claim.  Every benchmark
row used the production finite builder and `basis=ensembles=1` in fp64.

All times below are microseconds per steady-state evaluation. `builder` is
host-fp64 finite-tensor construction; `setup` is plan creation plus kernel
upload/transform.  The six phase columns are device-event times. `persist`,
`construct`, `workspace`, `current`, and `peak` are bytes from the estimator
and live tensor tracker.

| Shape | P | builder ms | setup ms | total us | pack/clear | forward | contract | inverse | scatter | energy | persist | construct | workspace | current | peak |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| thin 128x64x1 | 255x127x1 | 2.489 | 222.649 | 252.180 | 8.388 | 103.091 | 12.914 | 103.036 | 7.034 | 7.582 | 5,456,072 | 2,331,720 | 5,182,544 | 11,261,264 | 13,592,984 |
| thin 128x64x1 | 256x128x2 | 3.644 | 71.934 | 251.320 | 9.071 | 96.841 | 21.916 | 98.085 | 7.180 | 7.637 | 11,071,640 | 4,718,592 | 0 | 11,694,288 | 16,412,880 |
| thin 128x64x1 | 256x128x1 | 2.465 | 5.740 | 148.695 | 8.196 | 50.308 | 14.009 | 50.841 | 7.061 | 7.516 | 5,535,896 | 2,359,296 | 0 | 6,158,544 | 8,517,840 |
| racetrack 256x16x1 | 511x31x1 | 1.238 | 19.583 | 164.690 | 6.659 | 64.138 | 7.796 | 63.923 | 5.548 | 6.152 | 2,665,160 | 1,140,552 | 2,537,040 | 5,513,552 | 6,654,104 |
| racetrack 256x16x1 | 512x32x2 | 1.880 | 36.737 | 138.691 | 7.659 | 47.283 | 13.296 | 47.396 | 6.190 | 6.430 | 5,520,536 | 2,359,296 | 0 | 5,831,888 | 8,191,184 |
| racetrack 256x16x1 | 512x32x1 | 1.242 | 4.004 | 85.125 | 6.186 | 24.576 | 8.100 | 24.819 | 4.851 | 5.969 | 2,760,344 | 1,179,648 | 0 | 3,071,696 | 4,251,344 |
| 3D 32x24x16 | 63x47x31 | 6.996 | 27.040 | 442.151 | 11.895 | 188.595 | 35.656 | 178.716 | 8.031 | 8.948 | 15,595,880 | 6,608,952 | 0 | 16,529,824 | 23,138,776 |
| 3D 32x24x16 | 64x48x32 | 7.332 | 43.118 | 307.537 | 12.582 | 112.156 | 40.221 | 115.466 | 7.527 | 9.022 | 16,883,864 | 7,077,888 | 0 | 17,817,808 | 24,895,696 |
| 3D 32x24x16 | 64x48x36 | 7.676 | 32.365 | 347.625 | 13.834 | 126.525 | 47.657 | 132.421 | 7.536 | 8.830 | 18,994,328 | 7,962,624 | 0 | 19,928,272 | 27,890,896 |

The faster CUDA candidates were legal finite embeddings, and the existing
independent host suite already checks oversized-padding zero gaps and compares
the builder/convolution with direct finite point sums.  No new expected values
were derived from GPU output.

| Comparison | CUDA steady-state result |
|---|---:|
| Thin: 256x128x1 versus exact minimum | 41.0% faster; 37.3% lower tracker peak |
| Racetrack: 512x32x1 versus exact minimum | 48.3% faster; 36.1% lower tracker peak |
| 3D: 64x48x32 versus exact minimum | 30.4% faster; 7.6% higher tracker peak |

Correctness and allocation/safety checks on the CUDA host passed:

| Check | Result |
|---|---:|
| CUDA build: benchmark and layout executable | PASS |
| `ctest -R 'dipole-open-fft-layout|dipole-open-fft-coarse-e2e'` | PASS, 2/2 |
| Production field maximum | `2.2737367544323206e-13` |
| Production energy maximum | `2.7570656868647347e-10` |
| Coarse field / energy maximum | `9.9475983006414026e-14` / `3.4106051316484809e-13` |
| CUDA memcheck | PASS, `ERROR SUMMARY: 0 errors` |
| HIP build/launch | BLOCKED: no HIP toolchain/device |

## CUDA fp32 result

The fp32 device-storage seam was rerun against the independent finite
point/projection evaluator.  Its direct maxima were field
`3.20507819822069e-05`, energy absolute `3.8058345853642095e-04`, energy
relative `1.442401007589701e-07`, and field/energy derivative residual
`1.5934170747300413e-05`.  This meets the established WP9 mixed fp32 budget:
field and derivative are below `5e-5`; energy is accepted when either its
absolute error or its error relative to `max(1, |E|)` is below `5e-5`.

The actual fp32 `sd.f95.cuda` OPEN_FFT coarse E2E writer run also passed with
`max_field_tesla_error=4.5500000015329382e-36` and
`max_energy_mry_error=1.0999999992719863e-41`.  CUDA memcheck of the fp32
layout/production seam reported `ERROR SUMMARY: 0 errors`; that seam checks
the estimator against its allocation inventory.

## Explicit optimization decisions

| Item | Decision | Reason |
|---|---|---|
| Exact-minimum padding `P=2G-1` | Keep as the current default | It remains the portable accepted embedding. |
| `P=2G` or other backend-friendly padding | **Reject / not enabled** | CUDA has strong shape-dependent candidates, but HIP parity and repeated-run data are required before selecting a backend policy. |
| Immutable OPEN_FFT host-kernel/spectrum cache | **Reject / not enabled** | Builder time is only 1.24–7.68 ms in this sweep; no repeated-equivalent-initialization measurement justifies retained host memory. Any later key must include model/artifact revision, host/device storage format, `G`, selected `P`, block shape, exact primitive-vector bits, ordered basis offsets, and basis count. |
| Singleton-axis FFT-plan specialization | **Reject / not enabled** | No CUDA+HIP support/accuracy evidence. |
| Reduced construction allocation | **Reject / not enabled** | Current accounting remains the accepted construction inventory; no measured allocation benefit. |
| fp32 device storage/execution | **Enable for CUDA; reject for HIP** | CUDA passes the independent finite-oracle seam, fp32 production E2E, and memcheck under the `5e-5` mixed budget. HIP fp32 is explicitly rejected at preflight pending the identical oracle, E2E, and sanitizer matrix. |

No periodic images, real-space truncation, stale-padding shortcut, altered
coarse projection, fallback, or Hamiltonian change was introduced by this
work.  The next update must append repeated CUDA samples, the same HIP sweep,
and the HIP fp32 independent-oracle, production-E2E, and sanitizer matrix
before revisiting the HIP-specific rejection.

## WP10.8 completion and excluded work

**Decision: complete for the measured CUDA scope.** CUDA fp64 and fp32
OPEN_FFT have an explicit performance/precision record, production-E2E
coverage, allocation accounting, and sanitizer evidence.  This closure does
not claim cross-backend acceptance or silently enable an unmeasured shortcut.

The following are intentionally left out of the completed CUDA scope:

- HIP numerical, fp32, and sanitizer parity: no HIP compiler or device was
  available. HIP fp32 OPEN_FFT is rejected at preflight.
- A shape-specific padding policy: CUDA candidates are measured, but repeated
  samples and HIP results are required before selecting a portable policy.
- Immutable host-kernel/spectrum cache: measured builder time did not justify
  retained host memory, and no repeated-initialization gain was established.
- Singleton-axis FFT-plan specialization and reduced construction allocation:
  no CUDA/HIP safety or measured-benefit case supports enabling either.

These are deliberate non-enablement decisions, not changes to the finite
Hamiltonian. They may be reconsidered only in a separately measured follow-up.
