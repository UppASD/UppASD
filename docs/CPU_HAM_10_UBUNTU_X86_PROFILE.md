# CPU-HAM-10 — Ubuntu/x86 provider and hardware-counter profile

**Date:** 2026-09-02
**Source commit:** `fe72a43c7f7dcbf9abaa2d5f189010a56faf4c78`
**Status:** measured campaign complete with provider/counter limitations recorded

## Executive result

The current Ubuntu host confirms the workload crossover seen by CPU-HAM-05:
long-range Nd strongly favours CONVOLUTION, while the short-range control
favours DIRECT. Portable SPARSE did not beat DIRECT robustly in the measured
production controls. The available evidence does not prove whether Nd DIRECT
is DRAM-bandwidth-limited, cache/gather-latency-limited, or mixed because
hardware-counter access is disabled by the host policy.

FFTW is the only CPU-convolution provider wired through the current
`CPUFFTProvider` boundary. Its parity and threading were measured. oneMKL is
installed, but no oneMKL sparse or CPU-convolution provider was benchmarked:
the current `USE_MKL_FFT` switch selects the historical dipole provider, while
`CPUConvolution` is compiled only with `USE_FFTW`; the old sparse MKL path is
also explicitly disabled because it relies on obsolete entry points. A modern
provider would require a new optional provider layer, persistent handle
ownership, and parity/benchmark targets. No broad provider refactor was
introduced under HAM-10.

## Machine and build provenance

| Item | Value |
|---|---|
| Host | `alcazar` |
| OS/kernel | Ubuntu Linux, `6.8.0-137-generic`, x86-64 |
| CPU | 11th Gen Intel Core i9-11900 @ 2.50 GHz |
| Topology | 1 socket, 8 physical cores, 2 SMT threads/core, 16 logical CPUs |
| NUMA | 1 node; CPUs 0–15; approximately 64 GiB RAM |
| Affinity | `taskset -c 0-7`; one logical CPU per physical core |
| OpenMP | `OMP_DYNAMIC=FALSE`, `OMP_PLACES=cores`, `OMP_PROC_BIND=close` |
| Compiler | GNU Fortran 13.3.0 (`/bin/f95`), GCC 12.4.0 C/C++ |
| Build | CMake Release, CPU-only, double precision, OpenMP, FFTW enabled, MKL disabled |
| Fortran flags | `-ffree-line-length-0 -std=legacy -cpp -march=native -O3 -DNDEBUG -Ofast -funroll-loops -finline-functions` |
| FFTW | 3.3.10; `libfftw3_threads` available and used |
| oneMKL | 2026.1.0 installed at `/opt/intel/oneapi/mkl/2026.1`, not selected in measured build |
| VTune/LIKWID | Not installed/available as command-line tools |
| Git state at freeze | Dirty only from pre-existing untracked `build/`, `lib/`, and CPU-HAM prompt-pack files; no tracked source changes |

The measured binary was built in `/tmp/uppasd_cpu_ham10_build` from the source
commit above. The repository `build/` tree was not used for the final CPU
results because it is a CUDA/single-precision build.

## Method

Production measurements used
`benchmarks/harness/cpu_ham05_driver.py` and the real `sd.f95` executable.
Each fit used two independent complete-process samples and fit process wall
time to `T(nstep) = setup + nstep * steady_step`. DIRECT and SPARSE used
5/10/20 steps for the Fe and Nd controls; CONVOLUTION used 50/100/200 steps
where its setup cost made shorter fits noisy. The short-range control used
20/40/80 steps. The production metric is a complete ASD step, not a pair-only
timer. The interaction rate is

`2 * directed_interactions * ensembles / steady_step`

and is reported as million directed interactions per ASD step per second.

The HAM-09 stage benchmark was also run from the isolated CPU-only build with
five raw samples after warm-up. It reports persistent setup, pack, forward FFT,
spectral multiply, inverse FFT, unpack, total apply, and DIRECT parity.

## Production backend results

The table reports median setup and steady-state fit values. `raw steady` lists
the two independent fit slopes in seconds; large disagreement is retained as
dispersion evidence rather than hidden.

| Workload / size | Backend / threads | Setup (s) | Steady (s) | Raw steady (s) | M interactions/s |
|---|---:|---:|---:|---|---:|
| Nd 16³, 4 basis, 1338/atom | DIRECT / 8 | 0.5913 | 0.005437 | 0.005643, 0.005231 | 8064 |
| Nd 16³ | SPARSE / 8 | 0.6309 | 0.020934 | 0.027919, 0.013950 | 2094 |
| Nd 16³ | CONVOLUTION / 8 | 0.7541 | 0.001034 | 0.001077, 0.000990 | 42420 |
| Nd 20³, 4 basis, 1338/atom | DIRECT / 8 | 1.0003 | 0.012698 | 0.014679, 0.010718 | 6744 |
| Nd 20³ | SPARSE / 8 | 1.1114 | 0.025377 | 0.025267, 0.025488 | 3374 |
| Nd 20³ | CONVOLUTION / 8 | 1.3119 | 0.002205 | 0.002081, 0.002330 | 38835 |
| Fe 13³, 2 basis, 96/atom | DIRECT / 8 | 0.1324 | 0.003409 | 0.003664, 0.003154 | 2475 |
| Fe 13³ | SPARSE / 8 | 0.1393 | 0.004246 | 0.003941, 0.004551 | 1987 |
| Fe 20³, 2 basis, 96/atom | DIRECT / 8 | 0.4445 | 0.015207 | 0.015068, 0.015346 | 2020 |
| Fe 20³ | SPARSE / 4 | 0.4588 | 0.026876 | 0.027138, 0.026614 | 1143 |
| Short scalar-J 16³, 6/atom | DIRECT / 8 | 0.0147 | 0.000213 | 0.000220, 0.000205 | 231 |
| Short scalar-J 16³ | SPARSE / 8 | 0.0154 | 0.000235 | 0.000238, 0.000231 | 209 |
| Short scalar-J 16³ | CONVOLUTION / 8 | 0.0296 | 0.000559 | 0.000641, 0.000478 | 88 |
| Short scalar-J 32³, 6/atom | DIRECT / 8 | 0.0602 | 0.001453 | 0.001455, 0.001451 | 271 |
| Short scalar-J 32³ | SPARSE / 8 | 0.0613 | 0.001915 | 0.002127, 0.001704 | 205 |
| Short scalar-J 32³ | CONVOLUTION / 8 | 0.0646 | 0.002227 | 0.002224, 0.002229 | 177 |

The Fe 20³ SPARSE 8-thread fit had raw slopes 0.004108 and 0.021292 s with
fit RMSE 0.851 and 0.002 s respectively. It is treated as an unstable
outlier, not as a sparse crossover; the stable 4-thread comparison remains
slower than DIRECT.

REDUCED-DIRECT was not added to this production matrix. HAM-08 keeps it as an
explicit experimental control rather than a production backend, and the
current HAM-05 production harness exposes only the supported DIRECT, SPARSE,
and CONVOLUTION choices. The current-host matrix therefore does not claim a
REDUCED-DIRECT result or recommendation.

### DIRECT thread sweep

The current production harness completed the 1/2/4/8 sweep at Nd 16³, Fe
20³, and short scalar-J 32³. Median steady times were:

| Workload | 1 thread | 2 threads | 4 threads | 8 threads |
|---|---:|---:|---:|---:|
| Nd 16³ DIRECT | 37.991 ms | 18.749 ms | 10.231 ms | 5.437 ms |
| Fe 20³ DIRECT | 49.126 ms | 28.474 ms | 17.874 ms | 15.207 ms |
| Short 32³ DIRECT | 3.438 ms | 2.536 ms | 1.788 ms | 1.453 ms |

These are timing/scaling observations only. Without counters they do not
classify the hardware bottleneck.

### Setup amortization

For Nd 16³ at eight threads, CONVOLUTION pays approximately 0.163 s more
setup than DIRECT but saves approximately 4.403 ms per ASD step, giving a
measured-fit break-even near 37 steps. At Nd 20³ the corresponding estimate
is approximately 30 steps. Against SPARSE, CONVOLUTION breaks even near 6–9
steps in these two sizes. These values apply only to the measured fixtures;
they do not justify AUTO policy.

## HAM-09 convolution stage and provider measurements

At `OMP_NUM_THREADS=2`, `UPPASD_FFT_THREADS=1`, the five-sample median stage
times from the CPU-only HAM-09 benchmark were:

| Fixture | Pack (ms) | Forward (ms) | Spectral (ms) | Inverse (ms) | Unpack (ms) | Apply (ms) | Parity max |
|---|---:|---:|---:|---:|---:|---:|---:|
| Nd scalar-J, 25³, 4 basis | 0.107 | 0.545 | 0.179 | 0.439 | 0.106 | 1.378 | 1.5e-11 |
| Fe scalar-J, 20³, 2 basis | 0.059 | 0.170 | 0.045 | 0.157 | 0.054 | 0.489 | 5.0e-14 |
| Short scalar-J, 32³ | 0.077 | 0.168 | 0.049 | 0.160 | 0.068 | 0.539 | 5.6e-17 |
| 2D J+D, 32² | 0.005 | 0.011 | 0.005 | 0.012 | 0.005 | 0.041 | 2.8e-17 |
| 3D J+D, 16³, 2 basis | 0.051 | 0.079 | 0.106 | 0.075 | 0.041 | 0.365 | 1.2e-16 |

For the 25³ Nd scalar-J fixture, outer OpenMP thread raw apply samples in
milliseconds were:

| Outer threads | Raw samples | Median |
|---:|---|---:|
| 1 | 1.576, 1.557, 1.530, 1.530, 1.522 | 1.530 |
| 2 | 1.438, 1.424, 1.371, 1.420, 1.384 | 1.420 |
| 4 | 1.237, 1.238, 1.248, 1.240, 1.262 | 1.240 |
| 8 | 1.209, 1.170, 1.181, 1.179, 1.194 | 1.181 |

With outer OpenMP fixed at one, FFTW provider-thread raw samples for the same
Nd fixture were:

| FFTW threads | Raw samples (ms) | Median |
|---:|---|---:|
| 1 | 1.586, 1.551, 1.538, 1.531, 1.533 | 1.538 |
| 2 | 1.570, 1.536, 1.511, 1.502, 1.503 | 1.511 |
| 4 | 1.535, 1.507, 1.491, 1.493, 1.488 | 1.493 |
| 8 | 1.640, 1.607, 1.590, 1.594, 1.628 | 1.607 |

The other scalar-J medians from the provider sweep were Fe 0.512, 0.505,
0.506, and 0.582 ms for 1/2/4/8 FFTW threads, and short-range 0.523, 0.518,
0.558, and 0.593 ms. There is no robust reason to give FFTW ownership more
than one thread when the outer production path already owns OpenMP threads;
the best measured outer/provider combinations remain workload-dependent.

## Hardware-counter campaign

`perf` is installed, but both unprivileged and escalated attempts failed with
the host policy `perf_event_paranoid=4` and the message that CAP_PERFMON,
CAP_SYS_PTRACE, or CAP_SYS_ADMIN is required. VTune and LIKWID were not
available. Consequently the following were not claimed:

- DRAM bandwidth or memory-bandwidth saturation;
- L1/L2/LLC miss rates or cache-bound fraction;
- IPC/CPI, branch, TLB, vectorization, or frequency counters;
- a proof that Nd DIRECT is DRAM-bandwidth-limited versus gather-latency or
  mixed.

The observed DIRECT scaling is compatible with a memory/gather bottleneck,
but that is an inference from timing only and is not a HAM-10 hardware verdict.

## oneMKL provider decision

oneMKL 2026.1 exposes the modern `mkl_sparse_d_create_csr`,
`mkl_sparse_d_mm`, and `mkl_sparse_optimize` interfaces on this host. The
repository does not yet connect them to the production CPU sparse object:

- `source/System/sparse.f90` contains the historical `__INTEL_MKL__` path;
- the current CMake path sets `UPPASD_MKL_LEGACY_SPARSE=0` because its old
  `mkl_dbsrmv`/`mkl_dcsrmm` entry points are not exported;
- the production HAM-03B sparse implementation is a separate persistent
  portable CSR loop in `hamiltonianactions.f90`;
- a probe with `USE_MKL=ON; USE_MKL_FFT=ON; USE_FFTW=OFF` configures only
  `fftdipole_mkl.f90`; it does not add `CPUFFTProvider` or `CPUConvolution`.

This is an explicit STOP for the provider implementation portion of HAM-10.
The portable sparse backend remains the supported correctness reference, and
the measured results do not promote it to a general performance default.
No oneMKL sparse parity or thread-ownership result is fabricated.

The same boundary prevents a clean oneMKL convolution comparison in this
task. A future implementation should add an optional DFTI-backed provider
under the existing provider interface, retain FFTW as a fallback, and then
rerun provider parity, threading, and crossover measurements. oneMKL must not
become a mandatory UppASD dependency.

## Correctness and final decisions

The seven relevant CPU tests pass in the isolated CPU-only build:

```text
cpu-ham-field-energy-contract
cpu-ham-backend-policy
cpu-ham-sparse-backend
cpu-ham-target-order
cpu-ham-reduced-stencil
cpu-ham-convolution
cpu-ham-j-plus-d
```

The HAM-09 scalar-J and J+D convolution fixtures pass DIRECT parity, including
multi-basis/multi-ensemble cases, with the parity maxima shown above. No
unsupported AUTO policy was introduced.

1. **Is Nd DIRECT bandwidth/gather limited?** Not proven; counters were inaccessible. Timing suggests a memory/gather-sensitive or mixed bottleneck, but this remains a hypothesis.
2. **Does optimized sparse library execution add value?** Unanswered; oneMKL SpMM was not integrated. Portable SPARSE did not add value in these measured controls.
3. **Which sparse implementation remains supported?** The portable persistent CSR implementation remains supported as the optional correctness/reference backend; no oneMKL production recommendation is made.
4. **FFTW or oneMKL?** Only FFTW was measured. FFTW remains the supported CPU convolution provider for this revision.
5. **Should both remain?** FFTW is robust enough for the current path; oneMKL remains future optional work, not a recommendation from this campaign.
6. **Where are the crossovers?** CONVOLUTION wins the measured long-range Nd sizes; DIRECT wins the measured short-range sizes and the production Fe control. Fe scalar-J convolution timing is covered by the HAM-09 microbenchmark, but the production B01 Fe input is not convolution-eligible.
7. **Does evidence justify AUTO?** No. Backend selection remains explicit.

## Limitations

- Hardware counters and profiler analyses could not be collected under the host
  security policy.
- The production matrix uses two independent samples per fit cell; the
  HAM-09 stage benchmark uses five samples. This is enough to expose the large
  crossovers, but not a general statistical model.
- The complete all-backend, all-size, all-provider matrix requested by the
  profile is not possible without the unimplemented oneMKL providers and the
  experimental REDUCED-DIRECT control. Raw measurements are retained outside
  the repository under `/tmp` and the report includes the representative raw
  samples used for decisions.
- Production executable timing exports complete process/step timing rather
  than a dedicated pair-only DIRECT timer.
