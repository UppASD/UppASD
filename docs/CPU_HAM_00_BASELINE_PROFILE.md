# CPU-HAM-00 Production CPU Hamiltonian Baseline Profile

Date: 2026-09-02
Repository: `UppASD_gpu_hip_cu`
Branch: `gpu_hip_cu_ab_cg`

## Scope and outcome

This report follows `docs/cpu/UppASD_CPU_Hamiltonian_Luna_Prompt_Pack/CPU_HAMILTONIAN_OPTIMIZATION_BLUEPRINT.md` and the CPU-HAM-00 procedure in `docs/cpu/UppASD_CPU_Hamiltonian_Luna_Prompt_Pack/CPU-HAM-00_PROFILE_BASELINE.md`.

The measured path is the ordinary production CPU LLG path through
`source/Hamiltonian/hamiltonianactions.f90`. No Hamiltonian source optimization was
made or committed. A separate, user-requested configuration-only commit disabled
active `do_prnstruct 1` settings in tests and benchmarks: `efe1c6d`.

The baseline is complete for this host and toolchain, with the limitations recorded
below. Hardware performance counters and operating-system sampling were unavailable,
so memory-bandwidth saturation, IPC, cache-miss rates, and stall cycles are not
claimed as measured facts. The long-range Nd result is therefore classified as a
working hypothesis supported by source inspection, compiler diagnostics, and scaling;
it must be rechecked with counters on a Linux/VTune/LIKWID-capable host before
making implementation decisions.

## Build and host

| Item | Value |
|---|---|
| Host | Apple M1, Darwin 25.5, arm64 |
| Hardware reported | 8 physical/logical processors, 16 GiB memory (`hostinfo`) |
| Compiler | `/opt/homebrew/bin/gfortran` |
| Binary | `build/bin/sd.gfortran` |
| Build mode | CMake Release, CPU-only (`UPPASD_CMAKE_USE_CUDA=OFF`, HIP off) |
| Relevant flags | `-O3 -Ofast -funroll-loops -finline-functions -mtune=native -fopenmp` |
| OpenMP sweep | `OMP_NUM_THREADS=1,2,4,8`, `OMP_DYNAMIC=FALSE` |
| Standard placement request | `OMP_PLACES=cores`, `OMP_PROC_BIND=close` |
| Timing runs | 16,384 atoms, structure printing disabled, two repetitions per thread count |

The benchmark binary ran successfully for every recorded run. The source files
`source/Hamiltonian/hamiltonianactions.f90` and `source/sd_driver.f90` were unchanged
by this profiling task.

## Production call path

The benchmark cases use CPU mode and the normal measurement phase:

```text
uppasd/main
  -> run_measurement_phase
    -> sd_mphase
      -> per-step measurement and cadence-gated calc_energy
      -> effective_field (predictor)
        -> effective_field_full
          -> heisenberg_field and enabled term kernels
          -> final field assembly and scalar energy reduction
      -> evolve_first
      -> effective_field (corrector)
        -> effective_field_full
      -> evolve_second
      -> moment_update
```

The dispatch is the `effective_field` generic interface in
`source/Hamiltonian/hamiltonianactions.f90:24`, with the ordinary full-field
implementation beginning at `:59`. In `source/sd_driver.f90`, the normal non-sparse
branch calls `effective_field` for the predictor at `:539-541` and for the corrector
at `:580-599`, around the evolution calls at `:553-620`. Consequently, the Heun and
Depondt SDE paths represented here execute two Hamiltonian field evaluations per
production timestep.

The `do_sparse` branch was not active in these cases. The obsolete
`ApplyHamiltonian/heisge()` path was excluded from the audit, as required by the
blueprint; it is not the production owner-computes path.

## Workloads

| Case | System and size | Directed interactions | Ensemble | Active Hamiltonian terms | Integrator / phase notes |
|---|---:|---:|---:|---:|---|
| B01 | bcc Fe, 20x20x20, 16,000 atoms | 96 per atom; 1,536,000 total | 10 | Exchange; uniaxial/cubic onsite anisotropy | Depondt SDE, `SDEalgh=5`; no initial phase |
| B02 | 2-D skyrmion, 128x128, 16,384 atoms | 4 exchange + 4 DMI entries per atom | 1 | Short-range exchange + DMI | Heun, `SDEalgh=1`; `ip_mode=Y`, fixed initial q-search |
| B04 | dhcp Nd, 16x16x16, 16,384 atoms | median 1,338, maximum 1,340 per atom | 1 | Very long-range exchange | Default Heun; no initial phase |

The interaction counts distinguish the exchange-list metric used by the benchmark
harness from the total pair-list work in B02: the B02 Hamiltonian field evaluates
both its exchange and DMI neighbor loops.

## Wall-clock baseline

The following are medians of two complete-process runs at `Nstep=20`. Structure
printing and other cadence outputs were disabled for these timing runs. Speedup is
relative to the one-thread median; efficiency is speedup divided by thread count.

| Case | 1 thread (s) | 2 threads (s) | 4 threads (s) | 8 threads (s) | 8-thread speedup | 8-thread efficiency |
|---|---:|---:|---:|---:|---:|---:|
| B01 bcc Fe | 1.823914 | 1.324315 | 1.076214 | 1.013097 | 1.800x | 0.225 |
| B02 skyrmion | 0.624710 | 0.555494 | 0.542668 | 0.543595 | 1.149x | 0.144 |
| B04 dhcp Nd | 1.891668 | 1.329318 | 0.925703 | 0.802901 | 2.355x | 0.294 |

Raw measurements, in execution order, are preserved here rather than only reporting
the medians:

```text
case             threads  rep  wall_s
B01_bccFe        1        1    1.854964
B01_bccFe        1        2    1.792864
B01_bccFe        2        1    1.326382
B01_bccFe        2        2    1.322248
B01_bccFe        4        1    1.089371
B01_bccFe        4        2    1.063056
B01_bccFe        8        1    1.015462
B01_bccFe        8        2    1.010731
B02_skyrmion2D   1        1    0.663062
B02_skyrmion2D   1        2    0.586358
B02_skyrmion2D   2        1    0.555610
B02_skyrmion2D   2        2    0.555378
B02_skyrmion2D   4        1    0.544100
B02_skyrmion2D   4        2    0.541236
B02_skyrmion2D   8        1    0.540100
B02_skyrmion2D   8        2    0.547089
B04_dhcpNd       1        1    1.887384
B04_dhcpNd       1        2    1.895952
B04_dhcpNd       2        1    1.246534
B04_dhcpNd       2        2    1.412102
B04_dhcpNd       4        1    0.860415
B04_dhcpNd       4        2    0.990990
B04_dhcpNd       8        1    0.841776
B04_dhcpNd       8        2    0.764025
```

Built-in phase timing was also recorded with `Nstep=60`. It is rounded by the
executable, so it is useful for attribution but not a replacement for the wall-clock
measurements:

| Case | Measurement phase, 1 thread | Measurement phase, 8 threads | Approx. phase speedup | Interpretation |
|---|---:|---:|---:|---|
| B01 | 4.47 s | 1.80 s | 2.48x | Pair work plus onsite terms; moderate scaling |
| B02 | 0.15 s | 0.13 s | 1.15x | Steady phase is short; fixed initial q-search dominates whole process |
| B04 | 3.68 s | 0.94 s | 3.91x | Long neighbor list dominates steady-state work |

B02's initial phase was approximately 0.46-0.50 s at all thread counts, explaining
why its complete-process wall time is nearly flat even though the production
measurement phase changes modestly.

## Hamiltonian term attribution

The source-level term ownership in `effective_field_full` is:

1. `heisenberg_field` (`hamiltonianactions.f90:186-206`) gathers each neighbor
   moment through `emomM(:,ham%nlist(j,i,k))` and multiplies by the exchange
   coupling. This is the dominant pair kernel for B01 and B04.
2. `dzyaloshinskii_moriya_field` (`:286-314`) performs the analogous indirect
   neighbor gather through `ham%dmlist`; it is active in B02.
3. Onsite anisotropy and other optional terms execute in the same target-atom
   OpenMP loop (`effective_field_full:114-180`). B01 enables uniaxial/cubic
   anisotropy; B02 and B04 do not enable onsite anisotropy.
4. The finalization and energy expression are inside the same loop. The field
   writes `beff1` and `beff2`, assembles `beff`, computes the local energy
   contribution, and reduces it (`:170-181`). This scalar energy reduction is
   unconditional for every `effective_field` call, even when the caller does not
   use the returned value for integration.

The ordinary measurement phase therefore contains both Hamiltonian work and
non-Hamiltonian work: predictor/corrector evolution, thermal/RNG handling for the
selected integrator, moment updates, and measurement bookkeeping. The available
phase timer does not split those categories. `source/Tools/profiling.f90:149-264`
contains the category-timing routine, but it returns at `:182` before doing any
measurements, so it cannot provide valid term-level timing in this build.

## Scaling interpretation and bottleneck classification

| Case | Classification | Evidence | Confidence |
|---|---|---|---|
| B01 bcc Fe | `MIXED/UNRESOLVED` | 96 indirect exchange neighbors/atom, onsite anisotropy, and only 1.80x whole-process speedup at 8 threads; compiler reports no vectorization of the call-heavy neighbor loop | Medium-low without counters |
| B02 skyrmion | `MIXED/UNRESOLVED` | Short lists and both exchange/DMI gathers; fixed `ip_mode=Y` initial phase dominates complete-process timing; steady phase reaches only about 1.15x | Medium-low |
| B04 dhcp Nd | `CACHE_LATENCY/GATHER` working hypothesis | Approximately 1,338 indirect exchange gathers/atom and reduced parallel efficiency as the thread count rises; compiler rejects vectorization of the neighbor loop | Medium-low; counters unavailable |

The B04 classification is deliberately labelled a hypothesis. The source access
pattern strongly suggests irregular gather/cache pressure, and the scaling trend is
consistent with it, but no LLC-miss, bandwidth, or stall-cycle measurement was
available. It is not valid to state that B04 has reached a measured memory-bandwidth
roofline on this host.

## Hardware-counter and sampling status

| Requested evidence | Status |
|---|---|
| Memory bandwidth | Unavailable: `perf` and `likwid-perfctr` are not installed; Instruments/privileged performance tools were not available |
| IPC / retired instructions | Unavailable |
| LLC/cache misses | Unavailable |
| Memory-stall cycles | Unavailable |
| CPU frequency / residency | Unavailable: `powermetrics` requires superuser access and relevant `sysctl` access was blocked |
| Statistical call-stack sampling | Unavailable: macOS `sample` could not examine the process; DTrace listing was blocked by system policy |
| Compiler vectorization diagnostics | Available; compilation succeeded with `-fopt-info-vec-all` |
| OpenMP placement | Not verifiable: libgomp reported `Affinity not supported on this configuration` for `OMP_PLACES=cores`/`OMP_PROC_BIND=close` and `spread` |
| NUMA topology | Unavailable/not exposed by the permitted host interfaces |

The vectorization report is consistent with the source audit: the full-array clear
at `effective_field_full:99` was vectorized using 16-byte vectors, while the main
OpenMP target loop and the indirect neighbor loops in `heisenberg_field` and DMI
were not vectorized. The compiler specifically reports calls that clobber memory,
unprofitable loop vectorization, and unsuitable gather accesses. Some small
component/SLP blocks were vectorized; this does not make the indirect neighbor
reduction vectorized.

## OpenMP placement experiment

Three B04 runs at 8 threads and `Nstep=60` were made for each setting:

| Setting | Raw wall times (s) | Median (s) | Relative to `close` |
|---|---|---:|---:|
| `OMP_PLACES=cores`, `OMP_PROC_BIND=close` | 1.500105, 1.336854, 1.383836 | 1.383836 | baseline |
| `OMP_PLACES=cores`, `OMP_PROC_BIND=spread` | 1.408007, 1.325132, 1.399331 | 1.399331 | +1.1% |
| Unbound (`OMP_PROC_BIND=FALSE`) | 1.321539, 1.446579, 1.787475 | 1.446579 | +4.5% |

The run-to-run variation and libgomp affinity warning prevent a strong placement
conclusion. `close` is retained as the reproducible requested setting for the
baseline commands, but actual thread affinity was not established by this runtime.

## Energy-cadence overhead

The field kernel's local scalar energy reduction is always executed. Separately,
`sd_mphase` calls the term-resolved `calc_energy` pass only when the configured energy
cadence fires (`source/sd_driver.f90:431-455`). A B01 diagnostic with `Nstep=10`
compared cadence-off (`ene_step=11`) with `calc_energy` every step (`ene_step=1`):

| Threads | Cadence off median (s) | Every-step median (s) | Difference |
|---:|---:|---:|---:|
| 1 | 1.164342 | 1.408442 | +20.96% |
| 8 | 0.839453 | 0.844710 | +0.63% |

These are complete-process timings over a short run, so the 8-thread difference is
within the observed noise. They do establish that the separate energy pass can be a
material low-thread/short-run cost, while it is not the explanation for the
unconditional reduction inside `effective_field_full`.

## Full-array pass audit

| Operation | Location | Classification for this baseline |
|---|---|---|
| `beff=0.0` | `hamiltonianactions.f90:99` | Semantically required in the current implementation because optional dipole/precomputed fields are accumulated before per-term assembly; vectorized by the compiler |
| Local `beff_s`, `beff_q`, `beff_m` initialization | `hamiltonianactions.f90:118-120` | Required working-state initialization for each target atom |
| `energy=0.0` and local energy accumulation | `hamiltonianactions.f90:97`, `:175-181` | Required for the returned Hamiltonian energy; currently unconditional |
| `thermal_field=0.0` / evolution working buffers | `source/Evolution/evolution.f90` and `source/Evolution/depondt.f90` | Required by the current predictor/corrector buffer ownership; candidate for a later, separately verified optimization |
| `bdup(:,:,:)=0.0` | `source/Evolution/depondt.f90` | Working-buffer reset required by current Depondt accumulation semantics; candidate for later redesign |
| Time-dependent external-field clear | `source/Fields/calculatefields.f90:139-140` | Required to keep inactive/no-field inputs defined in the current interface; could be specialized only after semantic tests |
| `beff1`/`beff2` per-atom writes | `hamiltonianactions.f90:170-173` | No separate full-array clear is needed; each target element is overwritten before use |
| `emom`/`emomM`/`mmom` copies | `source/Evolution/updatemoments.f90` | Required state propagation for the selected integrator; not a redundant Hamiltonian clear |

No pass in this audit was removed or changed. The table records candidates and
semantic constraints for a future optimization pass.

## CPU-HAM-00 checklist

- [x] Read the CPU Hamiltonian optimization blueprint.
- [x] Read and follow `CPU-HAM-00_PROFILE_BASELINE.md`.
- [x] Identify the production CPU owner-computes Hamiltonian implementation.
- [x] Verify the ordinary production call path and predictor/corrector call count.
- [x] Profile bcc Fe medium-range exchange (B01).
- [x] Profile short-range exchange + DMI (B02).
- [x] Profile very long-range Nd exchange at the poor-scaling benchmark size (B04, 16,384 atoms).
- [x] Disable active `do_prnstruct 1` settings in tests and benchmarks before timing; committed separately as `efe1c6d`.
- [x] Record workload sizes, neighbor counts, ensemble sizes, integrators, and cadence assumptions.
- [x] Record raw wall-clock measurements and thread scaling.
- [x] Attempt hardware-counter and sampling attribution.
- [x] Record unavailable metrics explicitly rather than inferring them.
- [x] Check compiler vectorization diagnostics for the production Hamiltonian source.
- [x] Audit full-array clears and copies, classifying semantic requirements and candidates.
- [x] Compare energy-cadence overhead.
- [x] Check OpenMP placement behavior and record runtime limitations.
- [x] Classify bottlenecks without claiming unsupported counter evidence.
- [x] Make no CPU Hamiltonian optimization changes in this baseline task.

## Verification and handoff

All current-build benchmark runs recorded above returned exit code 0. The production
Hamiltonian source also compiled successfully with the normal optimization flags
plus `-fopt-info-vec-all`; the resulting report was inspected for the full-field,
exchange, and DMI loops. Repository whitespace validation was run with
`git diff --check` before commit.

The next optimization task should begin with a counter-capable Linux or profiling
host, then validate one hypothesis at a time against this report. The highest-value
candidate is the B04 indirect long-range exchange gather/reduction, followed by
carefully measured B01 pair-list and mixed onsite work. B02 should be measured with
its fixed initial phase separated from steady-state Hamiltonian timing.
