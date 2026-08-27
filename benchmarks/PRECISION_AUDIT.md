# Precision audit (WP-06 section A)

This document records what `UPPASD_PRECISION` actually changes, established
by reading the source rather than assumed from the CMake option name
(blueprint section 11: "Do not infer numerical precision solely from CMake
option names"). It backs `benchmarks/harness/precision_audit.py`, which
encodes the same table and functions as code; if the two ever disagree, the
source citations here (and in that module's docstring) are what to
re-verify against real `source/` files, not this document's prose.

## Headline finding

**`UPPASD_PRECISION=SINGLE` does not make UppASD's CPU arithmetic single
precision.** The CPU (Fortran) numerical path is double precision in every
supported build. Only the GPU (CUDA) numerical path responds to
`UPPASD_PRECISION`.

```
requested_precision=SINGLE  ->  effective_cpu_precision=DOUBLE, effective_gpu_precision=SINGLE
requested_precision=DOUBLE  ->  effective_cpu_precision=DOUBLE, effective_gpu_precision=DOUBLE
requested_precision=MIXED   ->  rejected at CMake configure time; no build exists
```

## Why

`source/Parameters/parameters.f90:8` defines the CPU working kind
unconditionally:

```fortran
integer, parameter :: dblprec = selected_real_kind(15, 307)
```

No `#ifdef`/CPP guard surrounds it, and no `.f90`/`.F90` file anywhere in
`source/` tests the `SINGLE_PREC` macro (a repo-wide search found zero
matches outside `source/gpu_files/` C++ headers and legacy `.make`
profiles). Every CPU numerical routine — the moment state, the Hamiltonian,
the Depondt integrator — declares its reals `real(dblprec)`.

`CMakeLists.txt:110-111` declares `UPPASD_PRECISION` (`DOUBLE`/`SINGLE`/
`MIXED`); `CMakeLists.txt:158-166` validates it and fatal-errors on `MIXED`
("not implemented yet; use DOUBLE or SINGLE") — there is no way to configure,
let alone build, a MIXED binary. `CMakeLists.txt:286-288` is the only place
`UPPASD_PRECISION=SINGLE` becomes a compiler define:

```cmake
if(UPPASD_PRECISION STREQUAL "SINGLE")
  target_compile_definitions(${UppASD_LIB} PRIVATE SINGLE_PREC)
endif()
```

`SINGLE_PREC` is applied to the one shared `asdlib` target (Fortran and
CUDA/HIP sources together), but as a C-preprocessor define it is inert for
`.f90` translation units — nothing there tests it. It is exclusively the GPU
C++ headers under `source/gpu_files/` that key their `real` typedef off it
(`real_type.h:50-58`), so `SINGLE_PREC` only ever changes GPU-side storage
and arithmetic.

The Fortran↔C++ ABI boundary is deliberately fixed at double regardless of
build precision — `fortranData.hpp:8-10`:

> The Fortran ABI always supplies `real(dblprec)` storage. Keep this type
> independent from the precision selected for device-side `real`.

The staging conversion between the two lives in `tensor.hpp:124-148`
(`Tensor::set()`/`copy_to()`): "the Fortran interface always owns double
precision buffers while a SINGLE_PREC build stores real values as float."

## Component table

| # | Component | DOUBLE build | SINGLE build | Citation |
|---|---|---|---|---|
| 1 | CPU magnetic moments / spin state | DOUBLE | DOUBLE (unchanged) | `source/Parameters/parameters.f90:8`; `source/System/momentdata.f90:10-16` |
| 2 | CPU Hamiltonian evaluation | DOUBLE | DOUBLE (unchanged) | `source/Hamiltonian/*.f90` (`real(dblprec)` throughout) |
| 3 | CPU Depondt integrator | DOUBLE | DOUBLE (unchanged) | `source/Evolution/depondt.f90:62-85` |
| 4 | GPU spin state (device emom/emomM/mmom) | DOUBLE | SINGLE | `gpu_files/real_type.h:50-58`; `gpuStructures.hpp:213-223` |
| 5 | GPU Hamiltonian kernels | DOUBLE | SINGLE | `gpuHamiltonianCalculations.hpp:26-65` |
| 6 | GPU Depondt kernels | DOUBLE | SINGLE | `gpuDepondtIntegrator.cpp:21-61` |
| 7 | Thermal field (curand) | DOUBLE (`curandGenerateNormalDouble`) | SINGLE (`curandGenerateNormal`) | `gpu_wrappers.h:124-131`; CPU-side RNG stays `real(dblprec)` unconditionally |
| 8 | Energy / measurement reductions | DOUBLE | SINGLE partial sums, widened to DOUBLE before returning to Fortran | `measurement/kernels.hpp:36,49,145`; `measurement/gpuMeasurement.cpp:835,972-973`; CPU side stays `real(dblprec)`: `source/Hamiltonian/energy.f90:49,108` |

Components 4-6 all key off the same `SINGLE_PREC`-gated `real` typedef, so
`effective_gpu_precision` is one value per `(requested_precision,
gpu_backend)` pair, not tracked per component
(`precision_audit.effective_gpu_precision`). Component 8's *device-side*
working precision still follows `SINGLE_PREC`; only the final host-bound
value is widened, so it is not modelled as a separate state either.

**HIP** shares the same `gpu_files/` headers textually, but per the
blueprint's "no HIP performance claim without real hardware" rule this audit
does not claim to have *verified* HIP's effective precision against real
execution — `precision_audit.effective_gpu_precision` reports `UNKNOWN` for
`gpu_backend="HIP"`, matching `gpu_provenance.py`'s already-unaudited HIP
device metadata.

## Consequence: no analogous CPU SINGLE/DOUBLE ratio

Since `effective_cpu_precision` is DOUBLE for every buildable configuration,
a "CPU SINGLE" build is numerically identical to CPU DOUBLE — there is no
second CPU precision mode to ratio against the first. Section F's `R_GPU_32_64`
is therefore GPU-only; `precision_metrics.compute_r_cpu_32_64` exists purely
to make that refusal explicit (see its docstring) rather than silently
absent.

## Comparison classification

`precision_audit.classify_comparison_precision_class` turns the above into
the schema's `comparison_precision_class` (VALIDITY.md):

* `MIXED` (`unsupported`) -> always `null`.
* Not `COMPLETED` -> always `null` (no execution evidence to back a claim).
* `unaudited` support state (HIP; a non-authoritative developer context) ->
  `UNAUDITED`.
* `supported` + `COMPLETED`:
  * CPU record -> `PRODUCTION_CONFIGURATION` (real production evidence, but
    "matched" is a claim about a pairing, not a standalone CPU record).
  * GPU record, `effective_gpu_precision == effective_cpu_precision` (only
    possible for a GPU DOUBLE build, since CPU is always DOUBLE) ->
    `PRECISION_MATCHED`.
  * GPU record otherwise (a GPU SINGLE build against the always-DOUBLE CPU
    host) -> `PRODUCTION_CONFIGURATION`, never `PRECISION_MATCHED`.

## Host OpenMP configuration for a GPU campaign (section D)

Real host-side, OpenMP-parallel work executes in a GPU production run, not
only the per-step GPU dispatch:

* RNG initialization always runs on the CPU:
  `source/uppasd.f90:823-826` (`!$omp parallel copyin(tseed)`),
  unconditional of `do_gpu`.
* Hamiltonian/neighbour-list and dipole-field setup (once, before GPU
  handoff) is OpenMP-parallel: `source/Hamiltonian/dipolecommon.f90:134,198`,
  `source/Hamiltonian/fftdipole_fftw.f90` (10+ `!$omp parallel do` sites).
- Spin-transfer-torque field construction, used by the GPU driver routines
  (`sd_iphaseGPU`/`sd_mphaseGPU`, `source/sd_driver.f90:951-1105`, which
  `use SpinTorques, only: btorque, stt`), is OpenMP-parallel:
  `source/Fields/spintorques.f90:86,136,205,239`.
* `source/uppasd.f90:1735-1740` prints an OpenMP thread-count banner
  unconditionally at program end, confirming OMP threading is live in every
  build/run regardless of `do_gpu`.

This audit did **not** establish whether the per-step GPU loop itself calls
back into host-OMP work every iteration for cases without STT (that would
need runtime/profiling evidence, not static reading) — only that setup-phase
and STT-conditional host OMP work is real. `GPU_HOST_OMP_THREADS_DEFAULT = 1`
(`gpu_campaign.py`) is therefore used as the blueprint's preferred initial
choice, not as a proven-safe conclusion; per section D, "Do not optimize
host OMP+GPU interaction in this task" — a case whose variant uses STT or
dipole/FFT setup is exactly where that choice most needs re-examining in a
later work package.

## Sanity methodology (section G)

`gpu_sanity.sanity_check_cpu_vs_gpu` compares a CPU and a GPU T=0 sample's
final `restart.<simid>.out` moment state for a *gross* physics failure
(NaN/Inf, a diverged magnitude, a flipped direction), reusing
`omp_sanity.compare_restart_moments` — the same mechanism WP-05 already
established for thread-count sanity, which does not care which run
configuration dimension differs between the two samples being compared. The
default magnitude tolerance is precision-aware: `1e-3` relative for a GPU
DOUBLE sample (matching the WP-05 thread-count bound), `1e-2` for a GPU
SINGLE sample (documented, not derived — fp32 accumulation over many
timesteps legitimately drifts further from an fp64 CPU reference). This is
not a substitute for UppASD's scientific validation suite.
