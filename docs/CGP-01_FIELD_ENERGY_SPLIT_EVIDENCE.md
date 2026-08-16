# CGP-01 — Field/energy split evidence

Task: `docs/CGP_work.md` CGP-01 ("Split field evaluation from energy
evaluation"). Branch `gpu_hip_cu_ab_cg`, commit base `1374c9bc` (CGP-00B).

## Environment

- GPU: 2x NVIDIA RTX A4000, driver visible via `nvidia-smi`; both idle
  (0% utilization, <110 MiB used) immediately before every timed run below
  — see the shared-host-contention note in project memory. Comparison runs
  still followed the memory's interleaving guidance (ABAB before/after
  alternation, no single-arm run trusted alone).
- CUDA 13.3.73 (nvcc), gfortran 13.3.0, GCC 12.4.0, CMake 3.28.3.
- "After" build: fresh out-of-tree `build_cgp01_cuda_fp32`
  (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`, `BUILD_TESTING=ON`,
  `CMAKE_BUILD_TYPE=Release`) at the working tree with this task's changes.
- "Before" build: identical CMake configuration, built from a `git worktree`
  pinned to commit `1374c9bc` (this task's parent, CGP-00B), so the
  comparison isolates exactly this task's diff.

## 1. Inventory (Part A)

Every call site of `GpuAdaptiveRuntime::evaluateHybrid` was traced (grep
across `source/gpu_files/` and `tests/coarse_graining/`) and classified:

| Call site | Classification |
|---|---|
| `gpuSimulation.cpp` "initial" field, every step | field-only |
| `gpuSimulation.cpp` "predictor" field, every step | field-only (energy only needed when `adaptiveDiagnostics>0`; see Part D) |
| `gpuSimulation.cpp` "corrector"/"reconstructed" field, `adaptiveDiagnostics>=3` only | field + diagnostic energy |
| `gpuSimulation.cpp::release()` teardown `last_energy_j` print | reads cached `lastEnergy()`, does not itself call `evaluateHybrid` |
| `test_gpu_adaptive_runtime.cpp` kernel-parity/continuum-oracle tests | field + user-visible measurement energy |
| `benchmark_gpu_adaptive_runtime.cpp` `adaptiveCoarseStep()` (the headline timing loop) | field-only |
| `benchmark_gpu_adaptive_runtime.cpp` `runHamiltonianOracle`/`reportEnergyStage` | field + measurement energy |

No GPU call site is transition-acceptance energy: CGP-00B established the
GPU adaptive path has no transition-energy evaluator at all (selector/
dwell/polarization-driven only), so that category does not exist on this
backend. The CPU (Fortran) transition-acceptance path
(`evaluate_static_hybrid_operator`) is a separate implementation untouched
by this task.

`lastEnergy_` was unconditionally overwritten by every `evaluateHybrid()`
call, including the field-only production calls whose return value was
immediately discarded via `(void)`.

## 2. Design implemented

`GpuAdaptiveRuntime::evaluateHybridImpl<bool MeasureEnergy>`
(`gpuAdaptiveRuntime.hpp`/`.cpp`) is now the single field implementation
behind two public entry points:

* `evaluateHybrid(...)` — unchanged signature/behaviour, instantiates
  `evaluateHybridImpl<true>`.
* `evaluateField(...)` — new, `void`-returning, instantiates
  `evaluateHybridImpl<false>`.

The eight device kernels that previously fused energy contribution
arithmetic and a `reduceAdaptiveEnergyBlock` block reduction onto field work
(`clearAdaptiveAtomistic`, `evaluateAdaptiveAtomisticBonds`,
`evaluateAdaptiveAtomisticOnsite`, `clearAdaptiveCoarse`,
`evaluateAdaptiveCoarseTensor`, `finalizeAdaptiveCoarseLocal`,
`addAdaptiveDipole`, `addAdaptiveBasisResolvedDipole`) are now templated on
the same `MeasureEnergy` parameter, gating every energy accumulation and
reduction call behind `if constexpr (MeasureEnergy)`. This matches the
`Heisge<..., Measure>` compile-time-specialization precedent already used by
`GpuHamiltonianCalculations` — no runtime `if(measureEnergy)` branch exists
inside any thread; the `MeasureEnergy=false` instantiation simply never
contains the energy code, and callers select the instantiation with an
ordinary `if` at the orchestration level (`evaluateHybridImpl`'s own launch
sites), exactly like `Heisge`'s call sites do.

Field-only launches:

* never populate `kernels.energyPartials`/`energyPartialBlocks` (left at
  the `AdaptiveKernelDevice` defaults, `nullptr`/`0`);
* never zero `kernels.energyTerms` (the `clearAdaptive{Atomistic,Coarse}`
  energy-term zeroing is behind `if constexpr (MeasureEnergy)`);
* never launch `reduceAdaptiveEnergyPartials` or `finalizeAdaptiveEnergy`
  (both calls are behind `if constexpr (MeasureEnergy)` in
  `evaluateHybridImpl`);
* never issue the `double terms[8]` D2H copy;
* never touch `lastEnergy_`.

A new counter, `energyEvaluationCount()`, increments exactly once inside
the `MeasureEnergy=true` tail of `evaluateHybridImpl` (immediately before
`lastEnergy_` is set and the D2H-populated result is returned) and is the
CGP-01 negative control: it must stay unchanged across any number of
`evaluateField()` calls.

### Production dispatch (`gpuSimulation.cpp`)

`evaluateAdaptiveFields` now takes an explicit `measureEnergy` argument
instead of always calling `evaluateHybrid`. Each of the four call sites in
`advanceAdaptiveStep` decides independently, preserving the pre-existing
observable diagnostics contract bit-for-bit rather than reinterpreting it:

* **initial** field: `measureEnergy = (adaptiveDiagnostics>=3 &&
  !adaptiveMaskEnabled)` — identical to `printStaticStageTrace`'s own
  internal gate, since this call's energy is consumed only by that stage
  trace.
* **predictor** field: `measureEnergy = (adaptiveDiagnostics>0)` — wider
  than a literal reading of "diagnostics>=3 requests energy" would suggest.
  This is deliberate: at `diagnostics` 1–2, the predictor call is normally
  the *last* `evaluateHybrid`-equivalent call of a step (corrector/
  reconstructed only fire at `>=3`), so its `lastEnergy_` is exactly what
  `release()`'s teardown `last_energy_j` print reports — a value several
  tracked e2e fixtures assert against numerically (see Part D below), not
  merely for presence. Preserving that existing, test-verified contract at
  `diagnostics` 1–2 took priority over a literal "only level 3" reading;
  this exception is intentional and documented here rather than left as an
  undocumented gap, matching the CGP-00B precedent for handling this kind
  of pre-existing behavioural dependency.
* **corrector**/**reconstructed** fields: unchanged gate
  (`adaptiveDiagnostics>=3 && !adaptiveMaskEnabled`); when they run they
  always request energy, exactly as before this task (every call used to
  request it unconditionally).

Net effect: at `adaptiveDiagnostics=0` (the production/benchmark default),
every field evaluation in a step is field-only — zero energy work per step.
At `diagnostics` 1–2, only the predictor call still requests energy (one
call/step, unchanged from before). At `diagnostics>=3`, all four calls
request energy, unchanged from before.

### Benchmark (`benchmark_gpu_adaptive_runtime.cpp`)

`adaptiveCoarseStep()` — the exact production-mimicking predictor/corrector
pair the headline PERF-CG-SWEEP numbers are timed from — now calls
`evaluateField()` instead of discarding `evaluateHybrid()`'s return value.
The pre-timing all-fine field-parity check was switched to `evaluateField()`
for the same reason: a field-only regression must not be masked by a
passing `evaluateHybrid()` field.

## 3. Field-parity and negative-control evidence (Part E, negative control)

New unit test `testFieldOnlyEvaluation`
(`tests/coarse_graining/test_gpu_adaptive_runtime.cpp`), added to the
`gpu_adaptive_runtime_tests` target:

* Calls `evaluateField()` three times on a non-uniform fixture and asserts
  `energyEvaluationCount() == 0` and `lastEnergy().totalJ == 0.0` after
  each call.
* Calls `evaluateHybrid()` once, asserts `energyEvaluationCount() == 1` and
  a non-trivial energy (`totalJ != 0`).
* Downloads the atomistic and coarse fields written by the field-only calls
  and by the energy call on **identical inputs** and asserts they are
  exactly equal component-by-component (`==`, not a tolerance) — the field
  arithmetic is untouched, so this must hold exactly, and does.
* Calls `evaluateField()` once more after the energy call and asserts the
  counter and `lastEnergy()` are both unchanged, proving a field-only call
  cannot silently republish a stale/cached energy.

`ctest` run (`coarse-graining-gpu-adaptive-runtime`): **PASS**. Full binary
run:

```
CG-09/CG-10 GPU adaptive runtime tests passed
```

## 4. Correctness evidence (`cg13-cuda` label)

`ctest -L cg13-cuda -j2` on `build_cgp01_cuda_fp32`: **32/33 passed**.

The one failure, `adaptive-cg-production-e2e`
(`run_production_e2e.py`, `gpu_fft_static_mixed` case,
`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), was
independently reproduced on the unmodified `1374c9bc` "before" build (same
fixture, same binary invocation): `coarse_dipole=0.0000000000000000e+00`
on **both** builds, byte-identical. This is a pre-existing, unrelated
defect in the open FFT dipole coupling for that specific fixture (already
flagged as pre-existing/unrelated in `docs/CGP-00B_ENERGY_GATE_EVIDENCE.md`
section 4.2) — not a CGP-01 regression. Diffing the full stdout of the
`gpu_fft_static_mixed` case between the before/after builds shows every
`last_energy_j`/`last_field_checksums_t`/stage-trace value agreeing to
~7-8 significant figures (this fixture runs at `cg_diagnostics 3`, so the
before/after code paths request energy identically on every stage; the
residual difference is ordinary FP32 atomic-add scheduling noise, not a
systematic change) and `coarse_dipole=0` identically on both.

`adaptive-cg-production-e2e` additionally exercises `cg_diagnostics 2`
fixtures (`gpu_adaptive_mixed`, `static_all_coarse`, `adaptive_mixed`,
`moving_all_fine`) whose `last_energy_j`/`last_total_energy_j` values are
checked against independent CPU/oracle references elsewhere in the suite
(`run_moving_*`, `run_energy_jump_threshold_sweep.py`); those all passed,
confirming the `diagnostics>0` predictor-energy-preservation design
decision above is correct in practice, not just in reasoning.

## 5. Performance evidence (fp32, primary target)

`benchmark_gpu_adaptive_runtime`'s `adaptiveCoarseStep()` PERF-CG-SWEEP loop
(`--warmup 20 --iterations 200 --repetitions 5`, `setPhaseTimingEnabled(false)`
for the headline number), ABAB before/after alternation, 2 samples per
arm/fraction (median of 2 reported); all-coarse, ~6.25% fine, ~25% fine, and
all-fine fractions are the tracked required points (the sweep also covers
50/75%). `atomistic-parity`/`adaptive-asd-parity` all `PASS` on every run
(no correctness regression at any fraction).

### 65,536 atoms, 4,096 blocks (16 atoms/block)

| fraction | active atoms | before (us) | after (us) | speedup |
|---:|---:|---:|---:|---:|
| 1.000 (all-fine) | 65536 | 204.4 | 106.7 | 1.92x |
| 0.750 | 49184 | 234.4 | 136.6 | 1.72x |
| 0.500 | 32800 | 211.9 | 129.7 | 1.63x |
| 0.250 (~25%) | 16416 | 184.1 | 118.7 | 1.55x |
| 0.125 | 8224 | 167.0 | 114.4 | 1.46x |
| 0.0625 (~6.25%) | 4128 | 167.3 | 114.9 | 1.46x |
| 0.000 (all-coarse) | 0 | 125.4 | 84.5 | 1.48x |

### 262,144 atoms, 16,384 blocks (16 atoms/block)

| fraction | active atoms | before (us) | after (us) | speedup |
|---:|---:|---:|---:|---:|
| 1.000 (all-fine) | 262144 | 757.1 | 317.9 | 2.38x |
| 0.750 | 196640 | 660.1 | 327.4 | 2.02x |
| 0.500 | 131104 | 536.1 | 298.5 | 1.80x |
| 0.250 (~25%) | 65568 | 406.3 | 274.0 | 1.48x |
| 0.125 | 32800 | 363.9 | 254.8 | 1.43x |
| 0.0625 (~6.25%) | 16416 | 333.6 | 246.2 | 1.36x |
| 0.000 (all-coarse) | 0 | 275.4 | 209.6 | 1.31x |

Performance **improved** at every measured fraction and both sizes. The
speedup grows with active fine fraction (largest at all-fine, where the
removed bond/onsite block reductions previously covered the most threads)
and grows with system size (2.38x vs 1.92x at all-fine between the two
sizes), consistent with the removed work being genuinely proportional to
active thread count rather than a fixed per-step constant.

### Instrumented phase breakdown (262,144-atom fixture, all-fine)

```
before: atomistic_us=712.75  coarse_us=18.17  interface_us=26.84  phase_sum_us=767.82
after:  atomistic_us=292.87  coarse_us=12.28  interface_us=24.48  phase_sum_us=339.53
```

The atomistic phase (which contains the bond/onsite energy reduction) drops
59% (712.75us -> 292.87us); this phase sum is instrumented-timing (device
event synchronized per phase) and is reported separately from the headline
`setPhaseTimingEnabled(false)` number above, per the governing protocol's
instrumented-vs-uninstrumented distinction.

### Kernel launch count (average launches/step, 262,144-atom fixture)

```
fraction=1.000 (all-fine):   atomistic 10.20 -> 8.16  (-2.04)   coarse 6.12 -> 2.04  (-4.08)
fraction=0.000 (all-coarse): atomistic  4.08 -> 2.04  (-2.04)   coarse 12.24 -> 8.16 (-4.08)
```

These deltas are exactly the removed `reduceAdaptiveEnergyPartials`
(atomistic-phase and coarse-phase calls) and `finalizeAdaptiveEnergy`
launches — no other launch count changed at any fraction, confirming no
field kernel was accidentally skipped or duplicated.

### D2H bytes (analytic, not separately instrumented by this benchmark)

Each `evaluateHybrid()` call previously issued one unconditional
`double terms[8]` (64-byte) device-to-host copy regardless of whether the
caller used the result. At `adaptiveDiagnostics=0` the production step made
two such calls (initial, predictor); both are now `evaluateField()` calls
with zero D2H copies — a reduction of 128 bytes and 2 synchronous D2H
transfers per production step, on top of the launch-count and phase-time
reductions above.

## 6. Checklist mapping

See `docs/CGP_work.md` CGP-01 checklist for the item-by-item disposition.

## Commit

`CGP-01: separate adaptive field and energy evaluation`
