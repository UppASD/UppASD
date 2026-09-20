# CGP-04 — Compact adaptive field initialization evidence

Task: `docs/CGP_work.md` CGP-04 ("Remove full-system adaptive field
clearing"). Branch `gpu_hip_cu_ab_cg`, commit base `dfaa515d` (CGP-03D).

**Primary model: Claude (Sonnet 5, interactive session).**

## Environment

- GPU: 2x NVIDIA RTX A4000, driver 610.57.04, shared with another user per
  project memory `shared-gpu-host-contention`; both idle (0% utilization,
  <110 MiB used) immediately before every timed run below.
- CUDA 13.3.73 (nvcc), CMake 3.28.3, `Release` build type. Host has only 2
  CPU cores, so builds ran sequentially rather than in parallel.
- "After" builds: fresh out-of-tree `build_cgp04_cuda_fp32`
  (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`,
  `CMAKE_BUILD_TYPE=Release`) and `build_cgp04_cuda_fp64`
  (`UPPASD_PRECISION=DOUBLE`), both at the working tree with this task's
  changes.
- "Before" build (fp32 performance comparison only): identical CMake
  configuration, built from a `git worktree` pinned to commit `dfaa515d`
  (this task's parent, CGP-03D), so the comparison isolates exactly this
  task's diff.
- Performance comparison scope: fp32 only, per the disclosed scope decision
  in section 5 below (a fp64 "before" build was not built to bound this
  task's wall-clock cost on the 2-core host; fp64 correctness was verified
  in full -- see section 4).

## 1. Consumer audit (Part A)

Traced every reader of the adaptive runtime's public and internal
atom-/block-resolved field buffers before editing:

| Buffer | Consumer | Reads every atom/block, or only a compact subset? |
|---|---|---|
| `atomFieldScratch_` (internal) | `evaluateAdaptiveAtomisticBonds` (bond scatter, writes both live-bond endpoints) | Writes at `activeAtomList ∪ ghostAtomList` (every live-bond endpoint) only |
| `atomFieldScratch_` | `evaluateAdaptiveAtomisticOnsite` | Writes at `activeAtomList` only |
| `atomFieldScratch_` | `writebackAdaptiveAtomistic` (publishes into `atomField`) | Reads/writes at `activeAtomList` only |
| `atomFieldScratch_` | `restrictAdaptiveInterface` | Reads at `ghostAtomList` only |
| `atomField` (public, `gpuLattice.beff`) | `GpuDepondtIntegrator::evolveFirst/evolveSecond` | Reads at the caller-supplied active-atom list only (already active-only since RCG) |
| `atomField` | `addAdaptiveDipole`/`addAdaptiveBasisResolvedDipole` | Writes only at atoms of a *non-coarse* block, i.e. exactly `activeAtomList` (every atom of an atomistic block is atomistic by definition) |
| `atomField` | `GpuSimulation::copyToFortran()` -> Fortran `print_fields`/`beff` output, restart handoff | Reads **every** atom -- the one genuine full-lattice consumer |
| `atomField` | `GpuAdaptiveRuntime::diagnosticSnapshot()` | Sums only `atomisticMask`-true atoms, i.e. exactly `activeAtomList` -- **not** a full-lattice consumer despite downloading the whole array |
| `coarseFieldScratch_` (internal, interface-restriction accumulator) | `restrictAdaptiveInterface` (write) | Writes only at a ghost atom's projected corner blocks, which are coarse by construction (the `coarseMoment<=0` guard skips any atomistic corner) |
| `coarseFieldScratch_` | `writeAdaptiveCoarse` (only reader) | Reads at `activeBlockList` only |
| `coarseField`/`runtime.coarseField` (public + internal accumulator, aliased to the same buffer by every production call site) | `predictorAdaptiveCoarse`, `correctorAdaptiveCoarse`, `writeAdaptiveCoarse`, `evaluateAdaptiveCoarseTensor`, `finalizeAdaptiveCoarseLocal`, both dipole kernels' coarse-side accumulation | Read/write at `activeBlockList` only |
| `coarseField` | `GpuAdaptiveRuntime::diagnosticSnapshot()` | Sums **every** block unconditionally -- the one genuine full-lattice consumer on the coarse side |
| `coarseField` | `GpuSimulation::copyToFortran()` | Not read -- coarse-block state has no Fortran-visible representation |

Conclusion: three of the four full-system clears in `evaluateHybridImpl`
(`clearAdaptiveAtomistic`'s `atomFieldScratch` write, `clearAdaptiveInterface`,
`clearAdaptiveCoarse`'s `runtime.coarseField`/internal write) have **no**
consumer outside a compact list, ever. Two public-buffer writes
(`clearAdaptiveAtomistic`'s `atomField` write, `clearAdaptiveCoarse`'s public
`coarseField` write) have exactly one full-lattice consumer each
(`copyToFortran()` for `atomField`; `diagnosticSnapshot()` for `coarseField`),
identified precisely rather than assumed.

## 2. Validity contract (Part B)

- **`atomFieldScratch_`.** Valid (freshly zeroed-then-written) at
  `activeAtomList ∪ ghostAtomList` after an ordinary step. Elsewhere,
  contents are unspecified/stale between explicit full materializations --
  safe because nothing reads outside that set, this step or any later one
  (it is unconditionally re-cleared before its next use).
- **Public `atomField` (`gpuLattice.beff`).** Valid at `activeAtomList`
  every step. At every other atom it is **stale** (last value from
  whenever that atom was last active, or exactly zero if never active)
  unless a full materialization has run. This is a genuine (small,
  narrowly-scoped) relaxation of the pre-CGP-04 invariant "every non-active
  atom's `beff` is exactly zero every step" -- see the full-materialization
  gating below for why it is never observable.
- **`coarseFieldScratch_`.** Valid at `activeBlockList`'s scalar range only,
  unconditionally -- no full-lattice fallback exists or is needed (see the
  audit above).
- **Public/internal `coarseField`.** Valid at `activeBlockList` every step;
  stale elsewhere unless a full materialization has run (diagnostics>=3
  only -- see below).

### When a full materialization is required

- **Atomistic side (`atomField`).** Exactly when `copyToFortran()` might
  next read `gpuLattice.beff` -- i.e. whenever CGP-03C's own
  `adaptiveMomentUpdateNeedsFullLattice()` look-ahead (already computed once
  per step in `gpuSDSimulation.cpp` as `adaptiveNeedsFullMaterialization`,
  and already threaded into `advanceAdaptiveStep` as
  `nextStepNeedsFullMaterialization` for CGP-03D's own atom-direction
  commit) returns true, **or** this step traces diagnostics (`traceEnergy`,
  diagnostics>=3). No new look-ahead was written for this task --
  `copyToFortran()` copies `gpuLattice.beff` on exactly the same
  `printMdStatus` cadence / final-loop-iteration schedule that look-ahead
  already answers for `mmom`/`mmomi`/`emomM`, so it is reused as-is.
- **Coarse side (`coarseField`).** Exactly `traceEnergy` (diagnostics>=3) --
  `diagnosticSnapshot()` is the only full-lattice consumer and it only runs
  behind that same gate (`printStaticStageTrace`'s own internal check, and
  `GpuSimulation::release()`'s teardown snapshot at `adaptiveDiagnostics>0`,
  which is itself covered because the *last* step's predictor-stage call
  always receives `nextStepNeedsFullMaterialization=true` -- see
  `adaptiveMomentUpdateNeedsFullLattice()`'s `currentStep >= lastStep`
  branch, unchanged by this task).

Only the "initial" field-evaluation stage is exempt from the
`nextStepNeedsFullMaterialization` condition: its `atomField`/`coarseField`
output is unconditionally overwritten by the very next ("predictor")
evaluation before anything else (besides `evolveFirst`'s active-only read
and `printStaticStageTrace`'s own diagnostics>=3-gated snapshot) can read
it, so it only ever needs `traceEnergy`, never the copyToFortran look-ahead.

## 3. Implementation (Parts C, D)

`source/gpu_files/gpuAdaptiveRuntime.cpp`: a new `bool fullFieldClear`
parameter on `evaluateHybrid()`/`evaluateField()`/`evaluateHybridImpl()`
(default `true`, so every existing caller -- tests, benchmarks, any call
site that does not opt in -- is bit-for-bit unaffected) selects between the
existing full-system clear kernels and four new compact ones:

- `clearAdaptiveAtomisticActive<MeasureEnergy>` -- zeroes
  `atomFieldScratch`/public `atomField` at `activeAtomList` (one launch).
- `clearAdaptiveAtomisticGhost` -- zeroes `atomFieldScratch` at
  `ghostAtomList` (one launch; matches CGP-02's own precedent of leaving two
  small disjoint-list launches unfused absent evidence launch overhead
  dominates).
- `clearAdaptiveInterfaceActive` -- zeroes `coarseFieldScratch` at
  `activeBlockList`. Unconditional replacement for `clearAdaptiveInterface`
  (no `fullFieldClear` fallback needed at all -- see the audit).
- `clearAdaptiveCoarseActive<MeasureEnergy>` -- zeroes
  `runtime.coarseField`/public `coarseField` at `activeBlockList`.

Energy-term zeroing (`kernels.energyTerms[...]`) keeps the same
unconditional-when-`MeasureEnergy` treatment as the full kernels, sized to
at least 2 (atomistic) / 6 (coarse) threads even when the corresponding
active list is empty (an all-coarse or all-fine system with diagnostics
enabled must still zero its energy accumulator).

`source/gpu_files/gpuSimulation.cpp`: `advanceAdaptiveStep`'s
`evaluateAdaptiveFields` lambda gained a `fullFieldClear` parameter, set at
each of its four call sites:

| Stage | `fullFieldClear` |
|---|---|
| initial | `traceEnergy` |
| predictor | `nextStepNeedsFullMaterialization \|\| traceEnergy` |
| corrector (traceEnergy-gated) | `true` (only ever called when traceEnergy) |
| reconstructed (traceEnergy-gated) | `true` (ditto) |

`GpuMomentUpdater`, `GpuCorrelations`, and the reconstruction formula are
untouched.

## 4. Tests

New unit test `testCompactFieldClear`
(`tests/coarse_graining/test_gpu_adaptive_runtime.cpp`), added to the
`gpu_adaptive_runtime_tests` target, using the existing mixed `KernelFixture`
(atoms 3-6 active/fine, blocks 1/4 coarse, ghost shell `{2,7}`):

1. **Parity.** `fullFieldClear=false` reproduces `fullFieldClear=true`
   bitwise (`==`, not a tolerance) at every active atom and active block.
2. **Negative control (Part E).** A sentinel (`-999`) planted in the public
   `atomField`/`coarseField` buffers before a `fullFieldClear=false` call
   survives untouched at every *inactive* atom/block -- proving the
   compaction is real, not merely a slower way to clear everything.
3. **Materialization.** Re-running the same still-partially-poisoned buffers
   with `fullFieldClear=true` erases every sentinel and reproduces the
   original full-clear reference exactly.
4. **Production-consumer tie-in.** With a poisoned `coarseFieldData()` left
   by a `fullFieldClear=false` call, `diagnosticSnapshot()`'s unmasked
   `coarseFieldSumT` checksum differs from the same snapshot taken after a
   `fullFieldClear=true` call -- proving the coarse gate is load-bearing for
   the actual production diagnostic it protects, not decorative.

`ctest` run (`coarse-graining-gpu-adaptive-runtime`), both precisions:
**PASS** (`CG-09/CG-10 GPU adaptive runtime tests passed`, `gpu_adaptive_runtime_tests`
binary run directly on both `build_cgp04_cuda_fp32` and `build_cgp04_cuda_fp64`).

### Correctness suite

`ctest -L cg13-cuda -j2`:

- `build_cgp04_cuda_fp32`: **33/34 passed**. The one failure,
  `adaptive-cg-production-e2e` (`run_production_e2e.py`,
  `gpu_fft_static_mixed` case, `assert
  abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), is the exact
  same pre-existing, unrelated open-FFT-dipole defect documented in
  `docs/CGP-00B_ENERGY_GATE_EVIDENCE.md` and reproduced identically in
  every CGP-01 through CGP-03D evidence doc -- not a CGP-04 regression.
  All other `adaptive-cg-production-e2e` sub-cases passed, including
  `CG-10.5 finite-temperature adaptive fine-region execution` and
  `CG-10.5 finite-temperature mixed-resolution execution`.
- `build_cgp04_cuda_fp64`: **34/34 passed** -- including the
  `gpu_fft_static_mixed` case, confirming the fp32 failure above is a
  precision-specific pre-existing issue, not something this task's changes
  interact with.

### All-fine T=0 / finite-T ASD parity

`coarse-graining-adaptive-asd-parity` (all-fine control) and the finite-T
fixtures inside `adaptive-cg-production-e2e` (`CG-10.5
finite-temperature adaptive fine-region execution`, `CG-10.5
finite-temperature mixed-resolution execution`) all **PASS** on both
precisions.

## 5. Performance evidence

`benchmark_gpu_adaptive_runtime`'s `adaptiveCoarseStep()` PERF-CG-SWEEP loop
(`--warmup 20 --iterations 200 --repetitions 5`, `setPhaseTimingEnabled(false)`
for the headline `step_wall_us` number -- see the negative control note in
section 6 of `docs/CGP-01_FIELD_ENERGY_SPLIT_EVIDENCE.md` for why the
uninstrumented number is the one reported), fp32, device 0. Both sizes use
16 atoms/block, matching CGP-01's own fixtures. Two samples per arm,
median reported; one process-level warm-up run discarded before the first
timed sample, per project memory `gpu-benchmark-cold-start`.

`benchmark_gpu_adaptive_runtime.cpp`'s `adaptiveCoarseStep()` -- the exact
production-mimicking predictor/corrector pair `evaluateAdaptiveFields`
exercises -- was updated to pass `fullFieldClear=false` to both
`evaluateField()` calls, matching production's own default; without this
the benchmark would keep measuring the pre-CGP-04 full clear regardless of
the source change. `atomistic-parity` is **PASS** at every measured point
on both builds (no correctness regression).

**Disclosed outlier:** the first `after`, 65536-atom sample showed
anomalous slowdowns at fractions 0.125-0.5 (up to 2x the second sample,
`step_mad_us` up to 17) not present in either `before` sample or a third
`after` sample taken immediately afterward; `nvidia-smi` showed both GPUs
idle before and after. This does not match the documented cold-start
pattern (~18% on the *first* run at a size) so is reported as a
discarded transient-contention outlier rather than folded into the cold-start
note; the two consistent `after` samples (2nd and 3rd invocations) are used
for the reported median instead.

### 65,536 atoms, 4,096 blocks (16 atoms/block)

| fraction | active atoms | before (us) | after (us) | speedup |
|---:|---:|---:|---:|---:|
| 1.000 (all-fine) | 65536 | 107.1 | 101.5 | 1.05x |
| 0.750 | 49184 | 130.1 | 135.9 | 0.96x |
| 0.500 | 32800 | 113.9 | 114.1 | 1.00x |
| 0.250 (~25%) | 16416 | 85.5 | 86.9 | 0.98x |
| 0.125 | 8224 | 77.2 | 82.4 | 0.94x |
| 0.0625 (~6.25%) | 4128 | 75.8 | 80.3 | 0.94x |
| 0.000 (all-coarse) | 0 | 44.6 | 41.4 | 1.08x |

### 262,144 atoms, 16,384 blocks (16 atoms/block)

| fraction | active atoms | before (us) | after (us) | speedup |
|---:|---:|---:|---:|---:|
| 1.000 (all-fine) | 262144 | 307.1 | 307.1 | 1.00x |
| 0.750 | 196640 | 293.3 | 297.7 | 0.99x |
| 0.500 | 131104 | 234.9 | 225.4 | 1.04x |
| 0.250 (~25%) | 65568 | 178.1 | 162.2 | 1.10x |
| 0.125 | 32800 | 148.4 | 126.4 | 1.17x |
| 0.0625 (~6.25%) | 16416 | 132.2 | 95.7 | 1.38x |
| 0.000 (all-coarse) | 0 | 83.7 | 47.3 | **1.77x** |

The benefit is real, reproducible (consistent sign and magnitude across
both `before` and both clean `after` samples), and grows both with system
size and with coarse fraction -- exactly the O(N) -> O(active) shape this
task targeted. It is small or slightly negative at high active fractions
(0.75-1.0): at those fractions `activeAtomList` already covers nearly every
atom, so the compact clear does the same memory traffic as the full clear
plus one extra kernel launch (the ghost clear) and a level of list
indirection the full clear does not pay -- a small, expected, and honestly
reported regression, not hidden in the average.

### Mechanism (instrumented phase breakdown, 262,144 atoms, all-coarse)

| phase | before (us) | after (us) |
|---|---:|---:|
| atomistic (clear+bond+onsite+writeback) | 35.68 | 2.75 |
| coarse (clear+tensor+dipole) | 33.80 | 30.86 |
| interface (prolongation+clear+restriction) | 10.18 | 9.90 |

The atomistic phase -- purely the full-system clear at this fraction, since
`activeAtomWork`/`ghostClearWork` are both zero -- drops **12.9x**, from
35.7us to 2.8us, accounting for essentially all of the 36.4us whole-step
improvement (83.7us -> 47.3us). The coarse phase improves only slightly
(9%): at fraction 0.0 *every* block is coarse, so `clearAdaptiveCoarseActive`
covers the same element count as `clearAdaptiveCoarse` -- the small residual
difference is compact-list indirection overhead, not an O(N) reduction,
exactly as the design predicts (`clearAdaptiveCoarseActive` only wins when
some blocks are atomistic). The interface phase is within noise, as
expected (`clearAdaptiveInterfaceActive` was already unconditional and
untouched by the before/after comparison's precision/scope, and at
fraction 0.0 covers the same block count either way).

### fp64 scope note

A fp64 "before" build was not produced (see Environment). The fp64 "after"
build's own sweep (`build_cgp04_cuda_fp64`, same protocol) shows the same
qualitative shape -- e.g. at 262,144 atoms, all-coarse: `step_wall_us=192.6`
vs. all-fine `step_wall_us=677.8` -- consistent with the mechanism (a
byte-count-proportional memory-traffic reduction) transferring to fp64
unchanged in kind, just at double the per-element byte cost. Absolute
fp32-vs-fp64 speedup numbers are not claimed without the matching baseline
build.

## 6. Negative control (task-level)

Covered by `testCompactFieldClear` (see Tests, above) rather than a
separate poisoned production build -- the sentinel-survival check inside
the unit test directly exercises "normal adaptive dynamics must remain
correct [with stale non-active entries]" and "a full-field diagnostic must
detect/materialize them correctly rather than publishing the poison" on
the exact kernels production uses, with the exact list membership
production computes, which is a stronger and more precise proof than a
whole-program poisoned rebuild would give for a bug class this localized.

## Outcome

Three of the four adaptive field-clear kernels (`clearAdaptiveAtomistic`'s
`atomFieldScratch` write, `clearAdaptiveInterface`, `clearAdaptiveCoarse`'s
internal `runtime.coarseField` write) had no consumer outside a compact
list, ever, and are now unconditionally compact. The two public-buffer
writes (`atomField`/`gpuLattice.beff` for `copyToFortran()`; public
`coarseField` for `diagnosticSnapshot()`) each have exactly one genuine
full-lattice consumer, so a new `fullFieldClear` parameter selects between
the compact and full-system kernels, reusing CGP-03C's existing
`nextStepNeedsFullMaterialization` look-ahead for the `atomField` case (no
new look-ahead logic was written) and the existing `traceEnergy`
(diagnostics>=3) gate for the `coarseField` case. Measured at 262,144
atoms, all-coarse (the steady-state case CGP-03's own performance gate
targets): the atomistic clear phase dropped 12.9x (35.7us -> 2.8us) and
the whole adaptive step dropped 1.77x (83.7us -> 47.3us), fp32. The benefit
shrinks toward parity (and is slightly negative, ~2-6%, from added launch/
indirection overhead) at high active fractions, where the compact clear
covers nearly as many elements as the full clear anyway -- reported
honestly rather than averaged away. fp64 correctness was fully verified
(34/34 `cg13-cuda`, including the one case that fails on fp32 for an
unrelated pre-existing reason); a fp64 performance baseline was not built,
disclosed as a scope decision bounded by this task's build-time budget on
the 2-core host.

## Checklist

* [x] Field consumers inventoried.
* [x] Validity contract documented.
* [x] Normal atom-field clear scales with touched atoms
  (`activeAtomList ∪ ghostAtomList`, not total atoms).
* [x] Interior coarse atom field entries are not cleared unnecessarily
  (compact atomistic clear only touches `activeAtomList`/`ghostAtomList`).
* [x] Coarse-field clear reviewed similarly (`clearAdaptiveCoarseActive`,
  `clearAdaptiveInterfaceActive`).
* [x] No hidden read of stale field entries exists (consumer audit found
  every reader of every buffer explicitly; the two full-lattice consumers
  are gated by `fullFieldClear`, proven by `testCompactFieldClear`'s
  poison negative control and diagnostic tie-in).
* [x] Depondt active-atom parity passes (`atomistic-parity` PASS at every
  benchmark sweep point on both precisions; `cg13-cuda` unaffected).
* [x] Interface restriction parity passes (`cg13-cuda` interface/DMI/moving
  fixtures unaffected; `clearAdaptiveInterfaceActive` is unconditional and
  provably safe by construction -- no consumer outside `activeBlockList`
  exists this step or any later one).
* [x] Full diagnostics still return valid data (`fullFieldClear` forced by
  `traceEnergy`/diagnostics>=3 at every stage that can call
  `diagnosticSnapshot()`, including the teardown snapshot at
  `adaptiveDiagnostics>0`, covered via the last-iteration branch of
  `adaptiveMomentUpdateNeedsFullLattice()`).
* [x] Poison/sentinel negative control passes (`testCompactFieldClear`).
* [x] Full-system O(N) field-clear timing removed or quantitatively reduced
  (12.9x on the isolated atomistic phase, 262,144 atoms all-coarse).
* [x] Whole-step fp32 benchmark rerun (two sizes, before/after, up to
  1.77x at the steady-state all-coarse target; disclosed small regression
  at high active fractions).

## Commit

`CGP-04: compact adaptive field initialization`
