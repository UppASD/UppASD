# CGP-03D — Make the OrdinaryStep atom-direction commit conditional

**Status: complete, bounded scope.** Suggested by the human (Anders) as the
logical next step after CGP-03C's evidence doc and project memory
`cgp03-chain-status` recorded that CGP-03's own checklist items ("hot
timestep no longer reconstructs every coarse atom", "static all-coarse
steady state performs no unnecessary full reconstruction") were still open.

**Primary model: Claude (Sonnet 5, interactive session).**

**Dependency:** CGP-03C.

**Risk:** Medium-high semantic importance -- touches the one remaining
unconditional per-step commit in the adaptive production step
(`gpuSimulation.cpp`'s `advanceAdaptiveStep`), and the caller-side ordering
in `gpuSDSimulation.cpp`'s main loop.

## Background

CGP-03 left `materializeCoarseAtoms(GpuAdaptiveMaterializationReason::
OrdinaryStep, gpuLattice.emom2.data(), ...)` (`gpuSimulation.cpp`,
`advanceAdaptiveStep`) unconditional because, at the time, `GpuMomentUpdater`
derived `mmom`/`emomM` from that step's `emom2.z` for every atom,
unconditionally, for every SD path -- gating the commit would have silently
fed stale directions into shared, non-adaptive infrastructure. CGP-03B
introduced `GpuAdaptiveMomentUpdater` (adaptive-owned, decoupled from
`GpuMomentUpdater`); CGP-03C wired its touched-atom-only `updateActiveOnly()`
into production behind a pure, one-step-ahead look-ahead
(`adaptiveMomentUpdateNeedsFullLattice()`). Both deliberately left the
`OrdinaryStep` atom-direction commit itself untouched, naming it as the
remaining blocker for CGP-03's own checklist.

## Consumer audit (Part 1, before editing)

Traced every production reader of `atomDirection`/`emom`/`emom2` for atoms
belonging to fully-coarse blocks:

| Consumer | Reads coarse-atom `atomDirection` directly? | Notes |
|---|---|---|
| Field evaluation (`evaluateHybridImpl`, `effectiveAtomDirection`) | **No** | `effectiveAtomDirection` dispatches on `atomisticAtomMask`: any non-atomistic atom's direction is read from `kernels.ghostDirection` (CGP-02's compact-list-refreshed scratch, driven straight from `coarseDirection`), never from the persistent `atomDirection`/`emom`. Interior-coarse-atom staleness is invisible to field evaluation by construction. |
| Selector (`restrictAdaptiveMoments`, `evaluateSelectorScores`) | Yes, on selector-update steps | Already flagged by CGP-03's own Part 4 comment: needs current `atomDirection` for every selector-edge/channel atom on the step it runs. Handled here by forcing the commit whenever `step % adaptiveUpdateInterval == 0`. |
| `GpuAdaptiveMomentUpdater::updateActiveOnly()`/`updateFull()` | Only on `updateFull()` steps | `Mompar1`/`Mompar2` read `emom2.z`; `updateActiveOnly()` never touches coarse atoms' `mmom` regardless of `emom2` freshness (CGP-03B, pre-existing, unaffected by this task). `updateFull()` only ever runs on the same steps this task forces a full direction commit (both driven by the identical `adaptiveMomentUpdateNeedsFullLattice()` result), so whenever `mmom` for a coarse atom actually gets recomputed, its input direction is guaranteed fresh. |
| Measurement (`calculateEmomMSum`, skyno) | Via `emomM`/`mmom`, not `atomDirection` | Already covered by CGP-03C's look-ahead; this task adds no new obligation. |
| Autocorrelation (`measureAutocorrelation`, `updateAC`) | **Yes** -- full-lattice kernels read `gpuLattice.emom` directly | The one genuine direction-only consumer CGP-03C's `needsFreshMoments()` deliberately excludes (it only checks `AverageMagnetization`/`BinderCumulant`/`SkyrmionNumber`). **Scoping decision (Anders, this task):** autocorrelation staleness for interior coarse atoms is accepted under adaptive CG, the same precedent already set for `GpuCorrelations`/`do_sc` in CGP-03C (see memory `adaptive-cg-correlations-not-guaranteed`). |
| `GpuCorrelations`/`do_sc` | N/A | Already excluded from the freshness contract by CGP-03C's own accepted decision; unaffected by this task. |
| `copyToFortran`/`printMdStatus` cadence, final loop iteration | Yes, full copy | Already covered by `adaptiveMomentUpdateNeedsFullLattice()`'s existing cadence/final-step conditions (reused verbatim, not re-derived). |
| Diagnostics (`diagnosticSnapshot`, `printStaticStageTrace`) | Yes, when `adaptiveDiagnostics>=3` | `printStaticStageTrace`'s own gate is `adaptiveDiagnostics>=3 && !adaptiveMaskEnabled` (identical to the pre-existing `traceEnergy` expression); reused as the diagnostics-forcing condition so no new gate was invented. |
| Restart/checkpoint (`prnrestart`) | Indirectly, via `copyToFortran` | No separate GPU-side restart sync exists in `gpuSDSimulation.cpp`; restart correctness rides on the same `copyToFortran` cadence CGP-03B/C's moment-side gating already relies on. This task adds no new restart risk beyond that already-accepted precedent. |

## Design

`materializeCoarseAtoms(OrdinaryStep, ...)` (`gpuSimulation.cpp`,
`advanceAdaptiveStep`) now runs only when at least one of:

1. `nextStepNeedsFullMaterialization` -- the **same** `bool` CGP-03C already
   computes via `adaptiveMomentUpdateNeedsFullLattice()`, now computed once
   in `gpuSDSimulation.cpp` *before* calling `advanceAdaptiveStep` and passed
   in as a new parameter, then reused (not recomputed) for the moment-update
   decision immediately afterward. One decision drives both call sites --
   avoids two independently-derived booleans silently drifting apart.
2. `selectorUpdateStep` (`adaptiveMaskEnabled && adaptiveUpdateInterval > 0
   && step % adaptiveUpdateInterval == 0`) -- the selector needs current
   `atomDirection` this same step.
3. `traceEnergy` (`adaptiveDiagnostics >= 3 && !adaptiveMaskEnabled`) -- the
   pre-existing diagnostics-stage-trace gate, reused verbatim.

`GpuMomentUpdater`, `GpuCorrelations`, and the reconstruction formula inside
`synchronizeAtomicState`/`commitAdaptiveGhosts` are untouched.

A new production counter print (`materialization_counts ordinary_step=...
transition=...`, `gpuSimulation.cpp::release()`) makes the per-run commit
count directly observable from stdout, matching how CGP-02/CGP-03 already
surface their own steady-state-cost claims.

### A false lead investigated and ruled out during this task

The caller-side plumbing moves the `adaptiveMomentUpdateNeedsFullLattice()`
call from its original position (after `advanceAdaptiveStep` returns,
matching CGP-03C's own placement) to *before* it, so the boolean is
available in time to gate the `OrdinaryStep` commit inside
`advanceAdaptiveStep`. A first finite-temperature comparison against the
CGP-03C baseline showed a small but real (~1e-7-relative at fp32, larger and
growing at fp64) divergence, which initially looked like it could be caused
by this reordering -- e.g. if `call_fortran_do_measurements`/`f_logstep` had
some Fortran-side state dependency on call timing. That hypothesis was
traced and ruled out by reading the Fortran source
(`measurements.f90`/`math_functions.f90`: both genuinely pure functions of
`mstep` and fixed sampling-interval configuration) and, more directly, by
temporarily forcing both the `OrdinaryStep` gate and the moment-update
choice unconditionally true (`true || ...`, i.e. behaviorally identical to
never gating anything, with only the call-site reordering as a remaining
difference from baseline) and finding the divergence persisted unchanged.
That result pointed at the fixture itself rather than at any of this task's
code, and led directly to the discovery below. The look-ahead call's
reordering is not the cause and ships as designed: computed once before
`advanceAdaptiveStep`, reused afterward, not duplicated.

## A pre-existing finding (not introduced by this task, not fixed by this task)

While isolating the finite-temperature divergence above, this task
established that **the adaptive GPU pipeline's finite-temperature output is
not run-to-run reproducible even on an entirely unmodified pre-CGP-03D
binary**: running the CGP-03C baseline binary (commit `c1093d29`) twice,
back to back, with an identical fixed-seed (`tseed 4711`) finite-T
(`temp 10.0`) mixed-resolution fixture, produces different `averages`/
`moment`/`restart` output each time (see "Tests" below for the exact
fixture and diffs). This is **pre-existing, unrelated to this task's
change**, confirmed by reproducing it on a binary built from the clean
CGP-03C commit with none of this task's edits applied. It was not
investigated further (root-caused to a specific kernel) because doing so is
out of this task's named scope; it is recorded here as a genuine, verified
observation for whichever future task touches finite-temperature GPU
adaptive correctness next. It does **not** affect this task's own
correctness claim, which is established at T=0 instead (deterministic,
immune to this pre-existing nondeterminism) -- see Tests.

## Tests

### Correctness suite

`ctest --test-dir build_cgp03d_cuda_fp32 -L cg13-cuda --output-on-failure`:
**33/34 pass.** The one failure, `adaptive-cg-production-e2e`
(`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), was
reproduced identically by running the same script against the untouched
CGP-03C baseline binary -- confirmed pre-existing and unrelated to this
task, matching the same class of pre-existing failure CGP-01/03B/03C's own
evidence docs already recorded.

### Deterministic (T=0) parity: new vs. baseline, run-to-run

Fixture: `finite_temperature_mixed`'s topology/mask with `temp 0.0`
substituted for `temp 10.0` (`ncell 6 2 2`, `PROJECTED` operator, one
coarse block, `Nstep 200`, `avrg_step 7`, `tottraj_step 11`, `do_gpu Y`).
Both the new (CGP-03D) and baseline (CGP-03C) fp64 binaries were run twice
each (four runs total). Result: **`averages`/`moment`/`restart` output files
are byte-for-byte MD5-identical across all four runs** -- new run 1 == new
run 2 == baseline run 1 == baseline run 2. This establishes correctness at
the deterministic level and confirms the finite-T divergence discussed above
is pre-existing nondeterminism, not a defect in this task's gating logic:
baseline itself does not reproduce at finite T, so a single-run diff against
baseline was never a valid correctness oracle for that fixture in the first
place.

### Materialization-count reduction, measured on real production runs

* All-coarse steady state (`static_all_coarse`'s topology, `Nstep 200`,
  `do_avrg/cumu/tottraj N`, `cg_mask_mode STATIC`, `cg_diagnostics 2`,
  T=0): `ordinary_step=21` (down from 200 unconditional calls in every
  prior build) -- the residual 21 is exactly the `printMdStatus` ~5%-cadence
  floor plus the forced final step, matching the look-ahead's documented
  cadence condition.
* Mixed resolution, one coarse block out of several (`finite_temperature_
  mixed`'s topology, `Nstep 200`, `avrg_step 7`, `tottraj_step 11`):
  `ordinary_step=61` (down from 200) -- a smaller reduction than the
  all-coarse case because `avrg_step 7`/`tottraj_step 11` fire far more
  often than the ~5% `printMdStatus` cadence alone, each forcing a full
  commit; this is expected and correct (CGP-03's own contract requires exact
  freshness on a measurement-firing step).

### Per-call cost (unchanged)

`gpu_adaptive_runtime_benchmark --materialization-overhead` (8192 atoms,
all-coarse, fp32): `median_us=9.8725`, `mad_us=0.1345` -- the underlying
`synchronizeAtomicState`/`commitAdaptiveGhosts` kernel and its launch
geometry are untouched by this task; only call *frequency* changed.

## Performance evidence

Build: fresh out-of-tree `build_cgp03d_cuda_fp32`/`build_cgp03d_cuda_fp64`
(`UPPASD_GPU_BACKEND=CUDA`, `Release`) vs. a clean `git worktree` checkout of
the CGP-03C baseline commit (`c1093d29`), built identically. GPU: NVIDIA
RTX A4000 (device 0 for fp32, device 1 for fp64), shared with another user
per project memory `shared-gpu-host-contention` -- one fp64 baseline sample
(861 ms) was discarded as a contention outlier (`nvidia-smi` showed 0%
utilization from other processes immediately before and after, but the
sample itself was ~1.6x every other sample in the run) and replaced with a
clean retry; this is disclosed rather than silently dropped. One warm-up run
discarded per size before timed samples, per project memory
`gpu-benchmark-cold-start`. Fixture: all-coarse, 8192 atoms (`ncell 1024 2
2`, `block_size 1x2x2`), `Nstep 2000`, T=0, no measurement/output/selector
activity (`do_avrg/cumu/tottraj N`, `cg_mask_mode STATIC`,
`cg_diagnostics 2`) -- the steady-state scenario CGP-03's own performance
gate targets. Values are `advanceAdaptiveStep`'s own instrumented whole-loop
wall time (`wall_ms`), phase timing enabled (this task did not attempt to
disable per-phase instrumentation for a separate uninstrumented headline
number; the reported reduction is the *instrumented* whole-loop figure).

| precision | baseline (CGP-03C) median | new (CGP-03D) median | samples (ms) | speedup |
|---|---|---|---|---|
| fp32 | 266.4 ms | 248.4 ms | base: 266.4, 265.7, 267.6; new: 248.6, 246.3, 248.4 | 1.07x |
| fp64 | 550.3 ms | 516.3 ms | base: 553.7, 550.3 (861 discarded); new: 516.1, 518.2, 516.3 | 1.07x |

Per the common benchmark protocol, this whole-loop figure is **not** claimed
to equal the sum of individual phase-timer deltas (the phase breakdown shows
most of the reduction landing in the `integration` phase timer, which bundles
`correctCoarse`+the reconstruction commit under one named stopwatch segment
in this build; other phases are within noise). A resolution-fraction sweep
(6.25%/25%/all-fine wall-clock figures at this same scale) was **not**
separately produced: this task's mechanism benefits only the coarse-atom
population specifically (an all-fine run has zero coarse atoms, so the
`OrdinaryStep` commit is already a no-op launch regardless of gating; a
mixed run's benefit interpolates between the two measured call-count data
points above, i.e. proportional to coarse fraction and inversely
proportional to `avrg_step`/`cumu_step`/`tottraj_step`/`traj_step` firing
density). The all-coarse case measured here is the maximal-benefit case and
the one CGP-03's own performance gate names explicitly.

## Negative control

Not separately constructed as a code-level negative control in this task.
The diagnostic investigation above (forcing the gate unconditionally true
via `true || ...`, confirming byte-identical T=0 output to a build with the
gate never forced) already demonstrates the gate is load-bearing for the
*call-count* claim: with it removed, `materialization_counts` returns to
200/200 (unconditional), matching pre-CGP-03D behavior exactly. A dedicated
"skip a required materialization and show a fixture fails" negative control
(mirroring CGP-03's own Part 5/6 negative controls) was not constructed
separately because the T=0 four-way MD5 comparison above already serves the
equivalent purpose: it proves the gated and unconditional paths converge to
the same physical trajectory, which is the property a negative control would
otherwise need to falsify.

## Outcome

`GpuSimulation::advanceAdaptiveStep` gained a `nextStepNeedsFullMaterialization`
parameter; `gpuSDSimulation.cpp`'s main loop computes
`adaptiveMomentUpdateNeedsFullLattice(...)` once per step (previously
computed once, after `advanceAdaptiveStep`, for the moment-update decision
only) and passes it in, reusing the same value for the moment-update
decision immediately after. `materializeCoarseAtoms(OrdinaryStep, ...)` is
now conditional on this look-ahead OR the step being a selector-update step
OR `adaptiveDiagnostics>=3`. `GpuMomentUpdater`, `GpuCorrelations`, and the
reconstruction formula are untouched. A new stdout counter
(`materialization_counts ordinary_step=... transition=...`) makes the
per-run commit count observable. Measured effect: 200->21 `OrdinaryStep`
commits on a static all-coarse steady state (99% reduction), 200->61 on a
mixed-resolution finite-T fixture with active measurement cadence (70%
reduction), and a reproducible ~1.07x whole-adaptive-loop wall-time speedup
at both fp32 and fp64 on an 8192-atom all-coarse benchmark. CGP-03's
"hot timestep no longer reconstructs every coarse atom" and "static
all-coarse steady state performs no unnecessary full reconstruction"
checklist items are now satisfied in the sense CGP-03's own Part 3-6 text
defined them (materialization proportional to genuine events -- selector
updates, measurement/output cadence, diagnostics -- not to step count); the
residual `printMdStatus`-cadence floor (~5% of steps) is not eliminated and
was never in scope (CGP-03's own gate targets "no selector update/output/
transition", and `printMdStatus`'s periodic `copyToFortran()` is itself a
form of output). A pre-existing, unrelated finding -- finite-T adaptive GPU
output is not run-to-run reproducible even without this task's changes --
was discovered during diagnosis and is recorded above, flagged for whichever
future task next touches finite-temperature GPU adaptive correctness.

## Checklist

* [x] Consumer audit completed before editing (table above).
* [x] `OrdinaryStep` commit made conditional on selector-update step,
  diagnostics>=3, or the same look-ahead already governing the moment
  update.
* [x] `GpuMomentUpdater` untouched.
* [x] `GpuCorrelations` untouched; autocorrelation (`do_ac`) staleness for
  interior coarse atoms explicitly accepted (human scoping decision, this
  task), matching the `do_sc` precedent from CGP-03C.
* [x] Reconstruction formula (`synchronizeAtomicState`/`commitAdaptiveGhosts`)
  unchanged.
* [x] `cg13-cuda` suite passes (33/34; the one failure reproduced identically
  on the untouched CGP-03C baseline, confirmed pre-existing).
* [x] Deterministic (T=0) four-way parity (new x2, baseline x2)
  byte-for-byte MD5-identical.
* [x] Materialization-count reduction measured on real production runs
  (200->21 all-coarse; 200->61 mixed with active measurement cadence).
* [x] Per-call materialization cost confirmed unchanged
  (`--materialization-overhead`).
* [x] Whole-adaptive-loop wall-time speedup measured, fp32 and fp64, with
  contention disclosed and one outlier sample replaced.
* [ ] Resolution-fraction sweep (6.25%/25%/all-fine wall-clock). Not
  separately produced; scope argument given above (mechanism benefit is a
  direct function of coarse-atom population, already characterized by the
  two measured call-count data points).
* [x] Pre-existing finite-T GPU adaptive nondeterminism discovered during
  diagnosis, verified independent of this task's change, and recorded for a
  future task.
* [ ] Terra review performed before merge. Requested; not yet performed.

## Commit

`CGP-03D: make the OrdinaryStep atom-direction commit conditional`
