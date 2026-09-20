# CGP-03 — Lazy coarse-atom materialization: audit and design

Task: `docs/CGP_work.md` CGP-03 ("Define and implement lazy coarse-atom
materialization"). Branch `gpu_hip_cu_ab_cg`, commit base `2ea4277e` (CGP-02).

Status: **Complete, bounded scope** (Human-selected option from "Scope
decision needed" below: the bounded implementation, not the full invasive
one). Parts 1-2 are the audit/contract below. Parts 3-6 are implemented to
the extent the audit shows is safe without also making `GpuMomentUpdater`
adaptive-aware and adding measurement-schedule look-ahead (§2, §5); that
larger change is explicitly deferred, not attempted, matching the CGP-00B
precedent of a documented partial result plus an escalation clause rather
than an under-tested full implementation.

## 0. What is reconstructed today, and where

`GpuAdaptiveRuntime::synchronizeAtomicState(atomDirection, policy)`
(`source/gpu_files/gpuAdaptiveRuntime.cpp:3374-3436`) is the sole routine
that writes atom-resolved directions for non-atomistic atoms. Under the
default `Aligned` policy it launches `commitAdaptiveGhosts`
(`gpuAdaptiveRuntime.cpp:795-820`) over `atoms_ * ensembles_` threads; every
thread whose atom is not atomistic (`!runtime.atomisticAtomMask[atom]`)
overwrites `atomDirection[atom]` with its owning block's
`coarseDirection[block]` — an unconditional full-lattice write every time
it is called, regardless of how many atoms actually changed.

It is called from exactly two sites, both in
`GpuSimulation::advanceAdaptiveStep` (`source/gpu_files/gpuSimulation.cpp`):

* **line 1252** — unconditionally, once per production step, right after
  the Depondt corrector and `correctCoarse`.
* **line 1266** — conditionally, only inside
  `if(adaptiveMaskEnabled && adaptiveUpdateInterval > 0 && step %
  adaptiveUpdateInterval == 0)` (`:1257-1268`), after a selector
  update/publish cycle.

So today: **every** production step pays one full-lattice reconstruction
(line 1252); selector-update steps pay a second one (line 1266). The line
1266 call is legitimately needed whenever a transition was actually
published (newly-fine blocks need real reconstruction, newly-coarse blocks
need their member atoms flipped to track `coarseDirection`); it is currently
unconditional on the interval, not on whether a transition actually
happened.

## 1. Consumer audit

Traced by reading `gpuAdaptiveRuntime.cpp`/`gpuSimulation.cpp` directly, and
by two focused read-only sub-agent traces of (a) `GpuMomentUpdater`,
measurement/correlation and restart on the GPU/C++ side, and (b) the CPU
Fortran adaptive path (`source/CoarseGraining/adaptivecgproduction.f90`,
`adaptivehybridsolver.f90`) and Fortran restart. Every claim below is a
direct code citation, not an inference from documentation.

| Consumer | Every step? | Needs *every* coarse atom? | Can use block state? | Can materialize lazily? |
|---|---|---|---|---|
| Hybrid field assembly, live/interface bonds (`evaluateAdaptiveAtomisticBonds`, `restrictAdaptiveInterface`) | Every step | No — only ghost-shell (CGP-02 compact list) | Yes, interior coarse bonds use `evaluateAdaptiveCoarseTensor` on `coarseDirection` directly, never `atomDirection` | Already effectively lazy (CGP-02); no interior-coarse read at all |
| Selector edge scoring (`selectorAdaptiveScores`, `gpuAdaptiveRuntime.cpp:337-363`) | Only on `adaptiveUpdateInterval` steps (`gpuSimulation.cpp:1260`) | Only for atoms that are selector-edge endpoints, which is a fixed edge set over the *topology*, not filtered by ghost-shell membership — can include interior-coarse atoms if a selector edge crosses deep into a coarse block | No — needs the literal per-atom direction for both edge endpoints | Yes — materialize immediately before this call, gated by the same `adaptiveUpdateInterval` condition already in place |
| Channel/moment restriction (`restrictAdaptiveMoments`, `gpuAdaptiveRuntime.cpp:249-302`) | Only on `adaptiveUpdateInterval` steps (`restrictMoments()` at `gpuSimulation.cpp:1259`) | **Fully-coarse blocks: no** — the loop computes `resultant` from `atomDirection` for every block-member atom regardless of state, but for `blockState==coarseState` the function discards `resultant` and uses `momentSum * coarseDirection` instead (`:283-288`, "return" before `resultant` is ever used) — the atomDirection read for those atoms is dead. **Fine/buffer blocks: yes**, but this is a non-issue: block membership is atomistic/non-atomistic *as a whole block*, never mixed within one block (see next paragraph), and fine/buffer-block atoms are always authoritative every step regardless of adaptive commit. | Yes for coarse blocks; N/A for fine/buffer blocks (already fresh) | N/A beyond the existing gate |
| Polarization gate (`evaluateAdaptivePolarizationGate`) | Only on update-interval steps | Pure function of `restrictAdaptiveMoments`' own outputs (`channelMomentSum`/`coarseMoment`) | Yes, no direct atom read | N/A (already indirect) |
| Transition proposal/dilation (`proposeAdaptiveState`, `dilateAdaptiveState`) | Only on update-interval steps | No — operate on `selectorScores`/`pendingState`, block-indexed only | Yes | N/A |
| Transition publication (`publishAdaptiveState`, `gpuAdaptiveRuntime.cpp:546+`) | Only for blocks whose state actually changes, on update-interval steps | Only writes/reads atoms in *its own block's* CSR range (one thread per block; coarse→fine reconstruction reads `ghostDirection`/writes `atomDirection` for that block only) | Partially — already block-scoped, not lattice-wide | Already effectively per-block lazy; only needs the block's own prior `ghostDirection`/`coarseDirection` to be current, not the whole lattice |
| Constrained-cone reconstruction (`commitAdaptiveGhosts` under `ConstrainedCone`) | Every call to `synchronizeAtomicState` | Yes, by design (CGP-02 evidence doc §2: `ConstrainedCone` commits an interpolated direction for *every* non-atomistic atom, not only the ghost shell) | No | Not under `ConstrainedCone` without changing accepted physics — out of scope per CGP-02 precedent |
| Dipole, non-coarse-block branch (`addAdaptiveDipole`, `gpuAdaptiveRuntime.cpp:1378-1427`, and `addAdaptiveBasisResolvedDipole`) | Every step this kernel runs (i.e. every field evaluation when FFT/open dipole is active — **not** update-interval-gated) | **Fully-coarse blocks (`runtime.coarseBlockMask[block]`): no** — uses `coarseDirection`/`channelMomentSum` directly (`:1393-1397`). **Every other block: yes** — loops `blockAtomOffset[block]..[block+1]` and reads `atomDirection` directly for every member atom (`:1401-1405`), not through the ghost-aware `effectiveAtomDirection` indirection used elsewhere. Confirmed non-issue: block membership is atomistic/non-atomistic as a *whole block* (`initializeAdaptiveWorkScan`, `gpuAdaptiveRuntime.cpp:1712`, `atomistic = state != coarseState` at atom granularity, driven only by that atom's own block's state — there is no partially-atomistic block), so every atom this loop visits in a non-coarse block is a fine/buffer atom that is already authoritative every step regardless of adaptive commit. | Yes for coarse blocks; N/A otherwise (already fresh) | N/A |
| GPU moment update (`GpuMomentUpdater::update()`, `source/gpu_files/gpuMomentUpdater.cpp:148-177`) | Every step, unconditionally, adaptive-agnostic (`gpuSDSimulation.cpp:308`) | Yes for every atom in the system, adaptive or not — `Copy1`/`Copy2` (`:66-120`) recompute `emomM[atom]=mmom[atom]*emom[atom]` for the full atom range every step, and with `mompar∈{1,2}` (the common/default settings, not `mompar==0`) `Mompar1`/`Mompar2` (`:27-62`) derive the **new** `mmom` for every atom from `emom2[atom].z` — i.e. a real, value-dependent read, not a defensive one. This kernel is not adaptive-aware and is not restricted to `activeAtomList()`. | No — this is generic UppASD infrastructure, identical for adaptive and non-adaptive dynamics, and out of this task's file scope (`source/CoarseGraining/`, adaptive parts of `source/gpu_files/`) | **This is the hard blocker** — see §2 below |
| Measurement (`GpuMeasurement::measure`, `measurement/gpuMeasurement.cpp`) | Called every step, but each observable internally gated (`timeToMeasure()`, `:1021-1044`): `avrg_step`, `cumu_step`, `skyno_step`, `ac_step`, `ene_step` | On a firing step, `calculateEmomMSum()` (`:896-916`) reduces the **entire** `emomM` array over all atoms, unfiltered by adaptive state | No, on a firing step | Yes, but only if materialization happens *before* the firing step's `emomM` is built by the moment updater — see §2 |
| Correlations (`GpuCorrelations::measure`, `correlations/gpuCorrelations.cpp`) | Called every step, internally gated by `sc_step`/`sc_sep` | On a firing step, `GPUSqSum` reads the full `emomM` array | No, on a firing step | Same caveat as measurement |
| CPU-side `emomM` copy (`cpuRestMeasurement`) | Gated by `alwaysCopy \|\| call_fortran_do_measurements(mstep)`; `alwaysCopy` defaults false | Only touches arrays when the Fortran sampling schedule fires | No, on a firing step | Same caveat |
| Restart/checkpoint (`prnrestart`, `System/restart.f90:719-750`; `prn_mag_conf`) | **Not periodic on the GPU/adaptive path** — called only at end of initial phase (`sd_driver.f90:1004-1011`) and end of measurement phase (`uppasd.f90:362-370`), i.e. twice per run, not per step | Yes, writes `emom`/`mmom` for every atom, every ensemble, no block shortcut, no adaptive-CG special case | No | Yes — trivially, since it already happens only twice per run; one explicit `materializeCoarseAtoms(reason=Output)` call before each site is sufficient |
| `GpuSimulation::copyToFortran()` (`gpuSimulation.cpp:1332`) | Only from `printMdStatus` at ~5% cadence (`gpuSDSimulation.cpp:39-52`), and at phase end | Yes, copies full `emom`/`emom2`/`emomM`/`mmom*` | No | Yes, same as restart — materialize immediately before it runs |
| Diagnostics (`diagnosticSnapshot()`, `gpuAdaptiveRuntime.cpp:3486-3572`) | Only `adaptiveDiagnostics>0` (teardown, once) or `adaptiveDiagnostics>=3 && !adaptiveMaskEnabled` (static-parity stage traces; masking disabled means transitions can't occur anyway) | Yes, downloads/reduces the full atom-resolved buffer | No | Yes — diagnostics already know they want a full snapshot; materialize immediately before taking it |
| **CPU adaptive path** — `reconstruct_coarse_atoms` (`source/CoarseGraining/adaptivecgproduction.f90:1227`) | Every step, unconditionally (not gated by `update_interval`) | Yes — overwrites every coarse atom's `emom` from the coarse aggregate every step | Yes, by construction it's a pure write, not a stale-value-dependent read (`validate_reconstruction_inputs` only checks shape/finiteness) | The CPU path has the exact same structural problem as the GPU path — this audit covers it, but no CPU implementation change is proposed in this task; see "Scope decision" |
| CPU `emomM` rebuild (`adaptivecgproduction.f90:1247-1251`) | Every step, for every atom | Yes | No | Same CPU-path caveat |
| CPU `restrict_channel_moments` (`adaptivehybridsolver.f90:170-253`) | Only every `update_interval` steps (`adaptivecgproduction.f90:1324`, gate `mod(step,update_interval)/=0`), plus once more inside `apply_adaptive_transitions` for a proposed coarse transition on that same step | Sums over every dynamic-channel atom, no atomistic/coarse filter | No on a firing step | Already occasional; same "mixed block" caveat as the GPU `restrictAdaptiveMoments` case likely applies but was not verified kernel-by-kernel for the CPU path in this task |
| CPU statichybridoperator field pass (`statichybridoperator.f90:359-384`) | Every step (called from `evaluate_all_ensembles`) | Receives the whole `emom` array but only *reads* `fine_direction(:,atom)` where `atomistic_atom(atom)` is true (`:382-384`) — coarse atoms use `ghost_direction` instead | Yes, coarse atoms' `emom` values are unused here | Already effectively lazy for this consumer |
| MC path (`gpuMCSimulation.cpp`, `gpuMetropolis*.cpp`) | N/A | Not integrated with adaptive CG at all (no references found) | N/A | N/A — out of scope |

## 2. The hard coupling: why "zero full-lattice reconstruction" is not a
   free change

CGP-03's performance gate requires: *"For a static all-coarse state with no
selector update/output/transition, the steady-state target is: zero
full-lattice atom reconstruction."*

The audit shows the field/selector/transition/dipole/output/diagnostics
consumers can all be satisfied by an explicit, occasional
`materializeCoarseAtoms(reason, ...)` call at the right boundary (selector
update, output, transition) — that part of CGP-03 is tractable and low-risk,
and is exactly what Part 3/4/5/6 ask for.

The blocker is **`GpuMomentUpdater::update()`**, which:

1. runs **every step, unconditionally**, for every atom in the system —
   this is pre-existing, adaptive-agnostic UppASD infrastructure, not
   something CGP-03 introduced or can gate;
2. with `mompar` 1 or 2 (real production settings, not merely `mompar==0`),
   **derives that step's `mmom` for every atom from that step's `emom2.z`**
   (`gpuMomentUpdater.cpp:40-42,58-61`) — a genuine value dependency, not a
   defensive/redundant one;
3. unconditionally rebuilds `emomM = mmom * emom` for every atom
   (`Copy1`/`Copy2`, `:81-90,110-118`) immediately after the swap — and
   `emomM` is exactly what measurement's `calculateEmomMSum()` and
   correlation's `GPUSqSum` consume on a sampling-firing step.

Concretely: if `synchronizeAtomicState`'s full commit is skipped on an
ordinary step, `emom2` (and hence, after the swap, `emom`) for an
interior-coarse atom keeps whatever value it held from the last
materialization. `GpuMomentUpdater::update()` still runs its full O(N) sweep
regardless (it cannot be gated by this task — it is shared with
non-adaptive dynamics) and will compute `emomM`/`mmom` for that atom **from
the stale direction**, silently, every step, with no counter or flag
distinguishing "fresh" from "stale" input. Two concrete correctness
consequences:

* **Measurement/correlation on a sampling-firing step reads whatever
  `emomM` the moment updater last built**, which is built from the *prior*
  step's `emom`/`emom2` — so materializing "before measurement" requires
  materializing during the *previous* step's `advanceAdaptiveStep`, before
  that step's `GpuMomentUpdater::update()` call, which in turn requires
  `advanceAdaptiveStep` to know in advance whether the *next* step
  (`mstep+1`) will hit any of `avrg_step`/`cumu_step`/`skyno_step`/
  `ac_step`/`ene_step`/`sc_step`/`sc_sep`/`tottraj_step`/`traj_step(i)`.
  That is a new coupling from the adaptive GPU runtime to the full
  measurement/correlation sampling-interval configuration, not present
  today in any form, and getting the look-ahead condition wrong fails
  silently (wrong energies/averages on some steps, not a crash) rather than
  loudly.
* **`mompar∈{1,2}` longitudinal-fluctuation physics freezes for
  interior-coarse atoms** between materializations: today,
  `commitAdaptiveGhosts` re-copies `coarseDirection.z` (which the coarse
  integrator moves every step) into every coarse atom's `emom2.z` every
  step, so `mmom` for a coarse atom under `mompar 1/2` genuinely tracks the
  coarse block's motion every step. If reconstruction is skipped on
  ordinary steps, that per-step `mmom` tracking stops for coarse atoms
  between materializations — a real (likely small, unquantified) physics
  deviation specific to `mompar 1/2`, not merely a performance change.

An earlier draft of this audit suspected a second, "mixed-block" subtlety
for `restrictAdaptiveMoments` and `addAdaptiveDipole` (atoms in a block
that is only partially atomistic, not covered by the CGP-02 ghost list).
That concern does not apply: block membership is atomistic/non-atomistic
as a whole block (`initializeAdaptiveWorkScan`, `gpuAdaptiveRuntime.cpp:
1709-1727`, `atomistic = state != coarseState` computed once per block and
applied to every atom in it) — there is no block with some atomistic and
some non-atomistic member atoms. Both consumers are therefore already
satisfied by the existing coarse/non-coarse block split; no new compact
list is needed for them.

## 3. Authoritative-state contract (Part 2)

Given the audit above, the following contract is proposed (not yet
implemented):

* **Fine atoms**: `atomDirection` is authoritative and current every step
  (unchanged from today).
* **Interface ghost atoms** (CGP-02 compact list): their prolongated
  direction is valid only for the field evaluation that requested it;
  refreshed every step already (cheap, compact — CGP-02).
* **Interior atoms of fully-coarse blocks**: `coarseDirection` is
  authoritative; `atomDirection`/`emom`/`emom2` may be stale/unmaterialized
  between explicit materialization events, **except** that
  `GpuMomentUpdater` will still silently consume whatever stale value is
  present every step for `emomM`/`mmom` (with `mompar 1/2`) unless
  `GpuMomentUpdater` itself is made adaptive-aware — which this task's
  audit shows is required for full correctness, not merely for
  performance.
* Routines allowed to assume full atom materialization: `GpuMomentUpdater`
  (today, unconditionally — this is exactly the problem), restart/output
  writers (`prnrestart`, `copyToFortran`, `diagnosticSnapshot`), and
  measurement/correlation on a sampling-firing step (transitively, via
  `emomM`).

## 4. Scope decision needed before Part 3-6 implementation

This audit shows CGP-03, taken literally (the "zero full-lattice
reconstruction" performance gate on ordinary steps), is not achievable as a
self-contained change to `GpuAdaptiveRuntime`/`gpuSimulation.cpp` alone: it
requires either

**(a)** also making `GpuMomentUpdater` adaptive-aware (scope it to
`activeAtomList()` for the direction-dependent parts, and define what
`mmom`/`emomM` mean for a non-materialized coarse atom), **and** adding a
look-ahead into the measurement/correlation sampling schedule so
materialization happens one step before any sampling-firing step — both
materially larger and riskier than CGP-03's own text implies, verging into
CGP-04/CGP-05 territory (which explicitly plan to touch field/clearing
scope, not moment update); or

**(b)** a bounded implementation that: introduces the
`materializeCoarseAtoms(reason, ...)` helper Part 3 asks for; uses it to
replace the *ad hoc* two call sites with explicit, named reasons
(`OrdinaryStep`, `SelectorUpdate`, `Transition`, `Output`, `Diagnostics`);
keeps the **ordinary-step full commit exactly as expensive as it is today**
(because of the `GpuMomentUpdater`/measurement coupling above), but
eliminates the genuinely redundant second full commit at
`gpuSimulation.cpp:1266` when a selector-update step runs but **no
transition was actually published** (a real, provable no-op today); and
documents the `GpuMomentUpdater` coupling as the reason the headline "zero
full-lattice reconstruction" target is not reached, explicitly flagged for
Terra, matching the CGP-00B precedent (Part G deferral) already established
in this project.

(a) matches the letter of the CGP-03 checklist but is a substantially
larger, higher-risk change than a routine task slice — it touches shared
non-adaptive infrastructure (`GpuMomentUpdater`) and introduces a new
cross-cutting dependency (sampling-schedule look-ahead) with a
silent-failure mode if the look-ahead condition is ever wrong. (b) is safe
and still delivers a real, measurable win (removes one of two full-lattice
passes on selector-update steps, consolidates scattered reconstruction
calls behind one named helper per Part 3's own requirement) but will not
satisfy the letter of the "zero full-lattice reconstruction" checklist
item, and is recorded honestly as a partial result rather than checked off.

**Human decision: (b), the bounded implementation.** Implemented below.

---

## 5. Part 3 — `materializeCoarseAtoms(reason, ...)` (implemented)

`GpuAdaptiveRuntime::materializeCoarseAtoms(GpuAdaptiveMaterializationReason
reason, real* atomDirection, const GpuAdaptiveReconstructionPolicy& policy)`
(`gpuAdaptiveRuntime.hpp`/`.cpp`) is a thin wrapper: it calls the unchanged
`synchronizeAtomicState()` (no reconstruction formula, kernel, or launch
geometry changed by this task) and increments a per-reason counter,
`materializationCount(reason)`, reachable from tests. Two reasons exist,
matching the two real production call sites
(`GpuAdaptiveMaterializationReason::OrdinaryStep`, `::Transition`);
`Output`/`Diagnostics` reasons from the original design sketch were dropped
rather than added as dead code, because the bounded scope keeps
`OrdinaryStep` unconditional, which already leaves `copyToFortran()`,
`prnrestart`/`prn_mag_conf`, and `diagnosticSnapshot()` (§1's audit) exactly
as fresh as before this task — nothing new needed to call
`materializeCoarseAtoms` on their behalf.

`gpuSimulation.cpp:advanceAdaptiveStep` now calls
`materializeCoarseAtoms(OrdinaryStep, ...)` where it used to call
`synchronizeAtomicState(...)` directly (unconditional, every step, cost
unchanged — see §2 for why), and
`materializeCoarseAtoms(Transition, ...)` where it used to call
`synchronizeAtomicState(...)` unconditionally inside the
`adaptiveUpdateInterval` block — now gated (§7).

## 6. Part 4 — selector behaviour (resolved, no code change required)

`selectorAdaptiveScores`/`restrictAdaptiveMoments` (both only run inside the
`adaptiveUpdateInterval` block, after the `OrdinaryStep` commit has already
run earlier in the same step) need current `atomDirection` for whatever
atoms their edges/channels touch. Because the `OrdinaryStep` commit stays
unconditional in the bounded scope, this need is already satisfied by
construction — the "least invasive option" Part 4 asks for turns out to be
"no change", which is recorded explicitly (with a code comment at the call
site, `gpuSimulation.cpp`) rather than left implicit.

## 7. Part 6 — transition materialization (implemented)

A new one-element device counter, `GpuAdaptiveDeviceRuntime::
transitionsPublished`, is `atomicAdd`'d by `publishAdaptiveState`
(`gpuAdaptiveRuntime.cpp`) exactly where `blockState[block]` is actually
overwritten (`oldState != newState`, after the rollback/reject early
returns) — so a rejected (`acceptedMask[block]==0`) or numerically rolled
back (`reconstructionOk==false`) proposal is correctly never counted, only a
genuine, committed state change is. It is cleared on-stream immediately
before each `publishAdaptiveState` launch (`GPU_MEMSET_ASYNC`, no host wait)
and read back inside `publishProposedState()`'s existing
`refreshHostWorkCounts()` stream sync (piggy-backed, no new host wait
introduced) into `transitionsPublishedLastCall()`.

`gpuSimulation.cpp` now calls `materializeCoarseAtoms(Transition, ...)`
**only when `transitionsPublishedLastCall() > 0`**. This is safe in both
transition directions, for different reasons (documented in a code comment
at the call site and verified by
`testTransitionMaterializationNegativeControl`, §8):

* **coarse → fine**: `publishAdaptiveState` already writes the block's
  atoms directly (the `oldState == coarseState` reconstruction branch), so
  the follow-up commit is a harmless no-op for those atoms —
  `commitAdaptiveGhosts` skips them outright once compaction's
  `atomisticAtomMask` update (run inside `publishProposedState`, before
  this decision) takes effect.
* **fine/buffer → coarse**: `publishAdaptiveState` reclassifies the block
  (`blockState[block] = newState`) but — this was not obvious before
  reading the kernel closely — **never touches that block's `atomDirection`
  in this direction**. The atoms are left at their stale, individually
  varying pre-transition values until a commit aligns them to
  `coarseDirection`. Skipping the follow-up commit here would be a genuine
  correctness bug, not merely a missed optimization.

When `transitionsPublishedLastCall() == 0` (no block changed state — the
common case once an adaptive run reaches a stable partition), skipping the
second commit is exactly correct: the `OrdinaryStep` commit earlier in the
same step already left every coarse atom aligned to its (unchanged)
`coarseDirection`, so a second commit would recopy identical values.

## 8. Tests (Part 6 negative control, `tests/coarse_graining/test_gpu_adaptive_runtime.cpp`)

* `testTransitionCounterAndMaterialization()`: rejecting every proposed
  transition publishes exactly zero (`transitionsPublishedLastCall()==0`)
  and leaves `blockState` untouched; accepting every proposal publishes
  exactly as many transitions as an independently-computed before/after
  `blockState` diff (not a hard-coded expected count, so the test does not
  silently drift if the fixture's selector-score arithmetic ever changes);
  `materializeCoarseAtoms()` reproduces `synchronizeAtomicState()` bitwise
  and credits only the reason it was called with.
* `testTransitionMaterializationNegativeControl()`: isolates a single
  fine→coarse transition (accepted-mask selects only that block), then
  proves, by direct download of `atomDirection`, that (1) immediately after
  `publishProposedState` the transitioned block's atoms are still at their
  pre-transition value (proving the "coarse→fine writes directly,
  fine→coarse does not" claim above, not merely asserting it), and (2) after
  `materializeCoarseAtoms(Transition, ...)` those atoms are exactly aligned
  to `coarseDirection`. This is CGP-03's own required negative control
  ("skip reconstruction before a transition and confirm parity fails")
  turned into a passing/failing assertion pair rather than a manual
  demonstration.

Both new tests pass; the full existing `gpu_adaptive_runtime_tests` binary
(53 fixtures) passes unchanged.

## 9. Correctness suite

Fresh out-of-tree build: `build_cgp03_cuda_fp32`
(`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`,
`CMAKE_BUILD_TYPE=Release`), CUDA 13.3.73, driver 610.57.04, 2x NVIDIA RTX
A4000 (both idle, 0% utilization, <110 MiB used immediately before every
run below).

`ctest -L cg13-cuda` (33 tests, includes `coarse-graining-gpu-adaptive-runtime`,
`coarse-graining-adaptive-asd-parity`, all `adaptive-cg-moving-*` fixtures —
which repeatedly exercise real fine↔coarse transitions and, for
`adaptive-cg-moving-static-mixed`, a held-static interface that produces
many "selector update, zero transitions published" steps, i.e. the exact
skip path this task adds): **32/33 pass**. The one failure,
`adaptive-cg-production-e2e` (`assert abs(float_metric(fft.stdout,
"coarse_dipole")) > 0.0`), reproduces byte-for-byte identically on an
unmodified CGP-02 build (`git worktree` pinned to `2ea4277e`, same CMake
configuration) — confirmed pre-existing and unrelated to this task, the
same pattern CGP-01's evidence doc recorded for its own cg13-cuda run.

All-fine T=0 and finite-T ASD parity ("CG-10.5 finite-temperature adaptive
fine-region execution", "CG-10.5 finite-temperature mixed-resolution
execution") pass within `adaptive-cg-production-e2e`'s own output, before
that test's unrelated later failure.

## 10. Performance evidence

The bounded scope's only performance change is: on a selector-update step
that publishes zero transitions, one fewer full-lattice
`prolongateAdaptiveGhosts(Compact)` + `commitAdaptiveGhosts` launch pair
(previously always 2 full commits per selector-update step regardless of
outcome; now 1 when nothing transitions, 2 when something does). The
`OrdinaryStep` commit's own per-step cost is unchanged (§2) — this task does
not claim a whole-step speedup, only removal of one specific, provably
redundant call.

Measured directly (new `--materialization-overhead` mode,
`benchmark_gpu_adaptive_runtime`, all-coarse 65536-atom fixture — the size
CGP-07 uses as its historical reference point — 3 discarded warm-up calls,
20 timed `materializeCoarseAtoms(OrdinaryStep, ...)` calls per run, device
sync per sample):

| Run | median (us) | MAD (us) | notes |
|---|---|---|---|
| 1 | 12.798 | 0.149 | 3 samples (16,17,19 of 20) show 19–24us, consistent with the shared-host contention/cold-start notes in project memory even though `nvidia-smi` read 0% utilization on both GPUs immediately before the run |
| 2 | 12.720 | 0.120 | rerun immediately after; no outliers |

So each redundant commit removed on a no-transition selector-update step
saves ≈12.7-12.8us at 65536 atoms (both GPUs otherwise idle). This is a
per-call, not a whole-run, number: the aggregate benefit of this task scales
with how many selector-update steps in a real run publish zero transitions,
which is workload-dependent (frequent for a converged/steady partition,
rare during an actively-moving interface) and was not separately swept here
— reporting a whole-run speedup from this alone would overstate what a
single-call microbenchmark can support.

## 11. Outcome summary

Implemented (bounded scope, Parts 3, 4 (no-op), 6):

* `GpuAdaptiveRuntime::materializeCoarseAtoms(reason, ...)` — named,
  countable entry point replacing two `synchronizeAtomicState()` call sites
  hand-picked ad hoc in `gpuSimulation.cpp`.
* Exact, atomic, host-visible count of blocks whose state actually changed
  per `publishProposedState()` call
  (`transitionsPublishedLastCall()`), correctly excluding
  rejected/rolled-back proposals.
* The post-selector-update full commit is skipped when that count is zero —
  the one provable, safe reduction in full-lattice reconstruction work this
  task's audit supports.
* Two new unit tests, including the Part 6 negative control demonstrating
  that skipping materialization after a *genuine* fine→coarse transition
  produces stale, misaligned atom directions (i.e., that the guard is
  correctness-critical, not merely an optimization).

Explicitly **not** implemented (flagged for Terra, matching CGP-00B's Part G
precedent):

* Zero full-lattice reconstruction on an *ordinary* (non-selector-update)
  ADP step. Blocked on making `GpuMomentUpdater` adaptive-aware (its
  `mompar 1/2` paths read `emom2.z` per atom every step, unconditionally,
  for every atom in the system) and adding a measurement/correlation
  sampling-schedule look-ahead so materialization happens one step before
  any `avrg_step`/`cumu_step`/`skyno_step`/`ac_step`/`ene_step`/`sc_step`/
  `sc_sep`/`tottraj_step`/`traj_step` firing. See §2 for the full audit of
  why this is required for correctness, not merely desirable for
  performance, and why it is materially larger than a routine task slice.
* The CPU adaptive path (`source/CoarseGraining/adaptivecgproduction.f90`
  `reconstruct_coarse_atoms`, called every step unconditionally) has the
  same structural issue; audited in §1 but not touched by this task, which
  scoped its implementation to the GPU runtime only.

Performance improved (measurably, on the one call this task removes, under
the specific condition it removes it); neutral everywhere else (the
`OrdinaryStep` commit, the large majority of the removed-call's context, is
byte-for-byte the same cost as before, verified bitwise by
`testTransitionCounterAndMaterialization`); no regression (32/33 cg13-cuda,
matching the pre-existing baseline exactly, with the one failure confirmed
unrelated on the CGP-02 baseline).

## Checklist

* [x] Every coarse-atom direction consumer inventoried.
* [x] Authoritative-state contract documented.
* [x] Fine atom state remains authoritative atomically.
* [x] Coarse block state is authoritative for interior coarse regions.
* [ ] Hot timestep no longer reconstructs every coarse atom. (Only the
  redundant *second* commit on a no-op selector-update step is removed; the
  once-per-step `OrdinaryStep` commit is unchanged -- see §2/§11 for why
  removing it needs `GpuMomentUpdater` changes out of this task's bounded
  scope, flagged for Terra.)
* [x] Interface ghosts remain correct.
* [x] Selector never consumes stale atomic state.
* [x] Transition reconstruction remains correct.
* [x] Atom-resolved output explicitly materializes when necessary. (Already
  true before this task -- the unconditional `OrdinaryStep` commit means
  output/restart/diagnostics never see stale coarse-atom state; no new code
  was needed to preserve this under the bounded scope.)
* [x] Restart/checkpoint output remains correct. (Same reason: untouched by
  this task, and still fed by the unconditional `OrdinaryStep` commit.)
* [ ] Static all-coarse steady state performs no unnecessary full
  reconstruction. (Partially: a static all-coarse *selector-update* step now
  performs zero extra reconstruction beyond the one `OrdinaryStep` commit
  every step already pays; that per-step floor itself remains -- see the
  "Hot timestep" item above.)
* [x] Moving adaptive fixture passes (`adaptive-cg-moving-*`, `ctest -L
  cg13-cuda`).
* [x] Finite-T mixed fixture passes (`adaptive-cg-production-e2e`'s
  finite-temperature mixed-resolution stage, before that test's unrelated
  later failure).
* [ ] Output negative control passes. (Not attempted: the bounded scope
  makes no change to output/restart materialization timing, so there is
  nothing new to demonstrate a negative control against; see the "Atom-
  resolved output" item above.)
* [x] Transition negative control passes
  (`testTransitionMaterializationNegativeControl`).
* [x] Before/after integration/interface timing recorded (§10; per-call
  microbenchmark, not a whole-run figure -- see §10 for why).
* [ ] Terra review performed before merge. (Requested; not yet performed.)

See sections 0-4 above for the full Part 1/2 audit and contract, and
sections 5-11 for the Part 3/4/6 design, tests, and evidence.
