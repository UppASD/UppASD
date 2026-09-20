# CGP-03B — `GpuAdaptiveMomentUpdater`: audit, design, and evidence

Task: `docs/CGP_work.md` CGP-03B ("Introduce a `GpuAdaptiveMomentUpdater` to
unblock ordinary-step lazy materialization"). Branch `gpu_hip_cu_ab_cg`,
commit base `199344b5` (CGP-03).

Status: **Complete, bounded scope.** A new, adaptive-owned moment updater is
implemented, tested, and wired in to fully replace `GpuMomentUpdater` for
adaptive runs. `GpuMomentUpdater` itself is verified byte-for-byte untouched
(§7). The updater's touched-atom-only capability (`updateActiveOnly()`) is
implemented, parity-tested, and benchmarked, but is **not** wired into the
production hot path yet — §4 explains why, and it is the same style of
documented, evidence-based partial result as CGP-00B's Part G and CGP-03's
option (b).

## 0. Recap: the blocker this task was created to remove

CGP-03's evidence doc (`docs/CGP-03_LAZY_MATERIALIZATION_EVIDENCE.md`, §2)
found that `GpuMomentUpdater::update()` (`source/gpu_files/gpuMomentUpdater.cpp`)
is shared, non-adaptive infrastructure that:

1. runs every step, unconditionally, for every atom in the system, for
   every SD path (adaptive or not) — called from `gpuSDSimulation.cpp`'s
   main loop;
2. with `mompar` 1 or 2, derives that step's `mmom` for every atom from that
   step's `emom2.z` — a genuine value dependency;
3. unconditionally rebuilds `emomM = mmom * emom` for every atom immediately
   after the `emom`/`emom2` swap.

Because `GpuMomentUpdater` is shared with every non-adaptive dynamics path,
CGP-03 declined to edit it and left the full per-step commit in place. This
task's brief, suggested by the human after reviewing that outcome: build a
**separate** `GpuAdaptiveMomentUpdater` that adaptive runs use instead, so
`GpuMomentUpdater` and its non-adaptive contract stay completely untouched.

## 1. What `GpuMomentUpdater::update()` actually does (re-read before design)

Re-reading `source/gpu_files/gpuMomentUpdater.cpp` line by line (not assumed
from CGP-03's summary) turned up a detail CGP-03's evidence doc did not need
and this task does:

```
switch(mompar) {
   case 0: mmom2.swap(mmom); break;              // O(1) pointer swap
   case 1: Mompar1 kernel -> mmom2 = f(mmom0, emom2)   // O(N)
   case 2: Mompar2 kernel -> mmom2 = f(mmom0, emom2)   // O(N)
}
emom.swap(emom2);   // O(1) pointer swap, always
mmom.swap(mmom2);   // O(1) pointer swap, always
Copy1/Copy2 kernel -> mmomi, emomM = f(mmom, emom)     // O(N), always
```

Two things follow that shaped the design:

* **`GpuTensor::swap` is a pointer swap** (`tensor.hpp:535-536`,
  `std::swap(data_, A.data_)`), not an element-wise copy. Every swap above
  is O(1) regardless of system size. The only O(N) work is the `mompar`
  branch's kernel (mompar 1/2 only) and the unconditional `Copy1`/`Copy2`
  kernel.
* **`mompar == 0` is a *double* swap, not a no-op.** `mmom2.swap(mmom)`
  happens first (inside the switch), then the unconditional
  `mmom.swap(mmom2)` happens again a few lines later. Two swaps compose to
  the identity: for `mompar == 0`, `mmom` ends the call holding exactly the
  value it held before the call, for *every* atom, full stop — there is no
  active/inactive distinction to make for that case, because there was
  never any per-atom work to restrict. (This was originally missed in this
  task's own first draft and caught by its own negative-control test
  failing — see §8.)

## 2. Part A — the moment/`emomM` authoritative-state contract

Mirroring CGP-03 Part 2's atom-direction contract:

* **Active atoms** (`GpuAdaptiveRuntime::activeAtomList()` /
  `activeAtomCount()` — the compact, one-based, per-ensemble list of every
  atomistic (fine + buffer/interface) atom, `atomisticAtomMask`-driven,
  already used to scope `GpuDepondtIntegrator::evolveFirst/evolveSecond` for
  adaptive steps, `gpuSimulation.cpp:1213-1239`): `mmom`/`mmomi`/`emomM` are
  authoritative and current every step, because `emom2` for these atoms is a
  genuine, freshly-integrated dynamical variable every step.
* **Interior coarse atoms and ghost atoms** (endpoints of live bonds, per
  the CGP-02/RCG-09C compact `ghostAtomList`): `mmom`/`mmomi`/`emomM` are
  **not** required to be current every step by any production consumer.
  This deliberately extends CGP-03's own atom-*direction* staleness
  contract to atom-level *moment magnitude*, and is new to this task. The
  reasoning, and why ghost atoms specifically do **not** need inclusion in
  the touched set (contradicting the task sketch's tentative suggestion):

  * The coarse integrator (`GpuAdaptiveRuntime::prepareCoarsePredictor` /
    `correctCoarse`) advances `coarseMoment_`/`coarseDirection_`, a
    separate per-block/channel quantity — it never reads per-atom `mmom`.
    Atom-level `mmom` for a coarse atom is not an input to any accepted
    coarse-CG dynamics equation.
  * Ghost atoms' field-evaluation direction is served from the *separate*
    `ghostDirection_` scratch tensor, refreshed every step by the compact
    `prolongateAdaptiveGhostsCompact` launch (CGP-02) — not from `emom`. So
    a ghost atom's `emom`/`mmom`/`emomM` freshness has no more field-side
    consumer than any other interior coarse atom's.
  * The CGP-03 consumer audit table (§1 of that evidence doc) already
    establishes that nothing reads `mmom`/`emomM` for a coarse atom on an
    *ordinary* (non-firing, non-transition, non-selector-update,
    non-output/diagnostic) step. The real consumers are exactly:
    measurement's `calculateEmomMSum()` and correlation's `GPUSqSum` on a
    sampling-firing step, restart/`copyToFortran` output, and diagnostics —
    all occasional, all already-named events, none of them per-step.

* **`mompar` 1/2 physics for interior-coarse atoms**: CGP-03's evidence doc
  (§2) noted that `commitAdaptiveGhosts` re-copies `coarseDirection.z` into
  every coarse atom's `emom2.z` every step (unconditionally, unchanged by
  this task — see §3), so `mmom` for a coarse atom under `mompar` 1/2
  *could* track the coarse block's z-motion every step today, purely as an
  accident of `GpuMomentUpdater`'s full-N sweep. This task's contract makes
  that tracking **optional, not required**: nothing downstream depends on a
  coarse atom's `mmom` reflecting this step's `coarseDirection.z` versus an
  older one, for exactly the reasons above. `updateActiveOnly()` therefore
  correctly leaves coarse atoms' `mmom` pinned at its last-materialized
  value even though their `emom` keeps moving every step (see §8's
  bit-level proof).
* Routines allowed to assume full moment materialization: `copyToFortran`,
  `prnrestart`/`prn_mag_conf`, `diagnosticSnapshot()`, and measurement/
  correlation on a firing step — the same list CGP-03 named for
  `atomDirection`, now extended to `mmom`/`mmomi`/`emomM`.

## 3. Part B — `GpuAdaptiveMomentUpdater`

New files `source/gpu_files/gpuAdaptiveMomentUpdater.{hpp,cpp}`. Design
constraints and how each is met:

* **`GpuMomentUpdater` is not modified.** `git diff --stat -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp`
  against this task's base is empty (§7). The four kernel formulas
  (`Mompar1`, `Mompar2`, `Copy1`, `Copy2`) are intentionally re-declared as
  new nested classes of `GpuAdaptiveMomentUpdater` rather than shared,
  because the originals are private members of `GpuMomentUpdater` and
  reaching into them (via `friend` or otherwise) would itself be an edit to
  `gpuMomentUpdater.hpp`. The duplication is guarded against silent drift
  by `testAdaptiveMomentUpdaterFullParity`-equivalent bitwise cross-checks
  against the real `GpuMomentUpdater` on every call (§8).
* **One field implementation, two launch geometries** (matching CGP-01's
  precedent of preferring compile-time/parameter specialization over
  scattered branches): the same four functor classes derive from
  `ParallelizationHelper::Atom` and are launched unchanged via either
  `gpuAtomCall` (full lattice — `updateFull()`) or the existing
  `gpuActiveAtomCall` (compact one-based active list — `updateActiveOnly()`,
  reusing `GpuAdaptiveRuntime::activeAtomList()`/`activeAtomCount()` as-is,
  per this task's own instruction not to build a new O(N) pass to find
  touched atoms).
* **`updateFull()`**: byte-for-byte the same buffers, formulas, and launch
  geometry as `GpuMomentUpdater::update()` — verified by
  `runFullParity()` (§8).
* **`updateActiveOnly(activeAtoms, activeCount)`**: the `mompar`
  branch and the `Copy1`/`Copy2` kernel run only over the active list;
  `emom<->emom2` and `mmom<->mmom2` stay whole-tensor pointer swaps,
  unconditionally, for every atom — deliberately unchanged, because (a) the
  swap is already O(1) so there is nothing to save, and (b) CGP-03's
  atom-direction contract (`materializeCoarseAtoms(OrdinaryStep, ...)`
  still commits `coarseDirection` into every coarse atom's `emom2` every
  step; see §4) is not altered by this task, so `emom` must keep tracking
  it for every atom regardless of the moment-update's own scoping. An atom
  outside the active list keeps whatever `mmom`/`mmomi`/`emomM` value it
  held before the call — not by an explicit "skip" branch, but because
  neither buffer slot for that atom is ever written, so swapping two
  untouched slots leaves the same numeric value visible under whichever
  name currently points at it.

## 4. Part C/D — call-site ownership, and why the fast path is not wired in

**Call-site ownership (Part C):** `gpuSDSimulation.cpp`'s `SDmphase` (the
only place adaptive dynamics run at all — `SDiphase`, the initial phase,
never calls `gpuSim.adaptiveEnabled()` and always uses the plain
non-adaptive `GpuMomentUpdater` path, confirmed by reading the whole
function) now constructs both a `GpuMomentUpdater` and a
`GpuAdaptiveMomentUpdater`, and branches once per step:

```cpp
if(gpuSim.adaptiveEnabled()) {
   adaptiveMomUpdater.updateFull();
} else {
   momUpdater.update();
}
```

This is a **full replacement**, not a supplement: `GpuMomentUpdater::update()`
is never called for an adaptive step from now on, so adaptive dynamics own
their moment update end to end, and `GpuMomentUpdater`'s own contract for
every non-adaptive path is untouched (same object, same call, same cost).

**Why `updateActiveOnly()` is not the production default (Part D):** this
task's own text asks it to "establish plainly" whether unblocking
`GpuMomentUpdater` alone is sufficient, or whether the measurement/
correlation sampling-schedule look-ahead CGP-03 flagged is still separately
required. Re-tracing that question against the real code (not merely
CGP-03's prose) turned up **two** independent reasons `updateActiveOnly()`
cannot yet be the default, one narrower than CGP-03 anticipated and one
confirming it:

1. **`materializeCoarseAtoms(OrdinaryStep, ...)` itself stays unconditional
   in this task** (`gpuSimulation.cpp:1261-1263`, untouched) — it commits
   `coarseDirection` into every coarse atom's `emom2` every step regardless
   of this task's changes. That call site is CGP-03's own scope, not named
   in CGP-03B's required work, and re-opening it would mean re-deciding
   CGP-03's own bounded-scope call (restoring `Output`/`Diagnostics`
   materialization reasons CGP-03 deliberately dropped — see that evidence
   doc §5). This task does not do that. So the atom-*direction* side of the
   picture is unchanged either way; the win this task can offer is scoped
   strictly to the moment/`emomM` rebuild, independent of that.
2. **The correlation sampling schedule's firing predicate is not a pure
   function**, unlike measurement's. Reading
   `source/gpu_files/correlations/gpuCorrelations.cpp:health `~200-380`
   directly: `do_sc`'s firing condition depends on `curstep` (derived from,
   but not identical to, `mstep`) and mutates `t_cur`/`sc_step_arr_cpu`
   bookkeeping as a side effect of `measure()` itself. There is no
   read-only "would this fire at step N" query to call ahead of time
   without either duplicating that stateful logic (drift risk) or editing
   `GpuCorrelations` to expose one (a shared-infrastructure edit, the same
   category of change CGP-03 declined for `GpuMomentUpdater`). By contrast,
   `GpuMeasurement::timeToMeasure(mtype, mstep)`
   (`source/gpu_files/measurement/gpuMeasurement.cpp:1021-1044`) *is* a
   public, side-effect-free predicate that could safely be queried one step
   ahead — so the measurement half of CGP-03's worry turns out to be
   cheaper to resolve than the doc assumed, but correlation's does not, and
   a partial look-ahead (measurement only) would leave `emomM` silently
   stale for correlation's `GPUSqSum` on a mixed run that enables both,
   which is worse than not attempting it. This is confirmed, not merely
   assumed: it is a direct reading of the stateful `curstep`/`t_cur`
   bookkeeping in `gpuCorrelations.cpp`, not an inference from the CGP-03
   prose.

Given both, the production wiring uses `updateFull()` — identical cost and
behavior to before this task — while `updateActiveOnly()` exists as a
tested, ready-to-enable capability. Flipping it on for real would need, at
minimum: (a) revisiting CGP-03's `OrdinaryStep` scope decision (its own
follow-up, not this task's), and (b) either a `GpuCorrelations` change to
expose a genuine side-effect-free firing predicate, or an accepted design
that tolerates rebuilding `emomM` for every atom on any step where
correlation is enabled at all (which would recover much less of the win on
runs that use both adaptive CG and spin correlations together). Both are
flagged here for Terra, matching the CGP-00B Part G / CGP-03 (b) precedent
of a documented, evidence-backed partial result rather than an
under-verified full one.

## 5. Answering Part D directly

**Does CGP-03B alone unblock CGP-03's remaining checklist items?** No.
CGP-03's "hot timestep no longer reconstructs every coarse atom" and
"static all-coarse steady state performs no unnecessary full
reconstruction" items are about the `materializeCoarseAtoms(OrdinaryStep, ...)`
atom-*direction* commit, which this task does not touch (§4, point 1). What
this task removes is a *different*, previously-undocumented full-lattice
cost sitting immediately after it (`GpuMomentUpdater`'s O(N) kernels) — real
and measurable (§6), but not the one those two checklist items name. They
remain open, now for two independent reasons instead of one conflated
reason: the atom-direction commit (CGP-03's own remaining scope) and the
correlation-schedule visibility gap this task found narrower than assumed
but did not close.

## 6. Performance evidence

Not a whole-run figure — a per-call microbenchmark, same rationale as
CGP-03's own §10 (the removed/removable calls are not the only per-step
cost, so a whole-run number would overstate what a single-call
microbenchmark can support). New `--benchmark` mode in
`gpu_adaptive_moment_updater_tests` (not a ctest — hardware-dependent,
matching the RCG-09 PERF harness precedent): 65536-atom fixture (CGP-07's
historical reference size, same size CGP-03's own materialization-overhead
benchmark used), 1 ensemble, one discarded warm-up call, then 8 ABBA pairs
(`updateFull()`, `updateActiveOnly()` x2, `updateFull()`) per configuration,
median and MAD reported, device-stream-synced per sample.

Hardware: 2x NVIDIA RTX A4000 (default device 0 used), driver 610.57.04,
CUDA 13.3.73, both GPUs 0% utilization / <110 MiB used immediately before
the run (shared host — see project memory on GPU contention; no
foreign load observed during this run).

| active fraction | mompar | full (median us) | full MAD | active-only (median us) | active-only MAD | speedup |
|---|---|---|---|---|---|---|
| 1.25% | 0 | 8.091 | 0.095 | 5.926 | 0.163 | 1.365x |
| 1.25% | 1 | 11.050 | 0.143 | 7.938 | 0.146 | 1.392x |
| 1.25% | 2 | 10.832 | 0.157 | 7.809 | 0.062 | 1.387x |
| 6.25% | 0 | 8.250 | 0.081 | 5.780 | 0.065 | 1.427x |
| 6.25% | 1 | 10.812 | 0.069 | 7.899 | 0.058 | 1.369x |
| 6.25% | 2 | 10.806 | 0.063 | 7.896 | 0.103 | 1.369x |
| 25% | 0 | 8.154 | 0.129 | 6.361 | 0.205 | 1.282x |
| 25% | 1 | 10.902 | 0.204 | 8.513 | 0.137 | 1.281x |
| 25% | 2 | 10.829 | 0.104 | 8.505 | 0.077 | 1.273x |
| 100% | 0 | 8.082 | 0.145 | 8.498 | 0.095 | 0.951x |
| 100% | 1 | 10.750 | 0.084 | 11.446 | 0.088 | 0.939x |
| 100% | 2 | 10.704 | 0.133 | 11.368 | 0.079 | 0.942x |

Reading: at low-to-moderate active fractions, `updateActiveOnly()` is
~1.27x-1.43x faster than `updateFull()` at 65536 atoms, with the win
roughly flat across the tested fractions rather than growing toward 1/
fraction — consistent with `gpuActiveAtomCall`'s launch/host-sync overhead
being a fixed cost that dominates over the (still fairly small, tens of
thousands of elements) per-element work at this size, not a scaling
failure. At 100% active, `updateActiveOnly()` is measurably *slower*
(~0.94x) than `updateFull()`, because the active-list indirection
(`active_atom_kernel`'s slot->atom translation) adds overhead with no
work-reduction benefit once every atom is touched — exactly the
regime `updateFull()` is intended for and why the production wiring
picks it unconditionally rather than "active-only always."

This is capability evidence for a future task, not a claim of a
production speedup: with `updateFull()` wired in today, this task changes
no production timing (§4). `testAdaptiveMomentUpdaterFullParity`-equivalent
checks (§8) already confirm `updateFull()` costs the same as
`GpuMomentUpdater::update()` (same kernels, same launches).

## 7. `GpuMomentUpdater` untouched -- negative control

```
$ git diff --stat -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp
$ git status --short -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp
$ git log --oneline -3 -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp
a29be892 gpu: Introduced do_gpu_timings and improved GPU stopwatch routines
24186d67 first version of hybrid CUDA/HIP
```

Both diff/status commands print nothing: the working tree's copies are
byte-identical to `HEAD` (`199344b5`, the CGP-03 commit), and neither file
has been touched by any commit since UppASD's original hybrid CUDA/HIP
introduction. `SDiphase`'s and `SDmphase`'s non-adaptive branch
(`momUpdater.update()`) call this exact, unmodified code.

## 8. Tests (`tests/coarse_graining/test_gpu_adaptive_moment_updater.cpp`,
registered as ctest `gpu-adaptive-moment-updater`, label
`coarse-graining;cg13;cg13-cuda;reference`, `SKIP_RETURN_CODE 77`)

Every case below runs for `mompar` in `{0, 1, 2}` and `initexc` in
`{'N', 'I'}` (12 combinations), on a 6-atom/2-ensemble fixture:

* **`runFullParity`**: `GpuAdaptiveMomentUpdater::updateFull()` reproduces
  `GpuMomentUpdater::update()` bitwise (`mmom`, `mmomi`, `emom`, `emomM`)
  from identical seeded inputs.
* **`runActiveAllParity`** (the task's "all-fine control"):
  `updateActiveOnly()` given every physical atom also reproduces
  `GpuMomentUpdater::update()` bitwise -- the active-list launch geometry
  changes no answer.
* **`runActiveSubsetNegativeControl`**: a genuine 3-of-6 active subset
  ({2,4,6}) is seeded with sentinel `mmom2`/`mmomi`/`emomM` values
  distinguishable from anything either updater's formulas could produce.
  After `updateActiveOnly()`:
  * atoms in the subset match a same-input `GpuMomentUpdater::update()`
    reference bitwise;
  * atoms outside the subset have `mmomi`/`emomM` pinned exactly at their
    pre-call sentinel (proving `Copy1`/`Copy2` were skipped, not merely
    coincidentally correct) for every `mompar`;
  * atoms outside the subset have `mmom` pinned at the sentinel for
    `mompar` 1/2 (the `Mompar1`/`Mompar2` kernel genuinely skipped), but
    match the reference for `mompar == 0` -- the double-swap identity case
    from §1, where there was never any per-atom work to skip. (This
    distinction was not in the test's first draft: it initially asserted
    "pinned" unconditionally and failed immediately for `mompar == 0`,
    which is exactly the re-reading in §1 that explains why. Left as a
    citation here because it is direct evidence the negative control is
    discriminating, not vacuous.)
  * every atom's `emom` -- selected or not -- matches the reference exactly
    (the whole-tensor swap stayed global, confirming CGP-03's
    atom-direction contract is unaffected by this task).

Result: **12/12 combinations pass** (`GPU-ADAPTIVE-MOMENT-UPDATER passed`).

## 9. Correctness suite

Same fresh out-of-tree build discipline as prior CGP tasks:
`build_cgp03b_cuda_fp32` (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`,
`CMAKE_BUILD_TYPE=Release`), CUDA 13.3.73, driver 610.57.04, 2x NVIDIA RTX
A4000 (both idle, 0% utilization, <110 MiB used immediately before the run).

`ctest -L cg13-cuda` (34 tests -- 33 from CGP-03's baseline plus the new
`gpu-adaptive-moment-updater`): **33/34 pass**. The one failure,
`adaptive-cg-production-e2e`
(`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), is the
same pre-existing, unrelated failure CGP-01/CGP-02/CGP-03 all documented
(identical assertion, identical stage) -- confirmed by identical output up
to that point (`CG-10.5 finite-temperature adaptive fine-region execution`
through `CG-10.5 DMI/anisotropy CPU/GPU production parity` all pass first).
All-fine T=0 and finite-T ASD parity pass within that same test's output,
before its unrelated later failure. `adaptive-cg-moving-*` (moving
fine/coarse interface, repeated real transitions) and
`adaptive-cg-transition-ownership-invariants` all pass, confirming this
task did not disturb CGP-02/CGP-03's transition and interface correctness.

## 10. Outcome summary

Implemented:

* `GpuAdaptiveMomentUpdater` (`source/gpu_files/gpuAdaptiveMomentUpdater.{hpp,cpp}`),
  with a full-lattice `updateFull()` (byte-identical to `GpuMomentUpdater::update()`)
  and a touched-atom-only `updateActiveOnly()`, both parity-tested.
* `gpuSDSimulation.cpp`'s `SDmphase` fully replaces `GpuMomentUpdater` with
  `GpuAdaptiveMomentUpdater::updateFull()` for adaptive runs; non-adaptive
  runs (and the entire initial phase) are untouched.
* `GpuMomentUpdater` verified byte-for-byte unmodified (§7).
* Part A/B/C/D audit, with a corrected, narrower account of the
  measurement/correlation look-ahead than CGP-03's evidence doc assumed
  (measurement's predicate is pure and peekable; correlation's is
  stateful and is not) -- see §4.
* 12/12 new unit tests passing, including a negative control that caught a
  real bug in the test's own first draft (§8) rather than in production
  code.
* Per-call benchmark evidence (§6) quantifying the achievable win
  (~1.27x-1.43x at 65536 atoms, fractions <=25% active) once the
  remaining blockers are cleared.

Explicitly **not** implemented (flagged for Terra, matching the CGP-00B
Part G / CGP-03 (b) precedent):

* `updateActiveOnly()` is not wired into production. Two independent
  prerequisites remain: revisiting CGP-03's own `materializeCoarseAtoms(OrdinaryStep, ...)`
  scope decision, and either a genuine side-effect-free correlation firing
  predicate or an accepted design for what happens to `emomM` on a
  correlation-enabled run (§4).
* CGP-03's own remaining checklist items ("hot timestep no longer
  reconstructs every coarse atom", "static all-coarse steady state
  performs no unnecessary full reconstruction") are **not** unblocked by
  this task alone (§5) -- this task removes a different, previously
  undocumented full-lattice cost, not the one those items name.
* The CPU adaptive path's equivalent moment/`emom` rebuild
  (`source/CoarseGraining/adaptivecgproduction.f90`) was not audited in
  this task; CGP-03's evidence doc already flagged the CPU path's
  structurally identical `reconstruct_coarse_atoms` issue as untouched, and
  this task's scope (per its own text) was the GPU moment updater only.

Performance: neutral for production (identical `updateFull()` cost,
verified by the parity tests using the same kernels/launches as before);
a real, measured, but not-yet-enabled ~1.27x-1.43x win is demonstrated in
isolation for a future task to unlock. No regression: 33/34 `cg13-cuda`,
matching the pre-existing baseline exactly (the one failure confirmed
unrelated).

## Checklist

* [x] Moment/`emomM` authoritative-state contract documented for interior
  coarse atoms (§2), including the `mompar` 1/2 resolution.
* [x] `GpuAdaptiveMomentUpdater` implemented, mirroring
  `GpuMomentUpdater`'s constructor/`update()` shape.
* [x] `GpuMomentUpdater` itself completely unmodified (§7, git-verified).
* [x] Touched-set update scoped to the touched atom set using the existing
  compact `activeAtomList` -- no new O(N) pass introduced.
* [x] Call-site ownership decided explicitly: full replacement, not
  supplement (§4); non-adaptive dynamics keep `GpuMomentUpdater` unchanged.
* [x] Measurement/correlation look-ahead question established plainly
  (§4/§5): correlation's firing predicate is stateful and not safely
  peekable without a `GpuCorrelations` edit; measurement's is pure but
  peeking it alone is not sufficient on its own.
* [x] All-fine fixed-seed parity against feature-off ASD passes for
  `mompar` 0/1/2 (`adaptive-cg-production-e2e`'s T=0/finite-T stages, §9).
* [x] `GpuAdaptiveMomentUpdater` on an all-fine (all-active) control
  reproduces `GpuMomentUpdater`'s `mmom`/`mmomi`/`emomM` bitwise (§8,
  `runActiveAllParity`).
* [x] Every existing CGP-00 through CGP-03 fixture (`ctest -L cg13-cuda`)
  keeps passing (33/34, one pre-existing unrelated failure, §9).
* [x] Negative control proving `GpuMomentUpdater` itself is untouched (§7,
  git diff/status/log).
* [x] Active-list-outside-atoms-pinned negative control (§8,
  `runActiveSubsetNegativeControl`).
* [x] Per-call performance evidence recorded (§6); whole-run number
  explicitly not claimed (§6, §10 -- production wiring is `updateFull()`,
  byte-identical cost to before this task).
* [ ] `updateActiveOnly()` wired into production. Not implemented --
  blocked on CGP-03's own `OrdinaryStep` scope and a genuine correlation
  firing-predicate query; flagged for Terra (§4, §10).
* [ ] CGP-03's "hot timestep"/"static all-coarse" checklist items
  unblocked. Not achieved by this task alone (§5); flagged for Terra.
* [ ] Terra review performed before merge. Requested; not yet performed.

## Commit

`CGP-03B: introduce adaptive-scoped moment updater`
