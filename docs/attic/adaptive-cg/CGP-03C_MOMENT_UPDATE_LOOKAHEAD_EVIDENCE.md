# CGP-03C — Wiring `GpuAdaptiveMomentUpdater::updateActiveOnly()` into production via a measurement look-ahead

Task: `docs/CGP_work.md` CGP-03C ("Wire `GpuAdaptiveMomentUpdater::updateActiveOnly()`
into production via a measurement look-ahead"). Branch `gpu_hip_cu_ab_cg`,
commit base `1f8c2fa7` (CGP-03B).

Status: **Complete, bounded scope.** `updateActiveOnly()` is now the
production default for adaptive runs; a pure, side-effect-free look-ahead
falls back to `updateFull()` on any step where a real consumer would
otherwise see stale `mmom`/`mmomi`/`emomM`. `GpuMomentUpdater`,
`GpuCorrelations`, and CGP-03's `materializeCoarseAtoms(OrdinaryStep, ...)`
atom-direction commit are all untouched. CGP-03's own "hot timestep"/
"static all-coarse" checklist items are **still not** unblocked by this
task alone — see §7.

## 0. Recap and the scoping decision that made this task tractable

CGP-03B built `GpuAdaptiveMomentUpdater::updateActiveOnly()` but could not
safely enable it in production: `GpuCorrelations::measure()`'s firing
predicate mutates its own bookkeeping (`curstep`/`t_cur`) as a side effect,
so it cannot be "peeked" one step ahead without editing shared correlation
code.

Before starting this task's own Part-1-style audit, Anders decided that
spin correlations (`do_sc`) do not need to be numerically correct/fresh
under adaptive CG at all — atomistic correlations are already expected to
be more approximate once part of the lattice is coarse-grained, and their
staleness under a touched-only moment update is an accepted consequence,
not a defect. (Recorded in project memory
`adaptive-cg-correlations-not-guaranteed`.) This removes `GpuCorrelations`
from the correctness picture entirely and leaves only measurement/output
consumers to audit.

## 1. Consumer audit: who actually needs fresh `mmom`/`mmomi`/`emomM`

Read directly from `source/gpu_files/measurement/gpuMeasurement.cpp`,
`source/gpu_files/measurement/cpuRestMeasurement.cpp`,
`source/gpu_files/c_helper.h`, `source/chelper.f90`, and
`source/Measurement/measurements.f90`:

| Consumer | Reads | Fresh required? | Predicate available? |
|---|---|---|---|
| `measureAverageMagnetization`/`measureBinderCumulant` | `emomMEnsembleSums`, built by `calculateEmomMSum()` from the full `emomM` | Yes | `GpuMeasurement::timeToMeasure(AverageMagnetization\|BinderCumulant, mstep)` — private, already pure |
| `measureSkyrmionNumber` | `gpuLattice.emomM` directly | Yes | `GpuMeasurement::timeToMeasure(SkyrmionNumber, mstep)` — private, already pure |
| `measureEnergy` | `gpuEnergies.energyM` only | **No** — never touches `emomM`/`mmom` | N/A |
| `measureAutocorrelation` | `gpuLattice.emom` (not `emomM`/`mmom`) | Governed entirely by CGP-03's own `materializeCoarseAtoms(OrdinaryStep, ...)` atom-direction commit, unconditional and untouched by this task | N/A for this task |
| `CpuRestMeasurement`/`FortranMeasurement` | Copy `emomM`/`emom`/`mmom`/`beff` to Fortran arrays whenever `alwaysCopy \|\| call_fortran_do_measurements(mstep)` | Yes, whenever the copy fires | `call_fortran_do_measurements()` → Fortran `do_measurements` (`source/Measurement/measurements.f90:32-97`) — a **pure** function of `mstep` and already-known sampling parameters (`avrg_step`, `cumu_step`, `do_tottraj`/`tottraj_step`, per-trajectory `traj_step(:)`, `logsamp`); confirmed by reading the subroutine body: a chain of `if(...) then; do_copy=1; return; end if`, no state mutated, no counters advanced |
| `GpuSDSimulation::printMdStatus` | Calls `gpuSim.copyToFortran()` (full `emom`/`emom2`/`emomM`/`mmom*`) on its own cadence, `mstep % ((rstep+nstep)/20) == 0`, or every step when `nstep<=20` | Yes, whenever it fires | Local, pure, already in the same file |
| Post-loop tail of `SDmphase` (final `measurement->measure()`, `flushMeasurements()`, unconditional `copyToFortran()`) | Whatever the last loop iteration left behind | Yes, unconditionally, once | Not a modulo condition — simplest correct answer is to force a full update on the last loop iteration |
| `diagnosticSnapshot()` | `blockState`/`stateAge`/`selectorScores`/`polarizationRatio`/energy/direction sums — **no** `mmom`/`mmomi`/`emomM` field | No | N/A |
| `GpuAdaptiveRuntime::activeAtomList()`/`activeAtomCount()` after a transition | Refreshed synchronously by compaction inside `publishProposedState()` (called from `advanceAdaptiveStep`) | N/A — this is the touched-set input, not a moment consumer | Re-fetch after `advanceAdaptiveStep` returns rather than reuse the pre-transition snapshot cached at the top of the step |
| `GpuCorrelations`/`FortranCorrelation` (`do_sc`) | `emomM`/`emom`/`mmom` on `sc_step`/`sc_sep` | **Explicitly not guaranteed** (Anders's decision, §0) | Not attempted |

Two findings worth calling out explicitly because they narrow the scope
below what CGP-03B's evidence doc anticipated:

* `measureEnergy` and `measureAutocorrelation` do not depend on
  `mmom`/`mmomi`/`emomM` at all — CGP-03B's Background section had grouped
  "measurement/correlation" together as a single concern; on closer reading
  only three of five `GpuMeasurement` types actually matter here.
* `fortran_do_measurements`'s underlying Fortran subroutine,
  `do_measurements`, takes every input as an explicit, `intent(in)`
  argument — it is not a hidden module-state read. The *wrapper*
  (`fortran_do_measurements` in `source/chelper.f90:463-470`) does forward
  live `InputData` module state into it, but that state is exactly what a
  real, already-initialized simulation has set up from its input file —
  the same state `CpuRestMeasurement`/`FortranMeasurement` already trust
  for their own copy-out gate today. Querying it one step ahead adds no
  new trust requirement beyond what production already relies on.

## 2. Transition correctness without a special case

A concern raised while scoping this task: does a fine→coarse or
coarse→fine transition on a selector-update step need to force a full
update, since the touched set changes mid-step?

Re-reading `GpuAdaptiveRuntime::activeAtomList()`/`activeAtomCount()`
(`gpuAdaptiveRuntime.hpp:432-439`) and `publishProposedState()`
(`gpuAdaptiveRuntime.cpp`) shows compaction runs synchronously inside
`publishProposedState()`, called from `advanceAdaptiveStep`, and refreshes
the host-mirrored `activeAtomList_`/`hostWorkCounts_` before that call
returns. Both transition directions are already correctly handled by
touched-only + a **freshly re-fetched** active list, with no forced full
update needed:

* **Coarse → fine**: `publishAdaptiveState` writes the newly-fine block's
  `atomDirection` directly (CGP-03 evidence doc §7). After compaction, that
  block's atoms are now *in* `activeAtomList`, so `updateActiveOnly()`
  correctly rebuilds their `mmom`/`emomM` from the fresh direction.
* **Fine → coarse**: the block is reclassified and its atoms drop out of
  `activeAtomList`; their `mmom`/`mmomi`/`emomM` correctly go stale at
  whatever value they held while still fine — exactly the intended
  contract (CGP-03B §2), no different from any other interior-coarse atom.

This is why `gpuSDSimulation.cpp` re-fetches
`gpuSim.gpuAdaptiveRuntime.activeAtomList()`/`activeAtomCount()` at the
moment-update call site rather than reusing the pre-transition snapshot
`advanceAdaptiveStep` caches at the top of the step for the Depondt
integrator (`gpuSimulation.cpp:1213-1214`) — using the stale snapshot here
would silently miss a coarse→fine block's freshly-written atoms.

## 3. Design

### A. `GpuMeasurement::needsFreshMoments(mstep)`

New public, `const`, side-effect-free method
(`source/gpu_files/measurement/gpuMeasurement.{hpp,cpp}`):

```cpp
bool GpuMeasurement::needsFreshMoments(size_t mstep) const {
   --mstep;
   return timeToMeasure(MeasurementType::AverageMagnetization, mstep) ||
          timeToMeasure(MeasurementType::BinderCumulant, mstep) ||
          timeToMeasure(MeasurementType::SkyrmionNumber, mstep);
}
```

Mirrors `measure()`'s own `--mstep; timeToMeasure(...)` sequence exactly, so
a caller passing the same raw `mstep` it would pass to `measure()` gets the
same answer `measure()` would compute for those three types. No existing
method's behavior changes; this is purely additive.

### B. `adaptiveMomentUpdateNeedsFullLattice(...)`

New free function in `gpuAdaptiveMomentUpdater.{hpp,cpp}` (placed there,
not in `gpuSDSimulation.cpp`, specifically so it is unit-testable and
decoupled from `GpuSimulation`'s internals — see §5):

```cpp
bool adaptiveMomentUpdateNeedsFullLattice(GpuMeasurement* gpuMeasurement,
                                          std::size_t currentStep, std::size_t lastStep,
                                          std::size_t nstep, std::size_t rstep) {
   if(currentStep >= lastStep) return true;
   const std::size_t nextStep = currentStep + 1;
   if(call_fortran_do_measurements(static_cast<int>(nextStep))) return true;
   if(gpuMeasurement && gpuMeasurement->needsFreshMoments(nextStep)) return true;
   if(nstep <= 20) return true;
   return (nextStep % ((rstep + nstep) / 20)) == 0;
}
```

`gpuMeasurement` may be null (`MeasurementFactory` chose
`FortranMeasurement`, i.e. `do_cuda_measurements != 'Y'`); the Fortran
oracle alone then matches that backend's own existing freshness guarantee
(§1's audit — `FortranMeasurement`'s own copy-out gate is the same
`call_fortran_do_measurements` predicate, and nothing about that backend's
skyrmion path bypasses it the way `GpuMeasurement::measureSkyrmionNumber`
does), so no special-casing is needed for it.

`GpuMeasurement` is forward-declared in `gpuAdaptiveMomentUpdater.hpp`
rather than `#include`d, so every other consumer of that header keeps its
pre-CGP-03C dependency footprint.

### C. Call-site wiring (`gpuSDSimulation.cpp`, `SDmphase`)

```cpp
if(gpuSim.adaptiveEnabled()) {
   if(adaptiveMomentUpdateNeedsFullLattice(gpuMeasurementForLookahead, mstep, rstep + nstep,
                                           nstep, rstep)) {
      adaptiveMomUpdater.updateFull();
   } else {
      adaptiveMomUpdater.updateActiveOnly(gpuSim.gpuAdaptiveRuntime.activeAtomList(),
                                          gpuSim.gpuAdaptiveRuntime.activeAtomCount());
   }
} else {
   momUpdater.update();
}
```

`gpuMeasurementForLookahead` is `dynamic_cast<GpuMeasurement*>(measurement.get())`,
computed once, non-owning. Non-adaptive dynamics (`momUpdater.update()`)
and the entire initial phase (`SDiphase`, confirmed by reading the whole
function: it never calls `gpuSim.adaptiveEnabled()` and always uses the
plain, non-adaptive path) are untouched.

### D. What is deliberately not touched

`GpuMomentUpdater` (still byte-for-byte unmodified — see §6),
`GpuCorrelations` (per Anders's scoping decision, §0), and
`materializeCoarseAtoms(OrdinaryStep, ...)` in `gpuSimulation.cpp` (still
unconditional — CGP-03's own remaining scope, not this task's; see §7).

## 4. Correctness suite

Same out-of-tree build as CGP-03B, rebuilt incrementally:
`build_cgp03b_cuda_fp32` (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`,
`CMAKE_BUILD_TYPE=Release`), CUDA 13.3.73, driver 610.57.04, 2x NVIDIA RTX
A4000 (both idle, 0% utilization, <110 MiB used immediately before every
run below).

`ctest -L cg13-cuda` (34 tests, same set as CGP-03B): **33/34 pass**,
identical to the CGP-03B baseline. The one failure,
`adaptive-cg-production-e2e` (`assert abs(float_metric(fft.stdout,
"coarse_dipole")) > 0.0`), is the same pre-existing, unrelated failure
CGP-01 through CGP-03B all documented — confirmed by identical output up
to that point (every `CG-10.5`/`RCG-03` stage before it passes). Notably,
`adaptive-cg-moving-*` and `adaptive-cg-transition-ownership-invariants`
— which exercise real, repeated fine↔coarse transitions — now genuinely
run through `updateActiveOnly()` with the freshly-refreshed active list
(§2), not merely `updateFull()` as in CGP-03B, and still pass, which is
real evidence for §2's transition-correctness argument beyond the
citation-based reasoning alone.

## 5. `GpuMomentUpdater` still untouched — negative control

```
$ git diff --stat -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp
$ git status --short -- source/gpu_files/gpuMomentUpdater.cpp source/gpu_files/gpuMomentUpdater.hpp
```

Both empty, same as CGP-03B's own check — this task adds no further edits
to either file.

## 6. Tests

### A. `gpu_adaptive_moment_updater_tests` (unchanged from CGP-03B)

Still 12/12 combinations passing after this task's refactor of
`gpuAdaptiveMomentUpdater.cpp` moved code around it (rebuilt and rerun to
confirm — see §4).

### B. Why an isolated unit test for `adaptiveMomentUpdateNeedsFullLattice()`
was not written, and what was done instead

`adaptiveMomentUpdateNeedsFullLattice()` was deliberately relocated out of
`gpuSDSimulation.cpp` (where it started as a `static`, untestable function)
into `gpuAdaptiveMomentUpdater.{hpp,cpp}` specifically to make it callable
from a test binary. Attempting to actually write that test surfaced a real
constraint: `call_fortran_do_measurements()` reads live Fortran `InputData`
module state through `fortran_do_measurements` (`source/chelper.f90:463-470`),
which is only meaningfully populated by a real, input-file-driven
simulation setup — not something a bare host-only test binary can safely
fake (there is no C-interop setter for those module variables, only a
one-time, one-directional pointer *export* from Fortran to C++ via
`FortranData::setFlagPointers`, which a test cannot redirect to control
what `fortran_do_measurements` itself reads). Constructing a real
`GpuMeasurement` to test `needsFreshMoments()` directly is similarly heavy
(its constructor requires a working `FortranData` pointer set and, for
`do_avrg_proj`/`do_cumu_proj` combinations, per-atom type/channel arrays).

Given that, verification instead used a real, running simulation (§7 in
CGP-03B's style: a per-call microbenchmark was preferred there for the
analogous reason; here, an end-to-end run is the analogous choice) rather
than a synthetic host-only fixture:

* A copy of the existing `finite_temperature_mixed` e2e fixture
  (`tests/coarse_graining/e2e/finite_temperature_mixed/inpsd.dat` — already
  part of the passing `adaptive-cg-production-e2e` test, 48 atoms, 1 fine +
  4 interface + 1 coarse block, `mompar` at its default of 0) was run in a
  scratch directory with `do_avrg Y avrg_step 2/5`, `do_cumu Y cumu_step
  3/7`, `do_tottraj Y tottraj_step 5/11` added and `Nstep` raised to 12 and
  then 60, using the real `sd.f95.cuda` binary built from this branch.
* **Sanity check**: the real implementation (`adaptiveMomentUpdateNeedsFullLattice`
  as shipped) produced `averages.*.out`/`cumulants.*.out`/`moment.*.out`
  files **MD5-identical** to a build with the function's body replaced by
  `return true;` (i.e., forced to behave exactly like CGP-03B's
  `updateFull()`-always wiring) — `b7c66a5c.../0a93e49b.../40dd548b...`
  matched at `Nstep=12`; `7dc4e2f6.../9359ccbe.../6d2382a7...` matched at
  `Nstep=60`. This directly demonstrates the look-ahead's `true` branches
  fire correctly and `updateActiveOnly()` reconstructs bit-identical output
  on every step it is used, for a real production run exercising avrg,
  cumu, and trajectory output at staggered, non-aligned intervals.
* **Attempted negative control**: the function's body was also replaced
  with `return false;` (never force a full update — i.e. `updateActiveOnly()`
  unconditionally, exactly the configuration CGP-03B's evidence doc warned
  against shipping). At both `Nstep=12` and `Nstep=60` this **also**
  produced MD5-identical output to the real implementation. Investigating
  why (rather than treating this as a pass): `AdaptiveCG:
  reduced_atom_updates` in this fixture's own stdout confirms exactly 8 of
  48 atoms are genuinely excluded from `activeAtomList` every step (the
  single coarse block), so the exclusion is real; the non-discrimination is
  because that block's `coarseDirection` does not move measurably within
  this fixture's duration (`delta_t=0.1fs × 60 steps = 6fs` of physical
  time, `mompar=0` so moment magnitude is provably unaffected regardless
  per CGP-03B §1's double-swap-identity finding, and `damping=0.05` at
  `T=10K` starting from a near-saturated configuration) — not because
  the excluded atoms' staleness is harmless in general. This is recorded
  honestly as a **non-discriminating** empirical negative control, not a
  successful one: it does not by itself prove skipping the look-ahead is
  safe, only that this specific fixture's coarse-block dynamics are too
  slow to reveal the difference within the tested window. The real,
  discriminating negative control for the underlying exclusion mechanism
  is CGP-03B's own `runActiveSubsetNegativeControl` unit test (synthetic,
  sentinel-seeded, and does detect a skipped write) — this task's
  contribution is the audit-and-citation argument (§1) for *when* that
  mechanism is safe to invoke, not a second, independent empirical proof
  of the exclusion mechanism itself.
* The source was reverted to the real implementation and rebuilt; the
  reverted build was re-verified bit-identical to the saved `Nstep=60`
  reference (`7dc4e2f6...`/`9359ccbe...`/`6d2382a7...`), and
  `gpu_adaptive_moment_updater_tests` was rerun (still 12/12 pass) before
  proceeding.

## 7. CGP-03's checklist items: still not unblocked

Same conclusion as CGP-03B's evidence doc, now with the second precondition
(the correlation look-ahead) removed but the first unchanged:
`materializeCoarseAtoms(OrdinaryStep, ...)` in `gpuSimulation.cpp` still
runs unconditionally, every step, committing `coarseDirection` into every
interior-coarse atom's `emom2` — this task's audit and required work were
scoped to the moment-update side only (per CGP-03B's own Background, which
frames this task as unblocking `GpuMomentUpdater`, not
`materializeCoarseAtoms`). Reopening that call site is CGP-03's own
remaining scope: it would mean restoring the `Output`/`Diagnostics`
materialization reasons CGP-03 deliberately dropped (`docs/CGP-03_LAZY_MATERIALIZATION_EVIDENCE.md`
§5) once the free ride from the unconditional commit is gone, and is left
for a further task rather than folded in here.

## 8. Performance evidence

This task changes *when* `GpuAdaptiveMomentUpdater::updateActiveOnly()` is
called, not its cost — CGP-03B's own per-call benchmark (65536 atoms, 1
ensemble, discarded warm-up, ABBA sampling) is unchanged and is reused
here rather than re-measured:

| active fraction | mompar | full (median us) | active-only (median us) | speedup |
|---|---|---|---|---|
| 1.25% | 0/1/2 | 8.09 / 11.05 / 10.83 | 5.93 / 7.94 / 7.81 | 1.37-1.39x |
| 6.25% | 0/1/2 | 8.25 / 10.81 / 10.81 | 5.78 / 7.90 / 7.90 | 1.37-1.43x |
| 25% | 0/1/2 | 8.15 / 10.90 / 10.83 | 6.36 / 8.51 / 8.51 | 1.27-1.28x |

(Full table with MAD: `docs/CGP-03B_ADAPTIVE_MOMENT_UPDATER_EVIDENCE.md` §6.)

What CGP-03C adds is realizing this per-call win on whatever fraction of
ordinary steps `adaptiveMomentUpdateNeedsFullLattice()` allows — which, for
a production run with realistic (not every-step) `avrg_step`/`cumu_step`/
`skyno_step`/output cadences, is the large majority of steps. A whole-run
wall-clock number was not (re-)measured for this task: the tiny e2e
fixture used for §6's testing (48 atoms) is dominated by kernel-launch and
Fortran I/O overhead rather than the moment-update kernels, so a whole-run
timing there would not be representative, and re-running CGP-07's own
65536-atom whole-run sweep was judged out of scope for a task whose
performance claim is fully supported by an unchanged, already-measured
per-call number now realized on real production steps — matching CGP-03's
own §10 precedent of reporting a per-call number rather than inferring a
whole-run figure it cannot support.

## 9. Outcome summary

Implemented:

* `GpuMeasurement::needsFreshMoments(mstep)` — small, additive, pure query.
* `adaptiveMomentUpdateNeedsFullLattice(...)` — the production look-ahead,
  relocated to `gpuAdaptiveMomentUpdater.{hpp,cpp}` for testability.
* `gpuSDSimulation.cpp` now calls `GpuAdaptiveMomentUpdater::updateActiveOnly()`
  as the default for adaptive runs, falling back to `updateFull()` exactly
  when the look-ahead requires it.
* Deliberately excludes `GpuCorrelations`/`do_sc` from the correctness
  contract (Anders's decision, §0) and does not touch
  `materializeCoarseAtoms(OrdinaryStep, ...)` (§7, CGP-03's own remaining
  scope).
* `GpuMomentUpdater` verified still byte-for-byte unmodified (§5).
* 33/34 `cg13-cuda` (unchanged pre-existing failure); real transition/moving-
  interface fixtures now genuinely exercise `updateActiveOnly()` end to end
  and pass.
* A real-run sanity check (MD5-identical output vs. an always-full build)
  and an honestly-reported non-discriminating negative-control attempt
  (§6B), with the actual discriminating negative control for the exclusion
  mechanism itself remaining CGP-03B's own synthetic unit test.

Explicitly **not** implemented / not achieved:

* CGP-03's "hot timestep no longer reconstructs every coarse atom" and
  "static all-coarse steady state performs no unnecessary full
  reconstruction" checklist items remain open (§7) — this task closes the
  moment-update half of CGP-03B's blocker list, not the atom-direction
  half, which was never this task's scope.
* No whole-run wall-clock number for this specific change (§8) — the
  per-call number is unchanged from CGP-03B and is judged sufficient
  support for the claim being made.
* No isolated host-only unit test for the look-ahead's Fortran-oracle
  branch (§6B) — blocked on the same live-Fortran-module-state constraint
  that makes it a real integration concern, not a unit-testable one, with
  the real-run evidence in §6B substituting for it.

Performance: the per-call win demonstrated in CGP-03B (~1.27x-1.43x at
65536 atoms, <=25% active) is now realized in production on every ordinary
step; no whole-run figure claimed (§8). No regression: 33/34 `cg13-cuda`,
matching CGP-03B's baseline exactly, with the one failure confirmed
unrelated and pre-existing.

## Checklist

* [x] Every `GpuMeasurement`/output consumer's `mmom`/`mmomi`/`emomM`
  dependency audited (§1), narrower than CGP-03B's Background anticipated
  (energy and autocorrelation excluded).
* [x] `needsFreshMoments()` added, pure/additive, no existing behavior
  changed.
* [x] Look-ahead combinator implemented and unit-testable in principle
  (relocated for that reason), though the Fortran-oracle branch could not
  actually be isolated from live module state (§6B).
* [x] `updateActiveOnly()` wired into production as the default for
  adaptive runs, falling back to `updateFull()` per the look-ahead.
* [x] Transition correctness verified without a special case, by
  re-fetching a fresh post-transition active list (§2), and confirmed by
  real moving-interface/transition fixtures continuing to pass (§4).
* [x] `GpuMomentUpdater` still completely unmodified (§5, git-verified).
* [x] `GpuCorrelations`/`materializeCoarseAtoms(OrdinaryStep, ...)`
  deliberately not touched (§0, §7).
* [x] Every existing CGP-00 through CGP-03B fixture (`ctest -L cg13-cuda`)
  keeps passing (33/34, one pre-existing unrelated failure, §4).
* [x] Real-run sanity check: look-ahead output MD5-identical to an
  always-full build, at two run lengths (§6B).
* [ ] Negative control proving the look-ahead is load-bearing on a real
  fixture. Attempted; not discriminating for the fixture used, with a
  physically-grounded explanation recorded rather than claimed as a pass
  (§6B). The exclusion mechanism itself remains proven by CGP-03B's own
  synthetic unit test.
* [x] Performance evidence recorded; whole-run number explicitly not
  claimed, with reasoning (§8).
* [ ] CGP-03's "hot timestep"/"static all-coarse" checklist items
  unblocked. Not achieved by this task (§7); still flagged for a further
  task (the `materializeCoarseAtoms(OrdinaryStep, ...)` scope revisit).
* [ ] Terra review performed before merge. Requested; not yet performed.

## Commit

`CGP-03C: wire touched-only adaptive moment update into production`
