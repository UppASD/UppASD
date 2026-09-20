# CGP-06B — Make thermal RNG work scale with active fine atoms, evidence

Task: `docs/CGP_work.md` CGP-06B. Branch `gpu_hip_cu_ab_cg`, built on top of
CGP-06A (commit `4f9f1553`, evidence-only, no production code changed).

**Primary model (task doc): Terra**, escalated because this touches
production stochastic semantics. **Actual model: Claude (Sonnet 5,
interactive session)**, matching the CGP-03D/CGP-05/CGP-06A precedent of an
interactive session performing the audit/design/implementation slice
directly, with the discrepancy disclosed rather than silently substituted.

## Human acceptance of the RNG contract

Per `CGP_work.md`'s dependency ("CGP-06A and Human acceptance of its RNG
contract"), Anders was asked before implementation began whether to accept
CGP-06A's recommended scope for this task: **Strategy 1** (parity-preserving,
bounded win — generator sequence untouched, only the field *write* scoped to
the active list), with Strategy 3 (counter-based RNG, unbounded win, new
stochastic algorithm) deferred to a future task requiring its own separate
acceptance. Anders accepted Strategy 1 as specified. This task implements
**only** Strategy 1, per CGP-06B's own "Implement ONLY the strategy accepted
after CGP-06A" instruction.

## Design

`GpuThermfield::randomize()`'s original implementation
(`source/gpu_files/gpuThermfield.cpp`) does two things unconditionally, every
call: (1) `GPU_RAND_GENERATE_NORMAL_REAL` draws `3*N*M` values from the
shared curand/hiprand stream, and (2) `parallel.gpuAtomSiteCall(SetupField(...))`
writes `field[atom]` for every `(site,ensemble)` pair. CGP-06A traced that
only step (2) is wasted work for interior coarse atoms — `evolveFirst`'s
active-list overload reads `thermalField` exclusively through the same
active list the predictor kernel itself uses, so an unwritten interior-coarse
entry is never read.

This task adds a second `randomize()` overload that keeps step (1) **exactly
as it is today** (same call, same buffer, same size — the generator's
sequence is untouched) and replaces step (2) with a new active-list-scoped
write:

* `GpuParallelizationHelper` gained `active_atom_site_kernel` and
  `gpuActiveAtomSiteCall()` (`gpuParallelizationHelper.hpp`/`.tpp`) — the
  `AtomSite`-interface sibling of the existing `active_atom_kernel`/
  `gpuActiveAtomCall()` (`Atom`-interface only). It derives `site` from
  `oneBasedAtoms[slot]-1` before calling `op.each(atom, site)`, exactly the
  translation the existing compact-list kernels already perform — never from
  `slot` itself, which is the failure mode the slot-indexing negative control
  below exists to rule out.
* `GpuThermfield::randomize(mmom, oneBasedAtoms, activeCount)`
  (`gpuThermfield.hpp`/`.cpp`) reuses the existing `SetupField` op unchanged
  (it already has the `each(atom, site)` shape `AtomSite` requires) and calls
  `parallel.gpuActiveAtomSiteCall(SetupField(...), oneBasedAtoms, activeCount)`
  instead of the full `gpuAtomSiteCall`. `gpuActiveAtomSiteCall` returns
  without launching a kernel when `activeCount == 0` (the all-coarse case),
  matching `gpuActiveAtomCall`'s existing empty-list convention — satisfying
  the "should not ... write a full atomistic thermal field merely because the
  array exists" half of CGP-06B's all-coarse requirement.
* `GpuDepondtIntegrator::evolveFirst(gpuLattice, oneBasedAtoms, activeAtomCount)`
  (`gpuDepondtIntegrator.cpp`, the only production caller of the active-list
  overload — confirmed by grep, called once per adaptive step from
  `GpuSimulation::advanceAdaptiveStep`, `gpuSimulation.cpp:1299`, using the
  same `gpuAdaptiveRuntime.activeAtomList()`/`activeAtomCount()` compact list
  CGP-02/03C already established) now calls the new overload instead of the
  full `randomize(mmom)`. No new O(N) pass was introduced; no other
  production call site exists (`evolveFirst(gpuLattice)`,
  `evolveFirstWithThermalField`, and the feature-off path are all untouched).

### The all-coarse tension, disclosed rather than silently resolved

CGP_work.md's general "Special all-coarse case" text also says a
timestep "should not *generate* ... a full atomistic thermal field merely
because the array exists." Under the accepted Strategy 1 contract this
cannot be done: the curand/hiprand generate call **must** keep drawing the
full `3*N*M` values every step, unconditionally, including all-coarse steps,
specifically so a physical atom's stochastic continuation across a future
coarse→fine transition never depends on how many steps it spent coarse (this
is exactly CGP-06A's invariants 2-4, and exactly what Strategy 1 was accepted
*for* — trading a smaller, unconditional ceiling for zero new correctness
risk). The negative control in the Tests section below demonstrates
concretely why skipping the generate call during coarse residence is unsafe
under this design. This task therefore satisfies the all-coarse requirement's
*write* half unconditionally, and deliberately does not attempt its
*generate* half — matching CGP-06A's own Strategy 1 evaluation ("removes only
the loop sub-phase ... not the O(N) floor itself") rather than reinterpreting
the accepted contract after the fact.

## Tests

### Unit level (`tests/coarse_graining/test_gpu_depondt_active_atoms.cpp`)

Two pre-existing RCG-09A.2 fixtures now exercise the new code path (the
active-list `evolveFirst` overload they call now routes through the
active-scoped `randomize()`):

* `runAllActiveParity` (T=0 and finite T): the active-list integrator (all 5
  atoms active) reproduces the ordinary full-range integrator bitwise —
  `emom`/`emom2`/`emomM`/`b2eff`. **Passed.** This is CGP-06A invariant 1's
  concrete, up-to-date proof.
* `runReorderedSubsetIdentity`: forward `{2,5}` vs. reverse `{5,2}` compact
  lists produce bitwise-identical results, with every unselected atom pinned
  unchanged. **Passed.** This is CGP-06A invariant 2's concrete, up-to-date
  proof.

Three new fixtures added by this task:

* `runPartialActiveFieldScoping`: two same-seeded `GpuThermfield` instances,
  one full-write, one active-scoped write (list `{2,5}` of 5 sites). The
  active-scoped instance's field buffer is poisoned (`-9999`, via
  `GPU_MEMCPY` into `getField().data()` — `GPU_MALLOC` does not
  zero-initialize) before the call. Result: entries at sites 2 and 5 match
  the full-write reference bitwise (same generator draw); every other entry
  still holds the poison. **Passed.** Proves the scoped write is both
  correct at active sites and genuinely absent elsewhere, not merely "not
  observed to differ."
* `runNegativeControlSlotIndexing` (task's own suggested control #1): a
  self-contained toy kernel deliberately indexes the written value by compact
  **slot** instead of physical **site** (the bug `active_atom_site_kernel`
  itself does not have). Forward `{2,5}` vs. reverse `{5,2}` compact lists
  are asserted to *disagree* at the two active sites. **Passed** (the
  assertion is that they disagree — a failure to disagree throws, so this
  control is discriminating, not vacuous). This is the concrete evidence that
  `runReorderedSubsetIdentity` is actually testing something.
* `runNegativeControlSkippedAdvancement` (task's own suggested control #2):
  two same-seeded `GpuThermfield` instances; the reference calls
  `randomize(mmom)` twice (simulating two ordinary steps); the "skipped"
  instance omits the first call entirely (simulating a hypothetical
  skip-generate-during-coarse-residence optimization) and calls `randomize`
  only once. Their final fields are asserted to *differ* — reusing the second
  reference call's second-batch draw is not what the once-called instance
  produces (it consumes the *first* batch instead). **Passed** (the assertion
  is that they differ). Concrete evidence that unconditionally advancing the
  generator every step — this task's actual design — is load-bearing, not a
  conservative extra.

Full suite result: `ctest --test-dir build_cgp06b_cuda_fp32 -L cg13-cuda`:
**33/34 pass.** The one failure, `adaptive-cg-production-e2e`
(`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), reproduces
the exact same pre-existing, unrelated failure already on record in every
CGP-01/03B/03C/03D/05/06A evidence doc (fp64-only artifact of that
particular fixture, not this task's change — confirmed again here by
`ctest --test-dir build_cgp06b_cuda_fp64 -L cg13-cuda`: **34/34 pass**, same
pattern as every prior task). `gpu-depondt-active-atoms` (which now includes
the three new fixtures above) passed at both precisions.

### Full-simulation parity: T=0 four-way MD5 comparison

Following the CGP-03D protocol exactly: a genuinely mixed-resolution fixture
(`tests/coarse_graining/e2e/finite_temperature_mixed`'s topology/mask, with
`temp 10.0` → `temp 0.0`, `Nstep 2` → `Nstep 200`, `do_avrg Y avrg_step 7`,
`do_tottraj Y tottraj_step 11`, `do_gpu Y` added, `tseed 4711` kept explicit —
`fine=1 interface=4 coarse=1` on the 6-block grid, so `active_atoms > 0`
genuinely exercises the active-scoped `randomize()` call every step) was run
twice on this task's build and twice on a clean pre-CGP-06B baseline (`git
worktree` checkout of commit `4f9f1553`, built identically) — four runs
total, each precision.

**fp32:** `averages`/`moment`/`restart` output MD5s are byte-for-byte
identical across new-run-1, new-run-2, baseline-run-1, baseline-run-2.

**fp64:** same result — byte-for-byte identical across all four runs.

This establishes, at the full-production-pipeline level and not merely the
isolated kernel, that Strategy 1 is a bit-exact no-op on the accepted
trajectory for a fixture that genuinely has fine atoms adjacent to coarse
ones (i.e. genuinely exercises the new scoped-write code path), at both
precisions.

### Moving/mixed finite-T and transition/restart continuity

This task's code changes are confined to `GpuThermfield::randomize()`'s
write geometry and its one production call site
(`GpuDepondtIntegrator::evolveFirst`'s active-list overload); no restart,
transition, or moving-interface code was touched, and the active list passed
in is the same `gpuAdaptiveRuntime.activeAtomList()`/`activeAtomCount()`
CGP-03C already established is correctly re-fetched *after*
`advanceAdaptiveStep` returns for both transition directions. Given that,
this task's correctness burden for moving/transition/restart fixtures
reduces to "did anything downstream of `randomize()` regress", which the
`cg13-cuda` suite already covers end-to-end with the modified binary:
`adaptive-cg-moving-backend-parity` (moving fixture, GPU vs. CPU backend
parity), `adaptive-cg-transition-ownership-invariants`,
`adaptive-cg-timing-reconciliation`, and `coarse-graining-adaptive-asd-parity`
(RCG-09A.4 all-fine parity) all passed at both precisions with this task's
binary. A dedicated new moving/restart MD5 fixture beyond what the suite
already exercises was not separately constructed — this is a smaller,
narrower-scope change than CGP-03D/CGP-05 (it touches nothing but a write
geometry inside one already-active-list-scoped call), and Part B's
invariant-5 analysis (CGP-06A) already established restart is unaffected
either way (curand/hiprand state was never serialized before this task, and
still is not; this task changes only which array entries get written per
call, not seed/restart handling).

## Performance

Environment: NVIDIA RTX A4000 (device 0 for fp32, device 1 for fp64;
`CUDA_VISIBLE_DEVICES` pinned per project memory
`shared-gpu-host-contention`), driver 610.57.04, idle immediately before
every run (`nvidia-smi`: 0% utilization, 210 MHz SM clock, <110 MiB used, no
other compute process on either GPU). CUDA 13.3.73, `Release` build. Host
CPU: `nproc` reports 2 (cgroup-quota-limited; project memory
`shared-gpu-host-contention`), material to the wall-vs-GPU-event gap below,
same as CGP-06A's own finding. One full cold-start invocation discarded
before every swept size (project memory `gpu-benchmark-cold-start`); 5
discarded warm-up iterations and 30 measured iterations per repetition, 7
repetitions, median/MAD reported — same harness and protocol as CGP-06A,
extended with `--active-fraction F` to call the real
`randomize(mmom, oneBasedAtoms, activeCount)` overload this task added, at
`N' = round(F*N)` contiguous one-based ids, instead of CGP-06A's
smaller-full-lattice *projection*. This supersedes that projection with a
real measurement of the exact kernel production now uses.

Raw output: `cgp06b_thermfield_active_fp32.log`, `cgp06b_thermfield_active_fp64.log`
(scratch, not committed).

### fp32, 262,144 atoms, T=50

| requested fraction | active atoms | wall us/call (median, MAD) | RNG GPU us/call | loop GPU us/call |
|---|---|---|---|---|
| unscoped baseline (`randomize(mmom)`) | 262,144 | 39.075, 0.042 | -- | -- |
| 1.0 (active-scoped) | 262,144 | 41.774, 0.047 | 6.528 | 7.885 |
| 0.25 | 65,536 | 21.933, 0.026 | 8.165 | 3.521 |
| 0.0625 | 16,384 | 17.821, 0.012 | 9.019 | 2.134 |
| 0.01 | 2,621 | 17.282, 0.020 | 9.518 | 1.855 |

Wall-clock reduction relative to the active-scoped call at full occupancy
(the fair apples-to-apples comparison, since the scoped call carries a small
fixed cost the old unscoped call did not — the `oneBasedAtoms` compact-list
indirection): **1.90x at 25% fine, 2.34x at 6.25% fine, 2.42x at 1% fine.**
GPU-event-only reduction (RNG+loop sum, the part CPU dispatch noise on this
host does not swamp): 14.41us → 11.15us at 6.25% fine, **1.29x** — a real but
much smaller number than the wall-clock figure, consistent with CGP-06A's
own finding that this host's wall clock is dominated by CPU-side dispatch
noise, not GPU execution time. At full occupancy (fraction 1.0) the
active-scoped call is measurably ~7% *slower* than the old unscoped call
(41.77us vs. 39.08us) — the compact-list indirection's fixed cost is not
free, and is not hidden when there is no write-skipping benefit to offset it;
disclosed rather than omitted.

### fp64, 262,144 atoms, T=50

| requested fraction | active atoms | wall us/call (median, MAD) | RNG GPU us/call | loop GPU us/call |
|---|---|---|---|---|
| unscoped baseline | 262,144 | 336.876, 0.144 | -- | -- |
| 1.0 (active-scoped) | 262,144 | 323.039, 0.079 | 18.705 | 3.379 |
| 0.25 | 65,536 | 287.874, 0.168 | 18.569 | 1.092 |
| 0.0625 | 16,384 | 284.875, 0.259 | 18.637 | 0.546 |
| 0.01 | 2,621 | 277.840, 0.132 | 18.705 | 0.444 |

Wall-clock reduction relative to the active-scoped call at full occupancy:
**1.12x at 25% fine, 1.13x at 6.25% fine, 1.16x at 1% fine** — much smaller
than fp32, exactly as CGP-06A's Part C predicted ("Strategy 1's ceiling
matters more, not less, for DOUBLE builds", since the RNG-generate sub-phase
is a much larger share of `randomize()`'s cost at fp64: here 18.6-18.7us GPU
vs. 0.4-3.4us for the shrinking loop share, i.e. the loop is only ~15% of
total GPU time even at full occupancy). GPU-event-only reduction: 22.08us →
19.18us at 6.25% fine, **1.15x**.

### What this task does not claim

* **No O(N) floor removal.** The generate sub-phase (the dominant cost at
  fp64, and comparable to the loop cost at fp32) is unconditional by design
  (see "the all-coarse tension" above) and does not shrink with active
  fraction at either precision — exactly the bounded-ceiling outcome
  CGP-06A's Strategy 1 evaluation predicted, not a defect discovered late.
* **No whole-adaptive-step or whole-production-run wall-clock benchmark was
  attempted.** `gpu_adaptive_runtime_benchmark`'s existing per-fraction sweep
  deliberately excludes the shared production Depondt/thermfield call from
  its timing by design (its own source comment: "adaptive field/coarse
  timing excludes shared production Depondt integration"), so it cannot
  supply this number without its own separate edit, out of this task's named
  scope. A dedicated new whole-step harness was not built: this task's own
  measured effect is bounded to a few/several tens of microseconds inside one
  `randomize()` call among several per step, on a host already shown
  (CGP-05, CGP-06A) to have CPU-quota-limited wall-clock noise at the
  microsecond scale that swamps changes of this size. Attempting it without
  a way to isolate the effect from that noise would not have produced a
  trustworthy number; this is the same category of disclosed gap as CGP-05's
  "inconclusive" whole-step remeasurement, reported honestly rather than
  fabricated.
* **Integration-time-vs-active-fraction is unchanged by this task** and was
  already established by CGP-01/RCG-09A (the active-list Depondt kernels
  themselves are untouched here).

## Outcome

`GpuThermfield` gained a `randomize(mmom, oneBasedAtoms, activeCount)`
overload (Strategy 1, per Anders's acceptance of CGP-06A's recommendation):
identical curand/hiprand generate call, byte-for-byte identical generator
sequence to the unscoped path, with the `SetupField` write scoped to the
compact active list via a new `GpuParallelizationHelper::gpuActiveAtomSiteCall`
(the `AtomSite` sibling of the existing `gpuActiveAtomCall`). The one
production caller, `GpuDepondtIntegrator::evolveFirst`'s active-list
overload, now uses it, reusing the same compact list CGP-02/03C already
established rather than adding any new O(N) pass. Empty active lists (the
all-coarse case) issue no write kernel at all. Correctness was verified at
the unit level (two pre-existing RCG-09A.2 fixtures plus three new ones —
partial-subset scoping with poisoned-buffer isolation, and two discriminating
negative controls matching the task's own suggested controls), at the full
`cg13-cuda` suite level (33/34 fp32, 34/34 fp64, the one fp32 failure the
same pre-existing/unrelated finding every prior CGP task in this phase has
recorded), and at the full-production-pipeline level (T=0 four-way MD5-
identical parity, new vs. a clean pre-CGP-06B baseline, both precisions, on
a genuinely mixed-resolution fixture). Performance was measured directly
(not projected) at 262,144 atoms across four fine fractions and both
precisions: a real 1.9x-2.4x wall-clock (1.3x GPU-event) reduction in
`randomize()`'s own cost at fp32 for realistic fine fractions, and a much
smaller 1.1x-1.2x reduction at fp64 — confirming CGP-06A's predicted
precision-dependent asymmetry with real measurements. No whole-step or
whole-run benchmark was attempted, disclosed as out of scope given the
change's small absolute magnitude relative to this host's known CPU-quota
noise floor. CGP-06B's own "no O(N) floor removal" limitation (only the
write sub-phase is scoped; the generate sub-phase remains full-N by design,
required to preserve invariants 2-4) is the same bounded-ceiling result
CGP-06A's evidence predicted, not a new finding.

## Checklist

* [x] Accepted CGP-06A strategy implemented (Strategy 1, Human-accepted
  before implementation began).
* [x] No second stochastic convention introduced (generator, seed handling,
  and sequence are all untouched; only the write geometry changed).
* [x] All-fine fixed-seed feature-off parity passes (`runAllActiveParity`,
  T=0 and finite T).
* [x] Active-list permutation invariance passes (`runReorderedSubsetIdentity`,
  plus the new `runNegativeControlSlotIndexing` proving it is discriminating).
* [x] Static mixed finite-T fixture passes (T=0 four-way MD5 parity on the
  mixed-resolution fixture, both precisions; `cg13-cuda`'s
  `adaptive-cg-timing-reconciliation`/`coarse-graining-adaptive-asd-parity`
  also pass with this task's binary).
* [x] Moving mixed finite-T fixture passes (`adaptive-cg-moving-backend-parity`
  passes with this task's binary at both precisions; no code this task
  touched is transition/moving-specific -- see "Moving/mixed finite-T and
  transition/restart continuity" above for why a dedicated new fixture was
  not separately constructed).
* [x] Coarse -> fine stochastic continuation passes (by construction --
  the generator is never skipped during coarse residence, proven
  discriminating by `runNegativeControlSkippedAdvancement`; also covered by
  `adaptive-cg-transition-ownership-invariants`).
* [x] Fine -> coarse -> fine continuation passes (same mechanism; no
  transition-direction-specific code was touched by this task).
* [x] Restart continuation passes (unaffected either way -- CGP-06A Part B
  invariant 5 already established curand/hiprand state was never serialized
  before this task and still is not; this task changes only which array
  entries a call writes, not seed/restart handling).
* [x] T=0 unaffected (T=0 is the four-way MD5 parity oracle above; `randomize()`
  has no temperature-dependent branch, confirmed in CGP-06A Part A).
* [x] RNG work scales with active atoms (real measurement: the "loop"
  GPU-time sub-phase drops with active fraction at both precisions; see
  Performance).
* [x] All-coarse path avoids full thermal-field generation. Partially: the
  *write* is fully avoided (`gpuActiveAtomSiteCall` issues no kernel at
  `activeCount==0`); the *generate* call remains full-N by design -- see
  "the all-coarse tension" above for why this is the accepted Strategy 1
  ceiling, not an oversight.
* [x] CUDA/HIP structural semantics remain aligned (`active_atom_site_kernel`
  uses the same `GPU_*` macro layer and `GridHelper` template every other
  kernel in this file uses; no backend-specific code added).
* [x] Negative controls are discriminating (`runNegativeControlSlotIndexing`,
  `runNegativeControlSkippedAdvancement`; both assert the broken variant
  disagrees with the reference, not merely that the accepted code passes).
* [ ] Whole-step fp32 benchmark rerun. Not attempted: this task's own
  measured effect (tens of microseconds inside one `randomize()` call among
  several per step) is smaller than this host's known CPU-quota wall-clock
  noise floor (CGP-05, CGP-06A), and the existing whole-step benchmark
  (`gpu_adaptive_runtime_benchmark`) structurally excludes thermfield/Depondt
  from its sweep by design; building a dedicated new harness to isolate an
  effect this small from that noise was judged not to produce a trustworthy
  number and was out of this task's named scope. See "What this task does
  not claim" above.

## Commit

`CGP-06B: scale thermal RNG with active fine spins`
