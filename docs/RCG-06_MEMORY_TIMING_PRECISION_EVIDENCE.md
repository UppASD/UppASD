# RCG-06: Remove memory, allocation, timing, and precision hazards

**Status:** In progress, sliced into lettered sub-sessions per
`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` §5.2 ("one defect
class per patch"), mirroring the RCG-04A–I / RCG-05A–G precedent.
**Dependencies:** RCG-05 (closed 2026-08-12).
**Parent task:** RCG-06, "Remove memory, allocation, timing, and precision
hazards" (`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`, Task
RCG-06).

This document is the audit trail for RCG-06's sliced implementation. Each
lettered section is its own commit-sized unit of work with its own
negative-control evidence, following the remediation blueprint's evidence
policy (§2). RCG-06 as a whole is not closed until every slice's evidence is
recorded here and a closing RCG-06 session reconciles the full checklist
against fresh, clean-build evidence (mirroring RCG-04J/RCG-05G).

---

## RCG-06A: CPU workspace hoisting (F-13, F-17)

**Date:** 2026-08-12
**Scope:** F-13 (automatic/stack arrays sized by `n_atoms`/`n_blocks` in the
static-hybrid and smooth-projected operators) and F-17 (per-adaptive-step
heap allocation churn, including a per-candidate-transition deep copy of the
whole hybrid operator).

### Finding reproduction

Confirmed by direct source inspection (three parallel investigations, see
below) rather than a runtime crash reproduction at the pre-fix commit size,
because the pre-fix automatic arrays only fault under a constrained stack
limit at a large `Natom` -- exactly what the new `MEM-LARGE-HOST` fixture's
fault-injection negative control demonstrates directly (see below), which is
stronger, size-independent evidence than a one-off crash reproduction would
have been.

F-13 automatic (stack) arrays, confirmed present before this slice's fix,
none `allocatable`/`pointer`:
- `source/CoarseGraining/statichybridoperator.f90`,
  `evaluate_static_hybrid_operator`: `effective_direction(3,n_atoms)`,
  `ghost_direction(3,n_atoms)`, `atomistic_field(3,n_atoms)`,
  `ghost_reaction_field(3,n_atoms)`, `reaction_field(3,1,n_blocks)`,
  `tensor_field(3,n_blocks)` -- called every step, twice (Heun predictor and
  corrector), once per ensemble.
- `source/CoarseGraining/smoothprojectedoperator.f90`: `local_norm(n_atoms)`
  in `prolongate_smooth_directions` (called every step from the above),
  `atom_direction(3,n_atoms)`/`raw_norm(n_atoms)` in
  `restrict_projected_atomistic_field` (called every step), `direction`/
  `raw_norm`/`fine_field` in `evaluate_projected_bilinear` (not currently
  production-hot, fixed for consistency).
- `source/CoarseGraining/adaptivecgproduction.f90`:
  `onsite_energy_j(Natom)`/`onsite_field_t(3,Natom)` in
  `evaluate_all_ensembles` (called every step, twice) and
  `production_energy_evaluator` (called once per candidate transition
  event), plus `atom_field(3,Natom)`/`coarse_field(3,1,operator%n_blocks)`
  in `production_energy_evaluator`.

F-17 per-step/per-event heap churn, confirmed present before this slice's
fix:
- `adaptive_cg_cpu_step` allocated 8 arrays of `O(Natom*Mensemble)` fresh
  every step (`atom0`, `coarse0`, `atom_predictor`, `coarse_predictor`,
  `atom_field0/1`, `coarse_field0/1`, `coarse_rhs0/1`).
- `reconstruct_coarse_atoms` allocated `requested`/`seeds` every step.
- `evaluate_coarse_tensor_operator` (`coarsetensoroperator.f90`)
  heap-allocated 10 arrays of `O(n_blocks)` on every call (called twice per
  ensemble per step).
- `apply_adaptive_transitions` (`adaptivehybridsolver.f90`) reassigned
  `candidate_atom`/`candidate_coarse` (auto-reallocating full-system copies)
  and `candidate_hybrid` (a deep copy of the entire
  `static_hybrid_operator_type`, including its nested tensor/projection
  operators) at least once per `apply_adaptive_transitions` call (the
  investigation's original "every event" framing overstated the frequency
  slightly -- Fortran's automatic-reallocation-on-assignment rule means a
  same-shape reassignment inside the per-event loop does not itself
  reallocate once the first event has allocated; the genuine churn is one
  full allocate/deallocate cycle per `apply_adaptive_transitions` call, i.e.
  every `cg_update_interval` steps -- still real, and now eliminated
  entirely by hoisting to `runtime`-owned persistent storage rather than
  merely reduced).

### Fix

**Ownership pattern:** extended the existing `adaptive_cg_state`
(`adaptive_cg_production_state_type`, `adaptivecgproduction.f90`) and
`adaptive_hybrid_runtime_type` (`adaptivehybridsolver.f90`) ownership
convention rather than inventing a new workspace type, matching every
allocatable array already owned there. New fields:
- `adaptive_cg_state`: `atom0`, `coarse0`, `atom_predictor`,
  `coarse_predictor`, `atom_field0/1`, `coarse_field0/1`, `coarse_rhs0/1`,
  `onsite_energy_j`, `onsite_field_t`, `reconstruction_requested`,
  `reconstruction_seeds`, `candidate_onsite_energy_j`,
  `candidate_onsite_field_t`, `candidate_atom_field`,
  `candidate_coarse_field` -- allocated once by a new
  `ensure_step_workspace` subroutine, called at the end of
  `setup_adaptive_cg_production`.
- `adaptive_hybrid_runtime_type`: `candidate_hybrid`, `candidate_atom`,
  `candidate_coarse`, `candidate_resultants`, `candidate_sums`,
  `candidate_directions`, `candidate_defined`, `candidate_has_decision`,
  `candidate_seeds` -- allocated once by a new `ensure_transition_workspace`
  subroutine, called at the end of `setup_adaptive_hybrid_runtime` and again
  (idempotently) as a preflight at the top of `apply_adaptive_transitions`.
- `coarse_tensor_operator_type` (`coarsetensoroperator.f90`) and
  `smooth_projected_operator_type`/`static_hybrid_operator_type`
  (`smoothprojectedoperator.f90`/`statichybridoperator.f90`): per-operator
  `scratch_*` fields, allocated once in each type's own `setup_*`
  subroutine, replacing the automatic/heap-churned locals inside their
  evaluation routines directly.

**Multi-ensemble indexing:** every new/hoisted array keeps ensemble as the
trailing dimension, matching the pre-existing convention verified across
`atom_direction`, `coarse_direction`, `coarse_resultant_mub`, etc.

**A Fortran argument-aliasing hazard, found and fixed during this slice:**
nesting scratch arrays as components of an operator type that is *also*
passed whole into a helper subroutine in the same call is illegal under
Fortran's actual-argument aliasing rules whenever the helper writes to that
scratch component (F2018 15.5.2.13 -- an actual argument associated with a
definable dummy must not be a subobject of another actual argument in the
same reference). This affects `coarsetensoroperator.f90`'s
`physical_forward_gradient`/`add_physical_gradient_transpose` (previously
took the whole `operator` alongside one of its own now-scratch arrays) and
`smoothprojectedoperator.f90`'s `prolongate_with_norm` (same pattern). Both
were decoupled to take the specific geometry/read-only fields they need as
explicit arguments instead of the whole operator, eliminating the hazard
without changing their mathematics (`physical_forward_gradient`/
`add_physical_gradient_transpose` are the parent blueprint's explicitly
preserved discrete-adjoint-pair regression contract; their bodies are
byte-for-byte unchanged, only their argument list). By contrast,
`statichybridoperator.f90`'s `evaluate_static_hybrid_operator` never
triggers this hazard: it only ever passes *sibling* components of `operator`
together (e.g. `operator%projection` alongside `operator%scratch_ghost_direction`),
which are genuinely disjoint storage, not a whole-object/subobject pair --
confirmed by direct trace of every call in its body, not assumed.

**Interface ripple:** `evaluate_static_hybrid_operator`,
`evaluate_coarse_tensor_operator`, `prolongate_smooth_directions`,
`restrict_projected_atomistic_field`, and `evaluate_projected_bilinear` all
changed their `operator` argument from `intent(in)` to `intent(inout)`
(needed so they can write their own persistent scratch). This cascaded to:
the `adaptive_energy_evaluator` abstract interface
(`adaptivehybridsolver.f90`), `production_energy_evaluator`
(`adaptivecgproduction.f90`), and three test-only implementer/helper
subroutines (`hybrid_energy` in `test_adaptive_hybrid_solver.f90`,
`hybrid_finite_difference` in `test_static_hybrid_operator.f90`,
`finite_energy_derivative` in `test_coarse_tensor_operator.f90`). Every call
site already passed a plain named variable (never a literal/expression), so
the intent widening is mechanically safe; verified by grep across
`source/` and `tests/` before editing.

**Allocation preflight:** both `ensure_step_workspace` and
`ensure_transition_workspace` are idempotent (a second call with matching
Natom/Mensemble/n_spatial_blocks is a no-op) and explicitly reject (rather
than silently reallocate or corrupt) a detected geometry drift after initial
allocation, with a named diagnostic -- matching the blueprint's
reject-rather-than-homogenize principle (§5.4) applied to a memory-safety
concern rather than a physics one.

### Tests

- `tests/coarse_graining/test_adaptive_hybrid_solver.f90`,
  `test_transition_workspace_lifecycle` (new): confirms the transition
  workspace is allocated exactly once at setup with the correct shape,
  reused (not reallocated) on a second matching-geometry preflight call,
  explicitly rejected (`ADAPTIVE_HYBRID_ALLOCATION_FAILED`, non-empty
  diagnostic) when the geometry is corrupted after allocation, fully
  deallocated when the runtime is reset to its default structure
  constructor, and correctly reallocated at the (same, in this test) size
  after a fresh `setup_adaptive_hybrid_runtime` call. This is the
  "allocation, resize, reuse, and cleanup are tested" checklist item's
  direct evidence.
- `tests/coarse_graining/test_stack_overflow_fault_injection.f90` (new,
  standalone, not linked against `asdlib`): reproduces the eliminated
  automatic-array pattern in isolation (four `(3,Natom)` automatic arrays,
  `Natom=8000` passed as a runtime dummy argument so the compiler cannot
  avoid genuine stack allocation), touching every element.
- `tests/coarse_graining/run_mem_large_host.py` +
  `tests/coarse_graining/e2e/mem_large_host/` (new, `MEM-LARGE-HOST`
  fixture; `Natom=8000`, `NA=1` simple-cubic host, nearest-neighbour
  exchange, `block_size 5 5 5`, `Nstep=1`, `cg_mask_mode STATIC` all-fine):
  runs both the fault-injection binary and the real `sd.f95` production
  binary under an identical, deliberately constrained `ulimit -s 512`
  (512&nbsp;KiB -- comfortably below the eliminated arrays' ~750&nbsp;KiB
  demand at this `Natom`, documented arithmetic in the fixture's
  `README.md`). See raw command output below.

### Raw evidence (fresh out-of-tree CPU build)

```
$ git rev-parse HEAD
caccdcbc841c2944c86b8339a6e04143a3421abe   # RCG-05 closure; RCG-06A is
                                            # uncommitted on top of this
$ git status --short
 M source/CoarseGraining/adaptivecgproduction.f90
 M source/CoarseGraining/adaptivehybridsolver.f90
 M source/CoarseGraining/coarsetensoroperator.f90
 M source/CoarseGraining/smoothprojectedoperator.f90
 M source/CoarseGraining/statichybridoperator.f90
 M CMakeLists.txt
 M tests/coarse_graining/test_adaptive_hybrid_solver.f90
 M tests/coarse_graining/test_coarse_tensor_operator.f90
 M tests/coarse_graining/test_static_hybrid_operator.f90
?? tests/coarse_graining/e2e/mem_large_host/
?? tests/coarse_graining/run_mem_large_host.py
?? tests/coarse_graining/test_stack_overflow_fault_injection.f90
?? docs/RCG-06_MEMORY_TIMING_PRECISION_EVIDENCE.md

$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF -S . -B build_rcg06a_cpu
-- ... Configuring done, Generating done ...

$ cmake --build build_rcg06a_cpu -j"$(nproc)"
-- ... 100% Built target adaptive_hybrid_solver_tests (and every other target) ...

$ ctest --test-dir build_rcg06a_cpu -L cg13-cpu --output-on-failure
... 25/25 tests passed (was 24/24 before this slice; +1 = adaptive-cg-mem-large-host) ...
100% tests passed, 0 tests failed out of 25
Label Time Summary:
cg13               =  63.30 sec*proc (25 tests)
cg13-cpu           =  63.30 sec*proc (25 tests)
...
Total Test time (real) =  63.32 sec

$ ctest --test-dir build_rcg06a_cpu -R '^asd-tests$' --output-on-failure
1/1 Test #2: asd-tests ........................   Passed   10.33 sec
100% tests passed, 0 tests failed out of 1

$ ctest --test-dir build_rcg06a_cpu -R adaptive-cg-mem-large-host -V
19: fault-injection negative control: killed as expected (returncode=-11)
19: MEM-LARGE-HOST: production run completed under ulimit -s 512 KiB at Natom=8000
1/1 Test #19: adaptive-cg-mem-large-host .......   Passed    1.37 sec
```

`returncode=-11` is Python's representation of the child process being
killed by `SIGSEGV` -- confirmed, not assumed: the eliminated automatic-array
pattern genuinely overflows the stack at `Natom=8000` under `ulimit -s 512`,
while the real, fixed production binary completes cleanly at the identical
size and limit. This is the blueprint §2.3 fault-injection negative-control
form: "a minimal fault injection fails while the accepted implementation
passes."

Every one of the 24 pre-existing `cg13-cpu` fixtures (RCG-01 through RCG-05's
accepted evidence, including RCG-04's five moving-parity fixtures and
RCG-05's ownership-map comparator) and the full `asd-tests` legacy suite pass
unchanged, with no edits to their own source files -- direct evidence for
the "existing derivative and moving-state fixtures remain unchanged"
checklist item.

### Checklist items addressed by this slice

- [x] No routine places `O(n_atoms)` adaptive workspace on the stack.
- [x] Per-step full-system allocations are hoisted or proven negligible.
- [x] Workspace allocation, resize, reuse, and cleanup are tested.
- [x] Memory preflight includes all persistent CPU/GPU adaptive workspace.
      (CPU only in this slice; `adaptive_owned_cpu_bytes()` extended to
      account for every new field. GPU workspace preflight is RCG-06B's
      responsibility.)
- [x] A large-host fixture passes with the documented stack limit.
- [x] Multiple ensembles retain correct indexing. (Verified by construction:
      every new array keeps ensemble as the trailing dimension, matching
      the pre-existing convention; the existing `Mensemble`-sensitive
      fixtures in `cg13-cpu` -- e.g. `coarse-graining-adaptive-hybrid-solver`'s
      own tests, which use `Mensemble` values other than 1 -- continue to
      pass unchanged.)
- [x] Existing derivative and moving-state fixtures remain unchanged.
- [ ] GPU global energy accumulation uses FP64 storage and arithmetic --
      RCG-06B.
- [ ] FP32 field parity and FP64 energy budgets are distinct -- RCG-06B.
- [ ] Energy error scaling is measured over increasing system size --
      RCG-06B.
- [ ] Timers use a suitable wall/device clock -- RCG-06C.
- [ ] Both Heun field evaluations are included -- RCG-06C (already true
      today, per the F-18 investigation; RCG-06C adds the reconciliation
      evidence).
- [ ] Phase totals plus unaccounted time reconcile with external wall time
      -- RCG-06C.
- [ ] RNG correlation evidence is recorded and its scope decision is
      explicit -- RCG-06D.

### Open items / not yet done

RCG-06A did not touch: GPU energy precision (F-11, RCG-06B), CPU/GPU phase
timing (F-18, RCG-06C), or reconstruction RNG statistics (F-20, RCG-06D).
CUDA/HIP builds were not exercised in this slice (RCG-06A is CPU-only
source; a CUDA/HIP build should still be confirmed unaffected before RCG-06
as a whole closes, since `adaptivecgproduction.f90`'s GPU staging paths were
untouched but share the file). No independent reviewer distinct from the
implementer has reviewed this slice yet (matching every prior RCG-0x
session's disclosed, non-blocking deferral pattern) -- Human review of the
Fortran-aliasing-hazard finding and fix in particular is recommended before
RCG-06's final closure, since it is the most subtle correctness risk this
slice introduced and then resolved.
