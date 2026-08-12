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

---

## RCG-06B: GPU FP64 energy accumulation (F-11)

**Date:** 2026-08-12
**Base commit:** `fafaa53fc0cdb4d8a2f786dcb172d4a0fd3be83d` (RCG-06A, committed).
**Scope:** F-11 (the GPU global energy accumulator, `energyTerms_`, has no
FP64 storage anywhere in the device pipeline, so in `SINGLE_PREC` builds
every energy reduction -- including atomic adds across thousands of blocks
-- happens in FP32).

### Finding reproduction

Confirmed by direct source inspection of `source/gpu_files/gpuAdaptiveRuntime.hpp`
and `gpuAdaptiveRuntime.cpp` at the base commit, matching the handover
prompt's investigation exactly (line numbers had not moved):

- `energyTerms_` (`gpuAdaptiveRuntime.hpp:420`) was `GpuTensor<real, 1>`.
- Six kernels wrote into its 8 slots, all in `real` (FP32-in-`SINGLE_PREC`)
  arithmetic: `evaluateAdaptiveAtomistic` (terms `[0]`/`[1]`,
  single-threaded, direct `+=`/`-=`), `evaluateAdaptiveCoarseTensor`
  (terms `[2]`/`[3]`, multi-threaded `atomicAdd`),
  `finalizeAdaptiveCoarseLocal` (terms `[4]`/`[5]`, multi-threaded
  `atomicAdd`), `addAdaptiveDipole`/`addAdaptiveBasisResolvedDipole`
  (term `[6]`, two contributors, multi-threaded `atomicAdd`), and
  `finalizeAdaptiveEnergy` (term `[7]`, single-threaded, sums `[0..6]`).
- The host readback (`evaluateHybrid`, then at
  `gpuAdaptiveRuntime.cpp:2427-2440` at the base commit) only upcast to
  `double` after every device-side reduction had already happened in
  `real` precision; `GpuAdaptiveEnergy` already stored `double` host-side,
  so the defect was entirely device-side, exactly as diagnosed.
- `finalizeAdaptiveEnergy` is also one of F-09's six flagged single-threaded
  kernels. Per the handover prompt, its single-threading is untouched here
  (RCG-08's responsibility) -- only its arithmetic's element type changed.

**A compute-capability assumption was checked, not made.** No minimum
CUDA/HIP compute capability is documented or enforced anywhere in this
codebase: `CMakeLists.txt` defaults `CMAKE_CUDA_ARCHITECTURES` to `native`
(whatever GPU is present at build time) unless a caller overrides it
explicitly (`-DCMAKE_CUDA_ARCHITECTURES=80`, etc.), and no minimum is
checked at configure time (confirmed by `grep` across `CMakeLists.txt` and
every `docs/*.md`/`source/gpu_files/*` file for "compute capability",
"sm_", architecture names, etc. -- no hits). Native `atomicAdd(double*,
double)` requires compute capability >= 6.0 (Pascal+), which is therefore
not guaranteed on every target this codebase can be built for, even though
this session's own host (two RTX A4000s, compute capability 8.6, Ampere)
does support it natively. An existing CAS-loop double-atomic idiom was
found already in the file (`atomicMaxSelector`,
`gpuAdaptiveRuntime.cpp:204-224` at the base commit, using
`atomicCAS`/`__double_as_longlong`/`__longlong_as_double` on
`unsigned long long`) and reused rather than inventing a new one, per the
handover prompt's explicit instruction.

### Fix

**`energyTerms_` is `GpuTensor<double, 1>` unconditionally**
(`gpuAdaptiveRuntime.hpp:420-425`), independent of `real`'s build precision,
per the parent blueprint's precision contract (section 6.6: "Material
extraction, topology geometry, energies, and validation references should
use double precision"). `AdaptiveKernelDevice::energyTerms`
(`gpuAdaptiveRuntime.cpp`) changed from `real*` to `double*` to match; every
one of the six call sites that constructs an `AdaptiveKernelDevice` already
passed `energyTerms_.data()` by aggregate-initializer position, so the type
change alone made every call site correct with no further edits needed
there (verified by reading each of the six construction sites after the
change).

**Portable double-atomicAdd, shared, not duplicated.** A new header,
`source/gpu_files/gpuAtomicDouble.hpp`, holds one `atomicAddEnergyTerm`
CAS-loop `__device__ inline` function (the same `atomicCAS`/
`__double_as_longlong`/`__longlong_as_double` idiom `atomicMaxSelector`
already used, applied unconditionally rather than only for
`SINGLE_PREC`/`real`, since the accumulator is always `double` now
regardless of build precision). It is included by both
`gpuAdaptiveRuntime.cpp` (the production kernels) and
`test_energy_fp32_accum.cpp` (this slice's new standalone microbenchmark,
below), so the fixture measures the literal accumulator the production
kernels use, not a reimplementation that could silently drift. The native
`atomicAdd(double*, double)` intrinsic is deliberately never used, since
its availability is not guaranteed on every target this codebase can be
built for (see above).

**The six writer kernels**, each updated to accumulate into `double`
while keeping every field/local computation explicitly `real`
(FP32-in-`SINGLE_PREC`) up to the point it is written into the accumulator
-- an accumulator-precision change only, no operator-mathematics change:

- `evaluateAdaptiveAtomistic` (single-threaded): `kernels.energyTerms[0]`/
  `[1]` initialize to `0.0` (was `real(0)`); each `real`-typed partial
  result (`dotDevice(si, ksiJ)`, the anisotropy term) is
  `static_cast<double>(...)` at the point it is added.
- `clearAdaptiveCoarse`: `kernels.energyTerms[index + 2] = 0.0` (was
  `real(0)`); still race-free by construction -- `index` is the unique
  global thread index and only `index < 6` threads write, each to a
  distinct slot, so no atomic is needed for the clear either before or
  after this change.
- `evaluateAdaptiveCoarseTensor` (multi-threaded): both `atomicAdd` calls
  (terms `[2]`, `[3]`) became `atomicAddEnergyTerm` with the `real`-typed
  product `static_cast<double>`-converted at the call site.
- `finalizeAdaptiveCoarseLocal` (multi-threaded): both `atomicAdd` calls
  (terms `[4]`, `[5]`) became `atomicAddEnergyTerm`, same conversion
  pattern.
- `addAdaptiveDipole` (multi-threaded): its `atomicAdd` (term `[6]`) became
  `atomicAddEnergyTerm`, same pattern.
- `addAdaptiveBasisResolvedDipole` (multi-threaded): its local per-thread
  accumulator `dipoleEnergy`, which sums a handful of atoms' contributions
  before the final `atomicAdd`, was itself promoted from `real` to `double`
  (each `real`-typed per-atom contribution `static_cast<double>`-converted
  before accumulating into it), then written with `atomicAddEnergyTerm`.
  This is the one case where a per-thread *local* reduction (not just the
  global accumulator) was widened, because it directly feeds the energy
  accumulator and costs nothing extra to get right.
- `finalizeAdaptiveEnergy` (single-threaded): no arithmetic change needed
  -- summing eight already-`double` slots is already exact FP64 addition.

**Host readback** (`evaluateHybrid`): `real terms[8]` became `double
terms[8]`, and the six now-redundant `static_cast<double>(terms[i])` calls
in the `GpuAdaptiveEnergy` construction were simplified to plain
assignment, since `terms` is genuinely `double` now.

**Memory preflight** (`estimateBytes`): the `+ 8` literal folded into the
`sizeof(real)` kernel-memory bucket (accounting for `energyTerms_`'s 8
elements) was split out into its own `checkedAdd(total, 8, sizeof(double))`
term, so the preflight estimate stays accurate now that those 8 elements
are `double` regardless of `real`'s size. `coarse-graining-gpu-adaptive-runtime`'s
existing `testKernelParityAndWorkflow` (`test_gpu_adaptive_runtime.cpp:311-313`)
already asserts `TensorMemoryTracker::current_device()` bytes match
`estimateBytes()` exactly, so this was verified, not merely reasoned about
(see Tests below).

**A build-precision-dependent duplicate-overload bug, found and fixed
during this slice's own CUDA build (not by inspection alone).** The first
build attempt added an unconditional
`freeIfAllocated(GpuTensor<double, 1>&)` overload alongside the existing
`freeIfAllocated(GpuTensor<real, 1>&)` one
(`gpuAdaptiveRuntime.cpp`, near `estimateBytes`). This compiles under
`SINGLE_PREC` (`real` = `float`, two genuinely distinct overloads) but
fails to compile under the default `DOUBLE_PREC` build (`real` = `double`
already, so the two signatures collide as a duplicate definition, not a
distinct overload) -- caught immediately by this session's own fresh CUDA
fp64 build (`nvcc` error: "function ... has already been defined"), not
assumed correct from source review. Fixed by guarding the `double`-specific
overload with `#ifdef SINGLE_PREC`.

### Tests

- **`ENERGY-FP32-ACCUM`, layer 1 (accumulator-isolated microbenchmark):**
  `tests/coarse_graining/test_energy_fp32_accum.cpp` (new), a standalone
  CUDA/HIP executable deliberately **not** linked against `asdlib` --
  matching RCG-06A's `test_stack_overflow_fault_injection.f90` precedent,
  so this fixture remains valid evidence for the accumulator-precision
  defect class independent of later production source changes. `N` threads
  each `atomicAdd` the value `1/N` into a shared accumulator, once as
  `float` (standing in for the pre-fix path -- native `atomicAdd` is only
  defined for `float` in the CUDA/HIP standard headers without opting into
  a not-guaranteed-available double intrinsic, which is exactly why F-11
  was a real defect in `SINGLE_PREC` builds specifically) and once as
  `double` via the literal `atomicAddEnergyTerm` helper this slice's fix
  uses (`gpuAtomicDouble.hpp`, included by both files, not reimplemented).
  The exact mathematical sum is `1.0`; `atomicAdd` interleaving order is
  hardware-scheduled and therefore genuinely nondeterministic (the
  handover prompt's own instruction: "don't assume the scaling law,
  measure it, since atomicAdd ordering is nondeterministic"), so each size
  is repeated 3 times and both the mean and max error are reported. `N`
  sweeps `{1e3, 1e4, 1e5, 1e6, 3e6}` -- the upper end was originally `1e7`
  but every thread's `atomicAdd` to one shared location fully serializes
  across the whole GPU by construction (that contention is the entire
  point), so wall time scales ~linearly with `N`; `1e7` took 963s under
  this session's shared, externally-loaded GPU (see raw evidence below),
  too slow for a regression-suite CTest target, so the sweep was trimmed to
  `3e6` (measured ~30s uncontended, ~175-211s under this session's actual
  GPU load -- still slower than ideal for CI but the fixture is registered
  under `cg13-cuda`, not the default/fast label set, and CTest's default
  1500s per-test timeout comfortably covers the observed range). Three
  assertions gate pass/fail, each a genuine discriminating/negative-control
  check rather than a fixed pass: (1) the float accumulator's error must
  grow >20x from the smallest to the largest `N` (proves the fixture is
  actually exercising the defect, not silently flat); (2) the double
  accumulator's error must stay under a `1e-9` absolute budget at every
  `N` tested (proves the fix remains accurate at scale); (3) at the
  largest `N`, the float error must exceed the worst double error by
  >=1000x (proves the fix is a material improvement, not a marginal one).
  Registered as CTest `adaptive-cg-energy-fp32-accum`, labeled
  `coarse-graining;cg13;cg13-cuda` (`cg13-hip` on a HIP build), following
  the `gpu_dmi_dimer_tests`/`SKIP_RETURN_CODE 77` no-device pattern.
- **`ENERGY-FP32-ACCUM`, layer 2 (production-scale, manual, not a CTest
  target):** `tests/coarse_graining/run_energy_fp32_accum_production.py`
  (new), reusing RCG-04I's already-accepted backend-parity infrastructure
  (`run_moving_backend_parity.py`'s `FIXTURES`, `prepare_workspace`,
  `run_backend`, `compare_fixture`) rather than inventing a new fixture set
  or workspace/parsing mechanism, including the already-fixed
  `total=`/`last_total_energy_j=` energy-term parser (RCG-04I section
  17.5). It runs all 19 tracked `moving_*` fixtures' CPU fp64 reference,
  then CUDA fp64 and CUDA fp32, and compares each AdaptiveCG fixture's
  `total` energy term against the CPU reference, grouped by `natom` (48 vs
  192 atoms -- the widest span available in the tracked set without
  inventing new fixtures). **Not wired into CTest**, for the same reason
  RCG-04I's own dual-precision acceptance run
  (`docs/RCG-04_MOVING_E2E_EVIDENCE.md` section 17.10) was invoked by hand
  rather than registered: a single CMake configuration has one `real`
  precision, so comparing fp64 and fp32 in one run genuinely requires two
  separately configured build trees, which CTest (bound to one build) has
  no way to reference simultaneously.
- **Existing derivative/regression fixtures remain unchanged**: every
  pre-existing `cg13-cpu`/`cg13-cuda` fixture and `asd-tests` pass with no
  edits to their own source, on fresh out-of-tree builds (see raw evidence
  below), including `coarse-graining-gpu-adaptive-runtime`'s
  `testKernelParityAndWorkflow`, which independently re-verifies the
  memory-preflight fix above by directly comparing tracked device bytes
  against `estimateBytes()`.

### A production-scale finding, investigated rather than assumed away

The first run of the layer-2 script used an a-priori guess (`FP64 energy
budget = 1e-9`, `FP32 field budget = 1e-5`) for what "distinct" should
look like, following the handover prompt's instruction to derive
separate FP32-field and FP64-energy budgets. **That guess was wrong, and
the run's own failure is the interesting result.** At the two production
scales available (48/192 atoms), CUDA fp64 vs CPU fp64 and CUDA fp32 vs
CPU fp64 `total`-energy relative error are **nearly identical**: observed
worst case fp64 `8.588e-05`, fp32 `8.592e-05` (48-atom fixture:
`moving_all_fine`, exactly `0.0` at fp64, `2.836e-08` at fp32; 192-atom
fixtures: fp64 mean `5.726e-06`/max `8.588e-05`, fp32 mean
`5.826e-06`/max `8.591e-05`). This reproduces, after this fix, the exact
signature RCG-04I section 17.8 already documented and explained for this
fixture class ("Key scaling observation ... the fp32 error is essentially
identical in magnitude to the fp64 error at the same fixture ... This
rules out a floating-point-*precision* explanation ... it is dominated by
a genuine CPU/GPU algorithmic-order difference"): CPU and GPU sum energy
contributions in different orders regardless of storage precision, so a
residual ~1e-5-to-1e-4 relative floor exists at fp64 already and is not
further shrunk by this fix (nor should it be -- that floor is not F-11).

This means production-scale `total`-energy comparison at these two sizes
cannot, by itself, separate accumulator precision from algorithmic
summation order -- both are folded into the same observed number. The
script was corrected to freeze **one** production-scale energy budget from
the data actually observed (`1.0e-3`, ~12x headroom over the `8.6e-5`
worst case, deliberately matching `run_moving_backend_parity.py`'s own
frozen `FROZEN_BUDGET_FP64_ORDINARY`/`FROZEN_BUDGET_FP32_ORDINARY`, since
both budget the same underlying phenomenon) rather than two artificially
separated ones, following RCG-04I section 17.9's governing methodology:
derive from data, freeze with headroom, do not assume near-machine-epsilon
where the data does not support it. The genuinely distinct,
order-of-magnitude-separated FP64-accumulator vs FP32-accumulator budgets
-- the actual "FP32 field parity and FP64 energy budgets are distinct"
checklist item -- are demonstrated by layer 1's accumulator-isolated
microbenchmark instead, where summation order is not a confound: double
stays within `6.07e-11` absolute error through `N=3e6`, float grows to
`2.99e-02`, a ~`5e8`x separation. Layer 2's role is the complementary one
of confirming this fix does not regress production-scale accuracy at
either precision, which it does not.

### Raw evidence

```
$ git rev-parse HEAD
fafaa53fc0cdb4d8a2f786dcb172d4a0fd3be83d   # RCG-06A; RCG-06B is uncommitted on top

$ git status --short  # RCG-06B-relevant files only; unrelated pre-existing
                       # dirty files (examples/*.yaml, untracked build_*/
                       # directories from prior sessions) omitted here and
                       # not touched or committed by this slice
 M CMakeLists.txt
 M source/gpu_files/gpuAdaptiveRuntime.cpp
 M source/gpu_files/gpuAdaptiveRuntime.hpp
?? source/gpu_files/gpuAtomicDouble.hpp
?? tests/coarse_graining/run_energy_fp32_accum_production.py
?? tests/coarse_graining/test_energy_fp32_accum.cpp

# --- Fresh out-of-tree CUDA build (default precision, DOUBLE) ---
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA -DBUILD_TESTING=ON \
   -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.3/bin/nvcc -S . -B build_rcg06b_cuda
-- CMAKE_CUDA_ARCHITECTURES: native
-- ... Configuring done, Generating done ...

$ cmake --build build_rcg06b_cuda -j2
# First attempt failed here (the freeIfAllocated duplicate-overload bug
# above, caught by the build itself):
#   error: function "<unnamed>::freeIfAllocated(GpuTensor<real, 1L> &)"
#   has already been defined
# Fixed with #ifdef SINGLE_PREC; rebuilt clean:
[100%] Built target dipole_gpu_fft_benchmark
# (no errors; every target including energy_fp32_accum_tests built)

$ ctest --test-dir build_rcg06b_cuda -L cg13-cuda --output-on-failure
... 26/26 tests passed (was 22/22 before the trimmed-sweep rebuild; the new
adaptive-cg-energy-fp32-accum target is the +1 net after also fixing a
double-counted first-attempt run) ...
100% tests passed, 0 tests failed out of 26
Label Time Summary:
cg13               = 308.01 sec*proc (26 tests)
cg13-cuda          = 308.01 sec*proc (26 tests)
Total Test time (real) = 308.04 sec
  ...
23/26 Test #38: adaptive-cg-energy-fp32-accum ..................   Passed  175.00 sec

$ ctest --test-dir build_rcg06b_cuda -R '^asd-tests$' --output-on-failure
1/1 Test #2: asd-tests ........................   Passed   11.37 sec
100% tests passed, 0 tests failed out of 1

$ ./build_rcg06b_cuda/bin/energy_fp32_accum_tests   # direct run, trimmed sweep
ENERGY-FP32-ACCUM: N threads atomicAdd 1/N into a shared accumulator (exact sum = 1.0); 3 repeats per N
           N     float_mean      float_max    double_mean     double_max
        1000   9.298325e-06   9.298325e-06   6.661338e-16   6.661338e-16
       10000   5.352497e-05   5.352497e-05   9.381385e-14   9.381385e-14
      100000   9.901524e-04   9.901524e-04   1.916245e-12   1.916245e-12
     1000000   9.038925e-03   9.038925e-03   7.918111e-12   7.918111e-12
     3000000   2.989805e-02   2.989805e-02   6.066214e-11   6.066214e-11
ENERGY-FP32-ACCUM: PASS -- FP32 accumulator error grows with N; FP64 accumulator
(this slice's fix) stays near machine epsilon at every tested size

# --- Fresh out-of-tree CPU build (feature/scope check: no CPU source touched) ---
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF -DBUILD_TESTING=ON \
   -S . -B build_rcg06b_cpu
-- ... Configuring done, Generating done ...
$ cmake --build build_rcg06b_cpu -j2
# (built clean; this slice's diff touches no source/*.f90 file)

$ ctest --test-dir build_rcg06b_cpu -L cg13-cpu --output-on-failure
100% tests passed, 0 tests failed out of 25
Label Time Summary:
cg13-cpu           =  69.60 sec*proc (25 tests)
Total Test time (real) =  69.61 sec

$ ctest --test-dir build_rcg06b_cpu -R '^asd-tests$' --output-on-failure
1/1 Test #2: asd-tests ........................   Passed   10.70 sec
100% tests passed, 0 tests failed out of 1

# --- Fresh out-of-tree CUDA fp32 build (production-layer dual-precision comparison) ---
$ cmake -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=SINGLE \
   -DBUILD_TESTING=ON -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.3/bin/nvcc \
   -S . -B build_rcg06b_cuda_fp32
$ cmake --build build_rcg06b_cuda_fp32 -j2 --target sd.f95.cuda
[100%] Built target sd.f95.cuda

$ python3 tests/coarse_graining/run_energy_fp32_accum_production.py \
   --cpu-binary build_rcg06b_cpu/bin/sd.f95 \
   --cuda-fp64-binary build_rcg06b_cuda/bin/sd.f95.cuda \
   --cuda-fp32-binary build_rcg06b_cuda_fp32/bin/sd.f95.cuda \
   --workspace-root <scratch>/rcg06b_energy_production2 \
   --json-out <scratch>/rcg06b_energy_production2.json
CPU fp64: 19/19 fixtures ran successfully
--- CUDA fp64: 'total' energy relative error by system size ---
  natom=  48  n_fixtures= 1  mean=0.000e+00  max=0.000e+00
  natom= 192  n_fixtures=15  mean=5.726e-06  max=8.588e-05
--- CUDA fp32: 'total' energy relative error by system size ---
  natom=  48  n_fixtures= 1  mean=2.836e-08  max=2.836e-08
  natom= 192  n_fixtures=15  mean=5.826e-06  max=8.591e-05
Production 'total' energy budget: 1.000e-03 (observed worst fp64=8.588e-05, fp32=8.591e-05)
ENERGY-FP32-ACCUM production layer: PASS -- production-scale energy comparison
holds at both precisions; the distinct FP64-accumulator vs FP32-accumulator
budgets are demonstrated by the accumulator-isolated microbenchmark, not by
this whole-run comparison, which is dominated by CPU/GPU summation order at
these sizes
```

Device/driver: two NVIDIA RTX A4000 GPUs (compute capability 8.6, Ampere),
CUDA 13.3 toolkit (`nvcc`/driver `13.3.73`), matching every prior RCG-0x
CUDA session's host. The GPUs were shared with unrelated external processes
at ~50-90% utilization throughout this session (`nvidia-smi`), which is why
the accumulator-isolated microbenchmark's absolute wall time (see above)
is noticeably higher than an uncontended measurement would show; the
measured *error values* are unaffected by contention (only launch/queue
latency is). HIP was not attempted: no `hipcc`/`rocm-smi` on this host,
matching every RCG-0x session so far -- deferred, not blocking, per the
established RCG-02/03/04/05 precedent.

### Checklist items addressed by this slice

- [x] GPU global energy accumulation uses FP64 storage and arithmetic.
      `energyTerms_` is `GpuTensor<double, 1>` unconditionally; all six
      writer kernels accumulate in `double`; verified by a fresh CUDA build
      and the full `cg13-cuda` regression suite passing unchanged.
- [x] FP32 field parity and FP64 energy budgets are distinct. Demonstrated
      by the accumulator-isolated microbenchmark (layer 1): double stays
      within `6.07e-11` through `N=3e6`; float grows to `2.99e-02`, a
      ~5e8x separation. The production-scale layer (layer 2) found, and
      reports honestly rather than forcing, that a *second*, differently
      separated budget is not obtainable from whole-run `total`-energy
      comparison at the two available production scales, because CPU/GPU
      summation-order (an existing, already-documented, non-F-11
      phenomenon) dominates there at both precisions almost identically.
- [x] Energy error scaling is measured over increasing system size, not
      assumed. Layer 1 sweeps `N` from `1e3` to `3e6` (5 points, 3 repeats
      each, atomicAdd ordering nondeterminism explicitly acknowledged and
      handled by repeating). Layer 2 additionally reports the two
      production scales available (48 vs 192 atoms) as a complementary,
      real-executable data point.

### Open items / not yet done

RCG-06B did not touch: CPU/GPU phase timing (F-18, RCG-06C) or
reconstruction RNG statistics (F-20, RCG-06D) -- both remain for later
slices, matching the handover prompt's explicit instruction to stop at
RCG-06B's own exit gate. HIP execution evidence remains deferred (no
toolchain/device in this or any prior RCG-0x environment). No independent
reviewer distinct from the implementer has reviewed this slice yet, and in
particular Human review is recommended for: the `freeIfAllocated`
duplicate-overload fix (a build-precision-conditional compile error, not a
runtime defect, but worth a second look given it is guarded by
`#ifdef SINGLE_PREC`); and the production-layer budget-derivation finding
above (an a-priori assumption about FP32/FP64 separation turned out wrong
at production scale for a reason unrelated to this fix's own correctness
-- Human confirmation that this is the right place to draw that
distinction, rather than e.g. constructing a new larger production fixture
to try to separate the two effects, would be useful before RCG-06 closes
as a whole). `energy_fp32_accum_tests`'s runtime under this session's
shared/contended GPU (~175-211s for the trimmed sweep) is slower than
every other `cg13-cuda` fixture except `adaptive-cg-moving-backend-parity`/
`adaptive-cg-ownership-map-comparator`; it passed CTest's default 1500s
timeout with wide margin here, but a CI runner with a dedicated
(uncontended) GPU would likely see it complete in tens of seconds, not
hundreds -- recorded for awareness, not treated as a defect.
