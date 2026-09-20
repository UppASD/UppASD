# CGP-00B — Transition-energy gating evidence

Branch: `gpu_hip_cu_ab_cg`. See `docs/CGP_work.md` CGP-00B for the task
definition this document reports against.

## 1. State machines as audited (before any code change)

### CPU (`source/CoarseGraining/adaptivehybridsolver.f90`,
`apply_adaptive_transitions`)

```
advance_selector_state (selector/dwell/hard-mask proposal)
        |
  no decisions%event allocated -> return (selector-only; unchanged)
        |
  [OLD] energy_evaluator(current hybrid) -> energy_before   <- unconditional
        |
  for each proposed event:
     restrict/reconstruct candidate
     rebuild_static_hybrid_ownership
     [OLD] energy_evaluator(candidate hybrid) -> energy_after   <- unconditional
     accepted = |energy_after-energy_before| <= energy_jump_limit_j
     accept -> publish candidate; energy_before = energy_after
     reject -> roll back resolution_state/age/epoch
```

`energy_jump_limit_j` defaulted to `huge(1.0_dblprec)`, which is a *finite*
number, so the existing `valid_reconstruction_configuration` finite/
nonnegative check passed trivially and the energy evaluator ran on every
proposed transition regardless of whether a user had ever configured a
threshold. This is the "implicit sentinel" contract CGP-00B replaces.

### GPU (`source/gpu_files/gpuAdaptiveRuntime.cpp`,
`proposeSelectorState`/`publishProposedState`, driven from
`gpuSimulation.cpp::advanceAdaptiveStep`)

```
restrictMoments / evaluateSelectorScores / evaluatePolarizationGate
        |
proposeSelectorState(policy, polarizationUnsafeMask)
        |
publishProposedState(direction, reconstructionPolicy, completeStep)
   -- called with acceptedBlockMask == nullptr (its default) at the one
      production call site (gpuSimulation.cpp:1220) --
        |
  every proposed transition is published; no energy is evaluated anywhere
  in this path
```

`adaptive_energy_jump_limit_j` reaches the GPU (`fortranData.cpp` ->
`gpuSimulation.cpp::adaptiveEnergyJumpLimitJ`) but, before this task, was
read only to be printed in the final diagnostic line
(`gpuSimulation.cpp:957`, historical); it was never compared against
anything. The final summary line also hardcoded
`rejected_transitions=0` (`gpuSimulation.cpp:983`, historical) — a
correct number (zero candidates were ever energy-rejected, because none
were ever energy-*evaluated*) presented with no indication that no
criterion had run at all.

`publishProposedState` does accept an `acceptedBlockMask` parameter, so the
publication step already has a hook for per-block accept/reject; nothing
upstream of it computes one.

**Conclusion confirmed by re-reading the code**: GPU resolution changes were,
before and after this task, unconditionally selector-driven. CGP-00B's Part F
"Gate OFF" contract (zero transition-energy work) was therefore already true
on GPU; the work this task adds there is making that fact an explicit,
honest, and load-bearing part of the contract (rather than an accident of
"nobody wired the mask yet"), plus a hard startup refusal for the one
configuration (`cg_energy_jump_gate Y` + `do_gpu`/`do_gpu_llg`) that would
otherwise silently ask for a safeguard that cannot be honoured.

## 2. Design implemented

* New keyword `cg_energy_jump_gate` (`Y`/`N`, default `N`), parsed in
  `source/Input/inputhandler.f90`, stored in
  `InputDataType::adaptive_cg_config_t` alongside a `..._set` flag recording
  whether the keyword was actually present in the input file (the same
  technique is applied to the pre-existing `cg_energy_jump_limit`).
* Resolution logic (`AdaptiveCGProduction::validate_configuration`):
  * `cg_energy_jump_gate` present -> use it verbatim.
  * absent, `cg_energy_jump_limit` present, **not** a GPU run -> gate=ON
    (preserves the historical CPU meaning) with a migration notice.
  * absent, `cg_energy_jump_limit` present, **GPU run** -> gate stays OFF
    (preserves the historical GPU meaning, where the value was already
    inert) with a diagnostic stating the value is ignored. This asymmetry
    is deliberate: applying the CPU inference rule uninformed of backend
    would have turned already-passing GPU fixtures
    (`tests/coarse_graining/e2e/gpu_adaptive_mixed`, which sets
    `cg_energy_jump_limit 1.0e-18` and has no `cg_energy_jump_gate`) into
    new hard failures, which is exactly the kind of surprise Part C asks to
    avoid.
  * neither present -> gate=OFF (the new default).
  * gate resolves ON -> `cg_energy_jump_limit` must have been supplied and
    be finite/nonnegative, and the run must not be a GPU run (see below).
* `AdaptiveHybridSolver::adaptive_reconstruction_configuration_type` gained
  `energy_jump_gate` (default `.false.`). `apply_adaptive_transitions` now
  evaluates `energy_before`/`energy_after` only inside
  `if (reconstruction_configuration%energy_jump_gate)`; with the gate off, a
  transition that survives selector/reconstruction/ownership-rebuild is
  published unconditionally, exactly as GPU already did.
* A new counter, `adaptive_hybrid_runtime_type%transition_energy_evaluations`
  (`int64`, monotonically increasing for the runtime's lifetime), increments
  once per `energy_evaluator` call. It is printed in
  `AdaptiveCG: energy_jump_gate=... transition_energy_evaluations=...`.
* Each logged transition event now carries `energy_gate_applied`
  (`adaptive_transition_event_type%energy_gate_applied`), so a `0.0`
  `energy_jump_j` in the log can be told apart from "gate was off, this
  field was never measured."
* GPU (`source/gpu_files/gpuSimulation.cpp`) diagnostics were reworded (no
  kernel/runtime code changed) to say `energy_jump_gate=disabled` plainly
  instead of printing a configured-but-inert limit next to "resolved", and
  to print `transition_energy_evaluations=0` alongside the pre-existing
  `rejected_transitions=0`.

## 3. Why GPU Part G (full-FP64 candidate energy) was not implemented

`GpuAdaptiveRuntime::evaluateHybrid` (used for ordinary measurement energy)
computes every local energy contribution in the kernels' native working
type (`real`, i.e. FP32 in a SINGLE build) and only widens the *result* to
`double` immediately before the deterministic FP64 reduction — e.g.
`evaluateAdaptiveAtomisticBonds`:

```cpp
real si[3], sj[3], ksiJ[3] = {}, ktSi[3] = {};
...
energy = -static_cast<double>(dotDevice(si, ksiJ));   // FP32 dot product, widened after
```

This is precisely the "local energies evaluated in FP32, then widened"
pattern CGP-00 found cannot recover a small energy difference (see
`docs/CGP-00_ENERGY_PRECISION_EVIDENCE.md`), and CGP-00B Part G explicitly
forbids using it for the transition-acceptance decision.

A correct implementation needs the *local* arithmetic itself — the bond
dot/cross products, the per-neighbour accumulation, the on-site
anisotropy/spiralization/dipole formulas — carried out in `double`, not
merely a `double`-typed accumulator fed FP32 products. That touches every
energy-producing kernel in `gpuAdaptiveRuntime.cpp`
(`evaluateAdaptiveAtomisticBonds`, `evaluateAdaptiveAtomisticOnsite`, the
coarse exchange/spiralization kernel, the coarse anisotropy/external kernel,
and both dipole kernels — roughly 700 lines spanning every accepted
Hamiltonian term CG-10.5 supports), each of which is presently fused with
its field-writeback side effect and hard-codes `real` for its local
temporaries.

Per this task's own escalation clause and the governing CGP contract's rule
against a second, independently maintained energy implementation, this
subpart is **not implemented** and is flagged for Terra review. Instead:

* `cg_energy_jump_gate Y` together with `do_gpu Y`/`do_gpu_llg Y` is
  rejected at startup (`validate_configuration`), with a diagnostic naming
  this exact gap, rather than silently running the FP32-widened
  `evaluateHybrid` path as if it were the safeguard.
* The GPU diagnostic output states plainly that the safeguard is not
  implemented rather than implying it might be active.

A future implementation should introduce a working-precision template
parameter on these kernels (`float`/`double`), mirroring the
`evaluateHybridImpl<MeasureEnergy=...>` compile-time-specialization pattern
already accepted for CGP-01, and instantiate it a second time with
`double` for exactly this candidate-only, energy-only, no-field-writeback
evaluation. That is a Terra-scale change to the accepted Hamiltonian kernels,
not a Luna-scale patch.

## 4. Correctness evidence

All commands below were run from a rebuilt `build_cpu` (`CMAKE_BUILD_TYPE
Release`, `UPPASD_PRECISION DOUBLE`, `UPPASD_GPU_BACKEND OFF`) and
`build_gpu_fp32` (`UPPASD_GPU_BACKEND CUDA`, `UPPASD_PRECISION SINGLE`) tree
already present in the working copy; both were incrementally rebuilt after
every source change described here and both compiled without warnings
attributable to this change.

### 4.1 Unit-level (`tests/coarse_graining/test_adaptive_hybrid_solver.f90`,
target `adaptive_hybrid_solver_tests`)

Extended `test_transactional_chain_transitions` with CGP-00B-specific
assertions (kept alongside its pre-existing transactional-rollback
coverage, which now explicitly sets `energy_jump_gate = .true.` where it
forces a synthetic energy-jump rejection):

* **Test 5 (selector authoritative)**: gate enabled, limit `0.0` (tightest
  possible), no selector proposal (`requests%refine/coarsen` all `.false.`)
  -> `status == ADAPTIVE_HYBRID_OK`, resolution state unchanged,
  `transition_energy_evaluations == 0`, no log event created.
* **Gate-ON veto**: gate enabled, limit `0.5d-18`, a forced synthetic
  energy jump -> transition rejected, state/epoch/directions/ownership
  rolled back exactly as before, and now additionally:
  `transition_log%event(1)%energy_gate_applied .and.
  transition_energy_evaluations == 2` (one `energy_before` + one
  `energy_after` call for the one proposed candidate).
* **Test 1 (gate-off zero cost)**: gate reset to `.false.` (the default),
  six further accepted refine/coarsen transitions over three cycles ->
  `transition_energy_evaluations` stays at exactly `2` (unchanged from the
  gate-ON call above) and none of those six events has
  `energy_gate_applied`.

`ctest -L coarse-graining` (30 tests, includes this target plus every other
CG-10.5/RCG unit and e2e fixture already tracked): **30/30 passed**, no
regressions. Full `ctest` (38 tests, includes non-CG dipole/regression
fixtures): 37/38 passed; the one failure
(`dipole-open-fft-oracle::test_energy_finite_difference_derivative`) is a
pre-existing, floating-point-tolerance finite-difference assertion in the
unrelated open-boundary dipole FFT kernel oracle
(`tests/dipole_validation/test_open_fft_oracle.py`), reproduced bit-for-bit
identically (`-1.155179279520103e-25 != -1.1551792784777866e-25`) against a
build from the unmodified `HEAD` tree (verified via `git stash`
before/after) -- it is untouched by, and unrelated to, this task's scope.

### 4.2 Backward compatibility (`tests/coarse_graining/run_production_e2e.py`,
`run_moving_backend_parity.py`, `run_energy_jump_threshold_sweep.py`)

* `run_production_e2e.py --binary <cpu> --gpu-binary <gpu fp32>`: every
  case up to and including `CG-10.5 DMI/anisotropy CPU/GPU production
  parity passed` succeeds unchanged; the only failure
  (`gpu_fft_static_mixed`'s `coarse_dipole` assertion) reproduces
  identically on the unmodified `HEAD` tree (same `git stash`
  before/after check) and is unrelated to adaptive-CG transitions.
* `run_moving_backend_parity.py --cpu-binary <cpu> --cuda-fp32-binary <gpu
  fp32>`: 19/19 fixtures on each backend, all fp32 parity budgets met,
  all four negative controls still discriminating. `moving_wall_adaptive`
  (the one tracked fixture with real physics-driven transitions; its
  `inpsd.dat` sets `cg_energy_jump_limit 1.0e-15` with no
  `cg_energy_jump_gate`, so CPU infers gate=ON while GPU stays
  selector-only) reports `accepted_transitions cpu=12 gpu=12` — unchanged
  from the pre-CGP-00B contract, because every real jump in that fixture
  (`run_energy_jump_threshold_sweep.py`: max `|energy_jump_j|` =
  `1.109571e-18`) sits far inside the configured `1.0e-15` limit.
* `run_energy_jump_threshold_sweep.py`: unchanged output
  (`accepted=12 rejected=0 (energy-jump-rejected=0)`), confirming the
  legacy-inference path reproduces the exact same transition history as
  before this task for every currently tracked fixture.

### 4.3 GPU gate-ON refusal

Manually amending `tests/coarse_graining/e2e/gpu_adaptive_mixed/inpsd.dat`
(which already carries `cg_energy_jump_limit 1.0e-18`) with an explicit
`cg_energy_jump_gate Y` and running the `build_gpu_fp32` binary produces:

```
AdaptiveCG setup rejected: cg_energy_jump_gate: the GPU adaptive-CG path
has no full-FP64 transition-energy evaluation yet (CGP-00B Part G is not
implemented; flagged for Terra review). Run this configuration on the CPU
path, or disable the gate.
```

Running the same file *without* that explicit line (i.e. exactly as
tracked) prints the legacy-inference diagnostic and proceeds normally:

```
AdaptiveCG: cg_energy_jump_limit was supplied but the GPU adaptive-CG path
has no transition-energy gate (CGP-00B Part G is not implemented); the
value is ignored, exactly as before this input became explicit.
Transitions remain selector-driven only.
```

## 5. Performance evidence (Test 1/3/4/H, "Required acceptance")

Fixture: `moving_wall_adaptive` (192 atoms, 900 steps, `build_cpu`
DOUBLE), the only tracked fixture with physics-driven transitions,
run three ways from the same starting state:

| Config | `cg_energy_jump_gate` | `cg_energy_jump_limit` (J) | `transition_energy_evaluations` | `accepted_transitions` | `rejected_transitions` (energy) | wall (s) |
|---|---|---:|---:|---:|---:|---:|
| gate off | `N` | (ignored) `1.0e-15` | **0** | 12 | n/a (0, disabled) | 0.2601 |
| gate on, generous | `Y` | `1.0e-15` | 14 | 12 | 0 | 0.2824 |
| gate on, tight | `Y` | `1.0e-19` | 2710 | 10 | 1800 |1.9164 |

(`wall` is the `AdaptiveCG: phase_wall_seconds wall=` figure, an independent
host wall-clock measurement of the whole adaptive step loop, not summed
phase timers; `transition_energy_evaluations` and the accepted/rejected
counts are the exact printed `AdaptiveCG:` summary values.)

Reading this evidence:

* **Gate off**: `transition_energy_evaluations=0` for the entire 900-step,
  12-real-transition run — the zero-cost path holds exactly as required,
  and this is the CPU DOUBLE case; the GPU case is architecturally
  incapable of transition-energy work at all in the current code (§1), so
  it is zero unconditionally, not merely by configuration.
* **Gate on, generous** (the limit this fixture has always shipped with):
  same 12 accepted transitions as gate-off — the gate changes nothing for
  a threshold well above the real jumps, at a cost of 14 extra
  `energy_evaluator` calls (root cause of the ~8.6% wall-time delta at this
  small scale: two calls are effectively free relative to a 192-atom
  hybrid-energy evaluation, but this system is small enough that even that
  shows up).
* **Gate on, tight** (deliberately set below the largest observed real
  jump, `1.109571e-18`, from
  `run_energy_jump_threshold_sweep.py`): 1800 candidates are proposed and
  rejected (`outcome=energy-jump-rejected` appears 1800 times in stdout)
  as the selector keeps re-proposing the same wall-adjacent block every
  dwell interval; only 10 of the 12 "successful" transitions from the
  generous run still clear the tighter bar. The run completes all 900
  steps and both preceding CPU/GPU parity/threshold-sweep scripts (§4.2)
  ran against binaries built from this same tree, i.e. later LLG dynamics
  are unaffected by the roundtrip through many rejected candidates. The
  ~7.4x wall-time increase here is real and expected: it is proportional to
  the actual number of *proposed* candidates (2710 evaluator calls across
  1800+10 proposals), not to the 900 timesteps, exactly the "amortize
  according to actual transition frequency" contract Required work item H
  and the Performance evidence section ask for.

No GPU-side performance evidence is reported for gate=ON, because that
configuration is refused at startup (§3/§4.3); GPU gate=OFF has no
transition-energy work to measure by construction (§1).

## 6. Checklist

* [x] CPU transition state machine audited (§1).
* [x] GPU transition state machine audited (§1).
* [x] CPU/GPU semantic difference documented before modification (§1).
* [x] Explicit transition-energy-gate setting introduced
  (`cg_energy_jump_gate`).
* [x] Gate defaults OFF.
* [x] `huge()` is no longer the semantic mechanism for disabling the gate.
* [x] Legacy input behaviour handled explicitly and documented (§2, backend
  asymmetry deliberate and justified).
* [x] Normal selector remains the sole source of transition proposals.
* [x] Energy cannot independently trigger refinement (§4.1 Test 5).
* [x] Energy cannot independently trigger coarsening (§4.1 Test 5, same
  code path).
* [x] Gate-OFF CPU path performs zero transition-energy evaluations (§4.1,
  §5).
* [x] Gate-OFF GPU path performs zero transition-energy evaluations (§1,
  architecturally -- no code path evaluates transition energy at all).
* [x] Ordinary measurement energy remains available independently (not
  touched by this task; `production_energy_evaluator`/`evaluateHybrid` are
  unmodified).
* [x] Measurement cadence does not alter gate enablement (gate is a static
  per-run configuration value, never read from measurement state).
* [x] Gate-ON energy is evaluated only for actual candidate transitions
  (§4.1, §5: 0 evaluations with 0 proposals, 2/14/2710 evaluations scaling
  with 1/12/1810 proposals).
* [x] CPU enabled-gate rollback remains correct (§4.1, pre-existing
  transactional test retained and now gate-explicit).
* [x] GPU enabled-gate accept/reject semantics match CPU: **not
  applicable** -- GPU gate=ON is refused at startup rather than
  implemented (§3); this is the one flagged, undone subpart.
* [ ] SINGLE GPU enabled gate uses full-FP64 local + reduction evaluation:
  **not implemented, flagged for Terra review** (§3).
* [x] No second independent Hamiltonian implementation introduced (none
  was written; §3 explains why one would have been required and why that
  was stopped instead).
* [x] Selector-only transition fixture passes (§4.1 Test 5; also every
  `active_block_updates=0`-selector-quiet step in §5's gate-off/on-generous
  runs).
* [x] Permissive-gate fixture passes (§5 "gate on, generous"; §4.2
  `moving_wall_adaptive` cpu=gpu=12).
* [x] Tight-gate rejection fixture passes (§5 "gate on, tight", 1800 real
  rejections, run completes).
* [x] No-proposal fixture performs zero transition-energy work (§4.1 Test
  5).
* [x] CPU/GPU transition-history fixture passes (§4.2; exact accept/reject
  count match for the only regime GPU can run, gate off/legacy-inferred
  vs. GPU's always-selector-only history).
* [ ] Near-threshold precision fixture passes: **not applicable to GPU**
  because GPU gate=ON does not exist (§3); the CPU path was already full
  `dblprec` before and after this task (confirmed in
  `run_energy_jump_threshold_sweep.py`'s own module docstring and by
  inspection of `apply_adaptive_transitions`, §1), so CGP-00's FP32-local
  precision-loss concern never applied to CPU transition acceptance in the
  first place.
* [ ] Negative precision control is discriminating: **not applicable**, no
  GPU gate=ON implementation exists to control against (§3).
* [x] Gate-OFF production timestep has no measurable energy-gate overhead
  (§5: `transition_energy_evaluations=0` for the entire run; the CPU
  wall-time noise between gate-off and gate-on-generous at this small
  192-atom scale is 14 evaluator calls total across 900 steps, not a
  per-step cost).
* [x] Documentation distinguishes selector, energy safeguard, and
  measurement (`docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` §7.7.1,
  `examples/AdaptiveCoarseGraining/README.md`).
* [x] Existing all-fine Depondt parity tests remain unchanged (`ctest -L
  coarse-graining`, 30/30; `run_moving_backend_parity.py`'s
  `moving_all_fine`/`moving_all_fine_wide` fixtures unaffected).
* [x] Existing finite-temperature adaptive tests remain unchanged (CG-10.5
  finite-temperature fixtures in `run_production_e2e.py` pass unchanged).

## 7. Remaining gaps

* GPU transition-energy gating (CGP-00B Part G) is not implemented. See §3
  for the audit of why, and the recommended template-specialization
  approach for whoever picks this up next (Terra-scale).
* The backward-compatibility inference is intentionally backend-asymmetric
  (§2); this is documented in
  `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` §7.7.1 and above so it is
  not mistaken for an oversight.
