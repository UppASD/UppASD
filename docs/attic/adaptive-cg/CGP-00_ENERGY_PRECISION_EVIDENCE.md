# CGP-00 — Adaptive-CG energy-precision evidence

Task: `docs/CGP_work.md` CGP-00 ("Re-evaluate the energy-precision contract").
Evidence-only task; no production precision behaviour was changed. The one
production source change (extracting the block-tree reduction into a
templated, shared helper) is mechanical and behaviour-preserving — see
"Production code change" below.

Branch: `gpu_hip_cu_ab_cg`, commit base `0c8fb2b8`.

## Environment

- GPU: 2x NVIDIA RTX A4000, driver 610.57.04. Both idle (0% utilization,
  <110 MiB used) for the full duration of this evidence run — see the
  shared-host-contention note in project memory; no interleaving/contention
  correction was needed here because CGP-00 makes no timing claims.
- CUDA 13.3.73 (nvcc), gfortran 13.3.0, GCC 12.4.0, CMake 3.28.3.
- Builds: fresh out-of-tree `build_cgp00_cuda_fp64` (DOUBLE_PREC, CUDA),
  `build_cgp00_cuda_fp32` (SINGLE_PREC, CUDA), `build_cgp00_cpu` (CPU-only),
  each configured with `BUILD_TESTING=ON`, `CMAKE_BUILD_TYPE=Release`.

## Production code change (Part A infrastructure, not a precision change)

`reduceAdaptiveEnergyBlock` (the RCG-09B per-block shared-memory tree
reduction, `gpuAdaptiveRuntime.cpp`) was extracted into a templated
`reduceAdaptiveEnergyBlockT<Acc>` in `gpuAtomicDouble.hpp` — the header
already shared between production and the RCG-06B `test_energy_fp32_accum.cpp`
fixture for the identical reason (a fixture must exercise the literal
production primitive, not a reimplementation that can drift). Production's
`reduceAdaptiveEnergyBlock(double, ...)` now forwards to
`reduceAdaptiveEnergyBlockT<double>` — byte-identical codegen, all 7 existing
call sites unchanged. This is what lets the new
`test_energy_hierarchical_precision.cpp` fixture (below) instantiate the
identical block-tree primitive at `Acc=float` without reimplementing it.

Verified non-regression: full `cg13-cuda` suite (33 tests) passes on the
DOUBLE_PREC build after this change; see "Correctness" below.

## Part A/B — Hierarchical reduction precision, synthetic sums

New fixture: `tests/coarse_graining/test_energy_hierarchical_precision.cpp`
(ctest `adaptive-cg-energy-hierarchical-precision`, labels
`cg13`/`cg13-cuda`/`cg13-hip`). Deliberately separate from
`test_energy_fp32_accum.cpp`, which stays untouched as CGP-00 instructs —
that fixture is protected historical evidence for the *older*, already-fixed
flat-atomicAdd defect class (a single globally-contended accumulator), and
comparing it against the current hierarchical topology would not be an
apples-to-apples FP32-vs-FP64 comparison, which CGP-00 explicitly warns
against.

The new fixture reproduces the current three-stage topology — per-thread
local contribution → block-tree reduction (the shared
`reduceAdaptiveEnergyBlockT`) → ordered serial sum of block partials — and
sweeps the accumulator type independently at each stage:

| Mode | local | block | final | matches |
| --- | --- | --- | --- | --- |
| 1 | F32 | F32 | F32 | (never shipped) |
| 2 | F32 | F32 | F64 | (never shipped) |
| 3 | F32 | F64 | F64 | today's SINGLE_PREC GPU build |
| 4 | F64 | F64 | F64 | today's DOUBLE_PREC GPU build / CPU reference |

Reference: Neumaier-compensated summation in `long double` on the host. The
reduction is fully deterministic (no atomicAdd anywhere in this topology), so
every (case, N, mode) result is exact and reproducible — no repeats needed.

Six term distributions, N swept 100 → 3,000,000 (matching
`test_energy_fp32_accum.cpp`'s existing top end):

**All-positive** (continuity check with the older fixture): mode 1 error
grows with N as expected (0 → 9.9e-5 relative at N=3e6); modes 2/3 stay flat
around 2.5e-9 absolute regardless of N; mode 4 stays ≤2.6e-13 relative.

**Alternating sign**: trivial (exact 0 error every mode) — a sanity baseline,
not informative on its own.

**Severe cancellation** (large near-equal ±1e6 pairs riding a tiny known
residual): **modes 1, 2, and 3 all show ~100% relative error at every N** —
not because the block/final reduction is bad, but because casting the
*local* term (1e6 + residual) down to `float` before it ever reaches the
reduction already rounds the residual away (float has ~7 decimal digits;
1e6 + 1e-5 rounds to 1e6 exactly). This is CGP-00's own background
distinction #1 ("precision of each local energy contribution") dominating
over distinction #2 ("precision of the reduction") for this case — a real,
measured finding, not an artifact. Mode 4 (F64 local) retains the residual,
with its own relative error set by ordinary cancellation loss-of-significance
(~1e-8 at N=100 down to ~3e-10 at N=3e6) — expected floating-point behaviour
for subtracting O(1e6) terms to get an O(1e-2..1e2) residual, not a defect.

**Wide dynamic range** (±10^[-6,6] terms): mode 1 relative error 1.8e-7 →
3.7e-7 (worsening with N); modes 2/3 improve on mode 1 by roughly 1–2 orders
of magnitude but are not always monotonic in N; mode 4 stays ≤3.6e-15
relative throughout.

**Randomized physical** (order-1 Gaussian terms, seeded): mode 1 relative
error grows to 5.2e-6 at N=3e6; mode 3 (today's SINGLE_PREC contract) stays
at 1.9e-8, roughly **280x better than mode 1** — this is the clearest
evidence that FP64 reduction over F32 local contributions is a real,
substantial, non-cosmetic win for ordinary (non-cancellation) measurement
sums; mode 4 stays ≤1.1e-15 relative.

**Small-difference-of-two-large-sums** (`energy_difference`, N=1,000,000,
delta = 0.1, 1e-3, 1e-6, 1e-9 against O(10)-magnitude Gaussian terms — the
quantity a transition accept/reject decision actually depends on, not either
total in isolation):

| delta | mode 1 rel err | mode 2 | mode 3 (today, SINGLE) | mode 4 (today, DOUBLE) |
| --- | --- | --- | --- | --- |
| 0.1 | 3.9e-3 | 6.1e-5 | 2.4e-7 | 1.5e-11 |
| 1e-3 | 1.00 (destroyed) | 7.1e-3 | 4.7e-5 | 2.0e-10 |
| 1e-6 | 1.00 (destroyed) | 1.00 (destroyed) | 4.6e-2 | 3.4e-7 |
| 1e-9 | 1.00 (destroyed) | 1.00 (destroyed) | 1.00 (destroyed) | 4.4e-4 |

The pattern is consistent with the cancellation case: **today's SINGLE_PREC
GPU contract (mode 3) does not protect small energy differences** — an
F32-local, F64-reduced evaluation loses essentially all information once the
delta drops much below ~1e-3 relative to the constituent term magnitudes,
because the loss already happened before the FP64 reduction ever ran. Only
mode 4 (F64 local *and* reduction) tracks the delta down to the
smallest tested value.

**Required negative controls (both satisfied):**
- FP64 reference mode stays within 1e-9 relative error on every
  well-conditioned case (all_positive, wide_dynamic_range,
  randomized_physical) at every N — worst observed 2.6e-13. (Cancellation and
  energy_difference are excluded from this specific gate: their relative
  error relative to the *residual* is inherently large under exact
  arithmetic too, for reasons explained above — the fixture's own comment
  block documents why blending that into one blanket relative-error budget
  would be a wrong test, not a real finding.)
- severe_cancellation and every energy_difference case discriminate
  F32/F32/F32 from F64/F64/F64 by absolute error at every tested point — the
  fixture asserts this and fails loudly if it ever stops holding.
- The older `test_energy_fp32_accum.cpp` (flat atomicAdd defect class) is
  untouched and still passes, still showing the previously-established poor
  FP32 scaling.

## Part C — Production energy terms

Two lower-risk substitutes for building four live precision-variant
instantiations of every production kernel (see the plan's reasoning: that
would be exactly the kind of production churn CGP-00's "low production-code
risk" framing warns against, and duplicates work CGP-01 needs to do properly
for an unrelated reason):

**1. New accessor, validated against real production kernel output.**
`GpuAdaptiveRuntime::downloadEnergyPartialsSnapshot()` (new, evidence-only,
not on the production timestep's data path) downloads the real per-block FP64
partials each production energy kernel wrote for the most recent
`evaluateHybrid()` call — atomistic exchange+DMI-folded bond energy, onsite
anisotropy, coarse exchange, coarse spiralization, coarse anisotropy, coarse
external, and dipole. A new check in `test_gpu_adaptive_runtime.cpp`
(`testKernelParityAndWorkflow`) downloads these partials right after the
existing, already-CPU-reference-checked `evaluateHybrid()` call and asserts
the ordered sum of each term's real partials reproduces `evaluateHybrid`'s
own cached total to the same tolerance used elsewhere in that test
(2e-12 fp64 / 2e-5 fp32) — passing on both DOUBLE_PREC and SINGLE_PREC
builds. This makes real per-term production magnitudes available to future,
larger-scale precision work without new production-code risk.

**2. Whole-run production-scale energy comparison, rerun fresh.**
`run_energy_fp32_accum_production.py` (RCG-06B, unmodified) was rerun against
freshly built CPU/CUDA-fp64/CUDA-fp32 binaries from this branch, across all
19 tracked `moving_*` fixtures (48 and 192 atoms):

```
CUDA fp64: natom=48  mean=1.210e-15  max=1.210e-15
           natom=192 mean=5.728e-06  max=8.592e-05
CUDA fp32: natom=48  mean=1.428e-07  max=1.428e-07
           natom=192 mean=6.068e-06  max=8.590e-05
```

fp64 worst (8.592e-05) and fp32 worst (8.590e-05) are essentially identical —
reproducing the "fp32 error tracks fp64 error" signature this script's own
docstring already documented from RCG-06B, now freshly re-confirmed
*including this task's `reduceAdaptiveEnergyBlockT` refactor*, i.e. no
regression at production scale. As that docstring already explains, this
whole-run comparison is dominated by CPU/GPU summation-order (associativity)
differences at these two sizes, not by accumulator precision — which is
exactly why Part A/B's controlled, same-topology comparison (above) is the
load-bearing evidence for this task, not this whole-run number.

**3. Real production term counts vs. the Part B sweep range.** The existing
RCG-09 evidence (`docs/RCG-09_PRODUCTION_PERFORMANCE_EVIDENCE.md` §3.5,
65,536 atoms) records up to 114,688 live bonds in a single energy-kernel
launch at 100% fine fraction. Part B's swept N range (100 to 3,000,000)
comfortably brackets this real production launch scale.

## Part D — Transition-decision stability

Two structural facts had to be established by reading the code before this
part could be scoped correctly:

1. The energy-jump accept/reject gate (`abs(energy_jump) <=
   energy_jump_limit_j`) exists in exactly one place in this codebase:
   `source/CoarseGraining/adaptivehybridsolver.f90`'s
   `apply_adaptive_transitions`, the **CPU-only** adaptive-CG path.
   `energy_before`/`energy_after`/`energy_jump` there are computed entirely
   in `dblprec` (Fortran double) — there is no FP32 anywhere in that call
   chain, so it is structurally unaffected by this task's local/reduction
   precision question.
2. The **GPU** adaptive-CG runtime (`gpuAdaptiveRuntime.cpp`, this phase's
   actual primary scope) has **no energy-jump gate at all**:
   `proposeSelectorState`/`publishProposedState` decide and publish
   transitions from selector/polarization thresholds only, and never
   evaluate or compare an energy jump. `adaptive_energy_jump_limit_j` is read
   GPU-side (`gpuSimulation.cpp`) only to print it in a diagnostic line.
   Independently corroborated by this repository's own RCG-04I docstring
   (`run_moving_backend_parity.py`): "the GPU AdaptiveCG diagnostic path ...
   prints only initial/final resolution-state snapshots, never a per-event
   transition log, and its one aggregate summary line hardcodes
   `rejected_transitions=0`" — exactly what "no gate" looks like from
   outside, and reproduced live in this evidence run (see below).

**So there is currently no live GPU decision for a reduction-precision
change to flip.** New script `tests/coarse_graining/run_energy_jump_threshold_sweep.py`
instead extracts real `energy_jump_j` values from the one tracked fixture
whose adaptive transitions are driven by real, physically-moving-wall
dynamics rather than a synthetic construction (`moving_wall_adaptive`, 192
atoms, 900 steps, `cg_energy_jump_limit_j=1e-15`), run through the CPU
binary (the only path that actually evaluates the gate):

```
transition events parsed: 12, accepted=12 rejected=0 (energy-jump-rejected=0)
|energy_jump_j|: min=0 median=0 max=1.110e-18
jump / limit ratio: max=1.110e-3; 0/12 events within 10x of the boundary
```

Every real observed jump in this fixture sits far below the configured
limit (all coarsen-requests; max ratio to the limit ~0.1%) — this run does
not itself stress the boundary, which is an honest evidence outcome, not a
gap to paper over. Combined with the finding above, Part A/B's synthetic
`energy_difference` sweep (which *does* sweep down through and past a tight
boundary) is the load-bearing evidence for this part: it shows that if/when
a GPU-side energy-jump gate is added, mode 3 (today's SINGLE_PREC contract)
would not reliably protect a decision once the jump drops much below
~1e-3 relative to the terms being differenced, while mode 4 (full double,
matching what the existing CPU-only gate already does) would.

## Correctness

- `ctest -L cg13-cuda` (DOUBLE_PREC build): **33/33 passed**, including the
  new `adaptive-cg-energy-hierarchical-precision` fixture (0.84s) and the
  untouched `adaptive-cg-energy-fp32-accum` (112.19s, still passing).
- `ctest -L cg13-cuda` (SINGLE_PREC build): one pre-existing, unrelated
  failure (`adaptive-cg-production-e2e`, `assert abs(float_metric(fft.stdout,
  "coarse_dipole")) > 0.0`) — **confirmed present on the unmodified base
  commit** (verified by stashing every CGP-00 change, rebuilding, and
  reproducing the identical failure before restoring the stash), so it is
  not a CGP-00 regression. Not investigated further; out of scope for an
  energy-precision evidence task. Every other test in the SINGLE_PREC
  `cg13-cuda` label set passed.
- CPU build (`build_cgp00_cpu`): used for Part C/D evidence generation;
  19/19 tracked fixtures ran successfully.

## Recommendation

Evidence supports **Contract B** — native-real hierarchical reduction is
sufficient for ordinary measurement energies, while transition-acceptance
energy retains a full FP64 evaluation — **not** Contract A (native
everywhere) and **not** Contract C (FP64 everywhere) as literally stated,
with one evidence-driven refinement to how Contract B is implemented:

- **Measurement energies**: FP64 *reduction* over native-real (F32 in
  SINGLE_PREC builds) local contributions is a real, substantial win for
  ordinary (non-cancellation-dominated) sums — up to ~280x lower relative
  error than a fully-native reduction in the randomized-physical case, and
  1–2 orders of magnitude in the wide-dynamic-range case. Moving measurement
  energy to Contract A (fully native) would measurably regress today's
  accuracy; nothing in this task's evidence justifies that regression for
  speed (this task made no timing measurements at all, deliberately).
  Today's FP64-reduction behaviour for measurement should be kept.
- **Transition-acceptance energy**: the evidence shows Contract B's literal
  wording ("retains a final FP64 reduction") is *necessary but not
  sufficient* on its own. An FP64 reduction over F32 *local* contributions
  (mode 3 — today's actual SINGLE_PREC GPU contract, the only kind of
  "FP64 reduction" a native-real field-evaluation kernel can produce) still
  loses essentially all signal once a jump drops below roughly 1e-3 relative
  to the terms being differenced (see the `energy_difference` table above).
  Protecting a transition decision near a tight threshold requires FP64
  *local* contributions too for the specific terms entering that decision —
  exactly what the existing CPU-only `apply_adaptive_transitions` path
  already does (dblprec throughout), and exactly why it should stay that
  way. This is a concrete design constraint for any future work that adds an
  energy-jump gate to the GPU transition path (there is none today, per Part
  D): such a gate must evaluate its candidate energies in full double
  precision, not merely reduce native-real per-thread contributions into a
  double accumulator.

No production precision behaviour was changed in this task. Adopting this
contract for the GPU field-only/energy-split work is CGP-01's task, subject
to Human acceptance of this evidence.

## Checklist (docs/CGP_work.md CGP-00)

- [x] Same-topology FP32/FP64 comparison implemented.
- [x] Positive-sum scaling tested.
- [x] Alternating-sign sum tested.
- [x] Severe cancellation tested.
- [x] Wide-dynamic-range sum tested.
- [x] Multi-million-term scale tested.
- [x] Production energy terms compared (accessor validated against real
      kernel output; whole-run comparison rerun fresh; real term counts
      cross-checked against the swept N range).
- [x] Energy-per-atom error reported (`per_term_err` column in the fixture's
      own output).
- [x] Small energy-difference accuracy reported.
- [x] Transition decisions tested near the threshold (synthetic sweep; real
      observed jumps in the one available real-transition fixture were far
      from the threshold, honestly reported as such; GPU path has no live
      gate today, also reported).
- [ ] fp32 and fp64 performance cost measured separately — **not attempted**;
      this task is evidence-only for precision/error, and CGP_work.md's own
      common protocol scopes performance measurement to tasks that change
      production runtime, which this one does not.
- [x] Recommended precision contract follows measured evidence.
- [x] No production precision behaviour changed in this task.

## Commit

`CGP-00: re-evaluate adaptive energy precision`
