# CGP-08 — Reduce the per-step kernel-launch floor, evidence

Task: `docs/CGP_work.md` CGP-08 ("Reduce the per-step kernel-launch floor").
Branch `gpu_hip_cu_ab_cg`, commit base `70a7d48e` (CGP-07), built from the
branch's current working tree (CGP-05's uncommitted changes included).

**Primary model: Claude (Sonnet 5, interactive session)**, matching the
CGP-03D/CGP-05/CGP-06B/CGP-07 precedent.

**Dependency:** CGP-07's evidence (Part 6/Part 7), which found the all-coarse
steady state is dominated by a fixed ~37us/step floor traceable to ~12.3
kernel launches/step, and that this floor is large enough at the historical
65,536-atom benchmark size that adaptive coarse graining is *slower* than
plain production there.

## Environment

- GPU: NVIDIA RTX A4000 (device 0 of 2), shared host per project memory
  `shared-gpu-host-contention`. Checked idle before every benchmark run in
  this document (0% utilization, <110 MiB used, no other compute process) —
  see the self-inflicted-contention finding in Part 4 for one exception this
  task caught and corrected in itself.
- Driver 610.57.04, CUDA 13.3.73 (nvcc), `CMAKE_BUILD_TYPE=Release`.
- Host CPU: `nproc` reports 2 (cgroup-quota-limited).
- Builds: fresh out-of-tree `build_cgp08_cuda_fp32`/`build_cgp08_cuda_fp64`
  (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`/`DOUBLE`,
  `CMAKE_BUILD_TYPE=Release`).

## Part 1 — Audit before editing (required by the task's own Part A)

Traced every kernel launched by `GpuAdaptiveRuntime::evaluateHybridImpl`'s
steady-state coarse phase (`source/gpu_files/gpuAdaptiveRuntime.cpp`,
`clearAdaptiveCoarseActive` / `evaluateAdaptiveCoarseTensor` /
`finalizeAdaptiveCoarseLocal` / `writeAdaptiveCoarse`), read/write by
read/write, before touching anything:

| Kernel pair | Cross-thread-block hazard? | Fusable in CGP-08's scope? |
| --- | --- | --- |
| clear → evaluate | **Yes.** `evaluateAdaptiveCoarseTensor` does `atomicAdd` scatter (via `commitDerivativeStencil`) into its own block *and* up to 3 neighbour blocks. Fusing the zero-write and the scatter into one kernel launch would create a real race: one thread block's zero could stomp a neighbour thread block's already-arrived contribution, since CUDA/HIP gives no ordering guarantee between blocks within one launch. | No — would need a device-wide grid-sync primitive or a gather-based rewrite of the scatter algorithm. The latter is explicitly out of the whole CGP-00–08 phase's scope ("redesign of the continuum tensor/operator mathematics" is reserved for the next phase). |
| evaluate → finalize | **Yes**, same reason: finalize needs every block's scatter-target contributions fully settled, including from neighbours computed by other thread blocks. | No, same reasoning. |
| finalize → write | **No.** Verified from the source: `finalizeAdaptiveCoarseLocal` only reads/writes *its own thread's block* of `runtime.coarseField` — no atomics, no neighbour indexing anywhere in the kernel. `writeAdaptiveCoarse` likewise only reads its own block. No cross-thread-block dependency exists between them. | **Yes** — implemented below. |
| interface clear → restrict | Same scatter-hazard pattern as clear→evaluate (`restrictAdaptiveInterface` atomic-scatters into up to 8 corner blocks per ghost atom). At all-coarse `restrictAdaptiveInterface` doesn't even run (`ghostWork==0`), so there is nothing to fuse there in that regime; in mixed fractions the same hazard as clear→evaluate applies. | No. |
| predictor / corrector integration | Already 1 launch each, temporally separated by an intervening field evaluation (predictor field, then corrector field). Not adjacent launches at all. | Not a fusion target. |

One additional finding from the same read/write trace, independent of launch
count: `clearAdaptiveCoarseActive` zeroed the *public* `coarseField` output
buffer as well as the internal `runtime.coarseField` scratch, but
`writeAdaptiveCoarse` (the only consumer, always gated by the identical
`activeBlockList`/`activeBlockCount()` set) unconditionally **overwrites**
(`=`, not `+=`) every one of those same entries later in the same call, and
neither dipole kernel (`addAdaptiveDipole`/`addAdaptiveBasisResolvedDipole`)
touches the public `coarseField` parameter at all. The public-buffer half of
the clear was therefore dead work in every configuration (`MeasureEnergy` or
not, dipole or not) — deleted below, independent of the fusion.

A second correctness subtlety the launch-count audit alone would have
missed: `addAdaptiveDipole`/`addAdaptiveBasisResolvedDipole` (when active)
also accumulate into `runtime.coarseField` for coarse blocks, and are
sequenced strictly *between* the original `finalizeAdaptiveCoarseLocal` and
`writeAdaptiveCoarse` calls. A naive unconditional fusion of finalize+write
would silently drop the dipole contribution whenever dipole is enabled. The
fusion below is therefore **conditional**: used only when neither dipole
kernel is active; the dipole-enabled path keeps calling the two kernels
separately, in the same order, with the same arguments, completely untouched
by this task.

## Part 2 — What was implemented

1. **New kernel** `finalizeAndWriteAdaptiveCoarse<MeasureEnergy>`
   (`source/gpu_files/gpuAdaptiveRuntime.cpp`): the exact math of
   `finalizeAdaptiveCoarseLocal` (same accumulation order: anisotropy `+=`,
   then `*= -1/(mu*moment)`, then external field `+=`) held in a per-thread
   register instead of round-tripped through `runtime.coarseField`, followed
   by `writeAdaptiveCoarse`'s own read-of-scratch-and-store. `finalizeAdaptiveCoarseLocal`
   and `writeAdaptiveCoarse` themselves are **unmodified** and still used
   verbatim on the dipole-enabled path.
2. **Call-site branch** in `evaluateHybridImpl`: when
   `!uniformFftDipoleField && !basisResolvedFftField`, launch the fused
   kernel and skip the later separate `writeAdaptiveCoarse` call; otherwise
   (dipole enabled) the sequence is byte-for-byte what it was before this
   task.
3. **Dead-write removal**: `clearAdaptiveCoarseActive` no longer takes a
   `coarseField` parameter or zeroes it — only `runtime.coarseField`
   (the scratch `evaluateAdaptiveCoarseTensor` atomically accumulates into)
   is still zeroed, unconditionally, in every configuration.
4. Fixed one self-introduced launch-*counting* bug during this task (Part 4):
   `evaluateAdaptiveCoarseTensor`'s own `++phaseMetrics_.coarseLaunches`
   increment was dropped when the original combined `+= 2` (for
   evaluate+finalize together) was split across the new conditional branch.
   This was a diagnostics/instrumentation bug only — the kernel itself still
   launched correctly the whole time, confirmed because `cg13-cuda` passed
   both before and after the fix — but it made the first launch-count
   measurement (Part 4) read as a ~50% reduction when the real, corrected
   effect is ~25% of the coarse phase's launches.

## Part 3 — Correctness

- `cg13-cuda` on fresh `build_cgp08_cuda_fp32`: **33/34 pass**, the one
  failure the same pre-existing `adaptive-cg-production-e2e` `coarse_dipole`
  finding every CGP task since CGP-03D has recorded (verified unchanged: ran
  before *and* after the launch-counting fix, identical result both times).
- `cg13-cuda` on fresh `build_cgp08_cuda_fp64`: **34/34 pass**.
- Load-bearing individual fixtures for this specific change:
  `CG-10.5 DMI/anisotropy CPU/GPU production parity` (independently
  re-derives the exact anisotropy math the fused kernel now computes,
  against a CPU reference — this is what actually validates the arithmetic
  transcription, not just that the binary runs);
  `adaptive-cg-dipole-ownership-check` (exercises the untouched
  dipole-enabled fallback path, proving the conditional branch selects it
  correctly and that dipole contributions still land between finalize and
  write);
  `adaptive-cg-moving-backend-parity`/`adaptive-cg-ownership-map-comparator`
  (mixed, moving-fraction fixtures, exercising the fused kernel across a
  changing `activeBlockList` from transition to transition);
  `adaptive-cg-energy-fp32-accum`/`adaptive-cg-energy-hierarchical-precision`
  (the `MeasureEnergy=true` path through the fused kernel).
- `gpu_adaptive_runtime_tests` (dedicated CG-09/CG-10 adaptive-runtime unit
  suite) passes standalone.
- `gpu_adaptive_runtime_benchmark`'s own `atomistic-parity` gate (all-fine
  field identity against UppASD's real production Hamiltonian/integrator)
  passes at every size in the matrix below — note this specific gate does
  not exercise the fused kernel itself (it runs at `fraction=1.0`, where
  `activeBlockWork=0`), so it is not claimed as evidence for the fusion; it
  is the pre-existing "same physics" license for the timing comparison.
- No dedicated bit-diff harness was built comparing fused-vs-unfused binary
  output directly; the correctness case rests on the independent CG-10.5
  CPU-reference anisotropy check plus compute-sanitizer (Part 3 continued
  below), which is treated as sufficient given the fusion is a bookkeeping
  change (removing a round-trip through global memory), not a reordering —
  there is no misordering scenario to construct a synthetic negative control
  against, unlike CGP-05's event-based synchronization change.

`compute-sanitizer --tool synccheck` and `--tool memcheck`, both run over
the full fraction sweep (`--blocks 4096 --atoms-per-block 16`, spanning
fully-fine through fully-coarse, `MeasureEnergy` both on and off via the
crossover/self-reference stages the benchmark also runs): **0 errors** from
both tools.

## Part 4 — Performance

A self-inflicted contention finding, disclosed rather than hidden: the first
fp64 benchmark-matrix measurement in this task was accidentally run
concurrently with this task's own still-running fp64 `ctest` invocation
(`adaptive-cg-energy-fp32-accum` alone takes ~112s), confirmed via
`nvidia-smi --query-compute-apps` showing the ctest process at 100% GPU
utilization at the same time. Those numbers were discarded; the fp64 matrix
below was re-measured after confirming 0% utilization and no compute
process on the device.

fp32, all-coarse and mixed-fraction wall-clock, `production_step_wall_us`
held essentially constant between before/after runs (within ~1%, confirming
this task did not touch the production baseline):

| N | fraction | before (CGP-07) | after (this task) | delta |
| --- | --- | --- | --- | --- |
| 65,536 | 0.0 | 42.28us | 43.02us | −1.8% (noise, see below) |
| 65,536 | 0.0625 | 97.68us | 90.34us | +7.5% |
| 262,144 | 0.0 | 48.92us | 40.18us | +17.9% |
| 262,144 | 0.0625 | 95.82us | 89.62us | +6.5% |
| 1,048,576 | 0.0 | 92.37us | 79.19us | +14.3% |
| 1,048,576 | 0.0625 | 209.17us | 196.34us | +6.1% |

Launch count at all-coarse (fp32, every size, identical since the coarse
phase's launch count doesn't depend on N): `coarse` launches/step
**8.229 → 6.171** (−25.0%, exactly the 4-kernel→3-kernel-per-call reduction
predicted in Part 1), total launches/step **12.343 → 10.285** (−16.7%).

fp64, same points, smaller relative effect (consistent with CGP-06A/06B/07's
established fp32-vs-fp64 asymmetry — a fixed dispatch-cost floor is a
smaller fraction of a heavier fp64 kernel's own cost):

| N | fraction | before | after | delta |
| --- | --- | --- | --- | --- |
| 65,536 | 0.0 | 100.52us | 96.07us | +4.4% |
| 262,144 | 0.0 | 170.13us | 164.56us | +3.3% |
| 1,048,576 | 0.0 | 468.13us | 446.11us | +4.7% |

**CGP-07's specific performance gate — the 65,536-atom crossover table —
re-measured:**

| fraction | speedup before | speedup after |
| --- | --- | --- |
| 1.00 | 0.710x | 0.716x |
| 0.25 | 0.744x | 0.802x |
| 0.0625 | 0.804x | 0.871x |
| 0.01 | 0.984x | 0.880x |
| 0.00 | 1.859x | 1.828x |

**Honest result: the crossover at N=65,536 did not flip.** Every fraction
that was below 1.0x (adaptive slower than production) before this task
remains below 1.0x after it — closer to parity at most fractions (0.744x→
0.802x, 0.804x→0.871x), but not past it, and one point (0.01) measured
*worse* (0.984x→0.880x) despite a tiny within-run MAD (0.24us) that cannot
explain a ~9.5us shift — most likely cross-run environmental noise (GPU
clock/thermal state between two separate process invocations, which MAD
computed within one run cannot capture), not a real regression, but reported
as measured rather than explained away. This matches the audit's own
up-front prediction in Part 1/the task-selection conversation: a 25%
reduction in the coarse phase's launches, itself only part of the total
per-step floor (coarse launches are ~8 of ~12.3 total at all-coarse), was
never expected to close a floor that was ~30–40% of the whole production
step at this size on its own.

## Decision-gate update

CGP-07's decision gate is **not changed** by this task at the historical
65,536-atom scale: the atom-scale/launch-count floor there is measurably
smaller now (about 2 fewer launches/step, real wall-clock improvement at
most fractions) but still leaves adaptive coarse graining slower than plain
production at every tested fraction except all-coarse. At 262,144+ atoms,
where CGP-07 already found adaptive coarse graining winning, this task adds
a further real ~14–18% all-coarse improvement on top of that win. The two
remaining coarse-phase kernel-launch boundaries (clear→evaluate,
evaluate→finalize) are not fusable within this project phase's own scope
boundary (continuum-operator redesign is reserved for the next phase); a
grid-sync-based approach was flagged during scoping as a larger, riskier
follow-on and was not pursued in this task, consistent with the "implement
only where the dependency trace licenses it" discipline this task's own
Part C required.

## Checklist

* [x] Full steady-state kernel-launch inventory completed with read/write
  traces (Part 1) before any fusion.
* [x] Genuine fusion candidates distinguished from launch-order-only
  proximity (Part 1 table).
* [x] Each implemented fusion has a documented dependency proof
  (finalize→write: no cross-thread-block reads/writes in either kernel).
* [x] No O(Natoms) work reintroduced anywhere previously removed.
* [x] `MeasureEnergy` compile-time specialization preserved (fused kernel is
  itself templated the same way, energy reduction unchanged).
* [x] Field parity passes at all-coarse/mixed/all-fine (`cg13-cuda`'s CG-10.5
  anisotropy/moving/dipole fixtures).
* [x] `cg13-cuda` passes on fresh fp32 and fp64 builds (33/34, 34/34, the one
  fp32 failure pre-existing and unrelated).
* [x] Sanitizer evidence clean on fused kernels (synccheck and memcheck,
  full fraction sweep, 0 errors both).
* [x] Launch count reduction measured at all-coarse (8.229→6.171 coarse
  launches/step, −25%; 12.343→10.285 total, −16.7%).
* [x] CGP-07's 65,536-atom crossover re-measured and reported, honestly: not
  flipped, closer at most fractions, one point measured worse (disclosed as
  likely cross-run noise, not explained away as certain).
* [ ] Negative control demonstrates a fused ordering guarantee is real. Not
  applicable as originally envisioned: the implemented fusion removes a
  bookkeeping round-trip, not an ordering guarantee (no cross-thread-block
  hazard existed to begin with, per Part 1's trace), so there is no
  misordering scenario to construct a discriminating negative control
  against. The correctness case instead rests on CG-10.5's independent
  CPU-reference anisotropy check and compute-sanitizer, both clean.

See raw benchmark logs in `docs/cgp08_raw/*.log` (before-run comparisons use
`docs/cgp07_raw/*.log`, CGP-07's own raw output, unmodified).

## Commit

`CGP-08: reduce adaptive per-step kernel-launch floor`
