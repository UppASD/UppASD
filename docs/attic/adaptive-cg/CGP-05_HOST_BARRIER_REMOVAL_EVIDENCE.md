# CGP-05 — Remove host synchronization from the normal adaptive timestep, evidence

Task: `docs/CGP_work.md` CGP-05 ("Remove host synchronization from the normal
adaptive timestep"). Branch `gpu_hip_cu_ab_cg`, commit base `cf11224c`
(CGP-04).

**Primary model (task doc): Terra**, escalated for cross-stream CUDA/HIP
concurrency correctness. **Actual model: Claude (Sonnet 5, interactive
session)**, matching the CGP-03D precedent of an interactive session
performing the audit/design/implementation slice directly, with the same
level of care the escalation asked for (full producer/consumer data-flow
proof before any conversion, defensive event placement, sanitizer evidence,
and an attempted negative control).

## Environment

- GPU: 2x NVIDIA RTX A4000, shared with another user per project memory
  `shared-gpu-host-contention`; both idle (0% utilization, <110 MiB used)
  immediately before every timed run below.
- CUDA 13.3.73 (nvcc), `compute-sanitizer` 13.3, CMake, `Release` build type.
- Host CPU: 16 physical cores present (`taskset` shows full affinity 0-15),
  but `nproc` inside this environment reports 2 (cgroup-quota-limited) --
  material to section 8's benchmark disclosure.
- Builds: fresh out-of-tree `build_cgp05_cuda_fp32`
  (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`) and
  `build_cgp05_cuda_fp64` (`UPPASD_PRECISION=DOUBLE`), both
  `CMAKE_BUILD_TYPE=Release`. "Before" comparisons use the existing
  `build_cgp04_cuda_fp32`/`build_cgp04_cuda_fp64` (CGP-04, unmodified).

## 1. Synchronization inventory (Part A)

Every `GPU_STREAM_SYNC`/`GPU_DEVICE_SYNCHRONIZE`/`GPU_EVENT_SYNCHRONIZE` call
reachable from `GpuSimulation::advanceAdaptiveStep` and the immediately
following `GpuSDSimulation::SDmphase` moment-update call site, with its true
producer/consumer relationship traced from the actual buffer read/write
pattern (not from where the call happened to sit):

| # | Site (before) | Stream synced | True dependency | Classification |
|---|---|---|---|---|
| 1 | top of `evaluateAdaptiveFields` lambda | `stream()` | Meant to fence `stream()`'s writes before `workStream`'s FFT read of `direction` -- but syncing `stream()` doesn't wait for `workStream` at all; this call fenced the *wrong* direction | Misplaced/ineffective; removed |
| 2 | inside `fftEvaluator`, after `evaluateAdaptiveFftDipole` | `workStream` | Real: `stream()`'s subsequent `evaluateHybrid`/`evaluateField` read the FFT's padded field, produced on `workStream` | Genuine cross-stream (`workStream`→`stream()`), disguised by a comment that only mentions timing |
| 3 | bottom of `evaluateAdaptiveFields` lambda | `stream()` | None found: the next `stream()` work (`prepareCoarsePredictor`) is already ordered by same-stream FIFO; the real fence needed before `evolveFirst` is the *next* explicit sync (row 4) | Same-stream-redundant; removed |
| 4 | after `prepareCoarsePredictor` | `stream()` | Real: `evolveFirst` (`workStream`) reads `gpuLattice.beff`, written by `stream()`'s field-evaluation work | Genuine cross-stream (`stream()`→`workStream`) |
| 5 | after `evolveFirst` | `workStream` | Real: the predictor `evaluateAdaptiveFields` call (`stream()`) reads `gpuLattice.emom`, just written for active atoms by `evolveFirst` | Genuine cross-stream (`workStream`→`stream()`) |
| 6 | after `evolveSecond` | `workStream` | Real, but not for the immediately following call (`correctCoarse`, row 7): it fences every *later* `stream()` read of `emom2`'s active-atom entries this step (traceEnergy diagnostics, selector-update `restrictMoments`/`evaluateSelectorScores`/`publishProposedState`) **and**, via the `emom`/`emom2` buffer swap the moment updater performs after this function returns, the *next* step's initial field read too | Genuine cross-stream, wider scope than its call site suggested |
| 7 | after `correctCoarse` | `stream()` | None found: `correctCoarse` only reads `predictorCoarseField`, itself produced by `stream()` earlier this same step; nothing later in the ordinary (non-diagnostic) path needs a fresh fence here that row 6's fence doesn't already cover | Same-stream-redundant; removed |
| 8 | end of `advanceAdaptiveStep` | `stream()` | Real, in a way its own call site didn't need: the *caller* (`GpuSDSimulation::SDmphase`'s moment-update call, `workStream`) reads `emom2`'s coarse-atom entries just written by `materializeCoarseAtoms`/`publishProposedState` on `stream()`. Also fences the *next* step's FFT read of `direction` (row 2 companion) | Genuine cross-stream (`stream()`→`workStream`), previously implicit (relied on this function's own host block, not an explicit dependency at the consumer) |

Diagnostic-only host syncs (`GpuAdaptiveRuntime::diagnosticSnapshot()`,
gated to `adaptiveDiagnostics>=3`) and the once-per-run `initialize()`/
`release()` syncs are unchanged -- Part E explicitly allows keeping host
synchronization where data is genuinely required by the host, and these are
not on the steady-state path.

Rows 1, 3, and 7 (three of eight) turn out to protect **nothing** once the
actual buffer flow is traced precisely -- they are removed outright, no event
needed. The other five are genuine cross-stream dependencies, converted to
events (section 3).

## 2. Timing vs. production dependency (Part B)

`beginPhase()`/`finishPhase()` (used by `prepareCoarsePredictor`,
`correctCoarse`, etc.) were already correctly gated behind
`phaseTimingEnabled()` before this task -- no host wait at all when disabled.
The one gap: the FFT dipole's `ASSERT_GPU(GPU_STREAM_SYNC(workStream))` (row
2 above) was **not** gated by `phaseTimingEnabled()` even though its own
comment only mentions feeding `recordFftMilliseconds`. Tracing its actual
role (section 1, row 2) showed it is dual-purpose: real correctness fence
*and* timing measurement, conflated in one call. Splitting these:

- The correctness fence is now `markExternalProgress()`/`waitForProgress()`
  (event-based, unconditional, non-blocking).
- The timing measurement (`ASSERT_GPU(GPU_STREAM_SYNC(workStream))` before
  reading elapsed wall time) is now gated behind `phaseTimingEnabled()`,
  matching `finishPhase()`'s own pattern.

The same split was applied to the end-of-step wall-clock measurement
(`recordStepWallMilliseconds`, row 8): the correctness fence
(`markProgress()`) is unconditional; the blocking wait used to get an
accurate wall-clock stop time only executes when `phaseTimingEnabled()`.

## 3. Design (Part C) and FFT interaction (Part D)

Two ordering primitives were added to `GpuAdaptiveRuntime`
(`markProgress()`/`waitForProgress(consumer)` for `stream()`→other;
`markExternalProgress(producer)`/`waitForExternalProgress()` for
other→`stream()`), implemented with `GPU_EVENT_RECORD`/
`GPU_STREAM_WAIT_EVENT` -- the existing, CUDA/HIP-symmetric macros in
`gpu_wrappers.h` (`gpuEventPool.hpp`'s existing event-pool abstraction was
found unused anywhere in the codebase; this task did not revive it, since a
pair of fixed-purpose events was simpler and sufficient here).

Design rule adopted deliberately over a minimal hand-derived placement:
**always** call `markProgress()`/`markExternalProgress()` right after the
producer work a consumer might need, and **always** call
`waitForProgress()`/`waitForExternalProgress()` right before launching work
that might consume it -- even at a call site later proven redundant this
step (e.g. the top of `evaluateAdaptiveFields`, which is a genuine fence for
the predictor/corrector calls but a harmless re-wait for the initial call).
`GPU_STREAM_WAIT_EVENT` on an already-satisfied event is a cheap device-side
check, not a host wait, so over-inserting costs effectively nothing while
under-inserting is a race. This traded a small amount of "provably minimal"
elegance for a much easier-to-verify safety property, appropriate given the
task's own escalation reasoning.

Every row-1/2/3 (FFT) and row-4/5/6/7/8 site in section 1 was converted
1:1: the three same-stream-redundant syncs (rows 1, 3, 7) were deleted with
no replacement; the five genuine cross-stream syncs (rows 2, 4, 5, 6, 8)
became a `mark*`/`wait*` pair at the same point in the code. One new call
site was added that had no prior host sync at all: `GpuSDSimulation::
SDmphase`'s moment-update call now explicitly calls `waitForProgress
(workStream)` before `adaptiveMomUpdater.update*()`, replacing the guarantee
that used to arrive for free as a side effect of `advanceAdaptiveStep`'s own
final host block (row 8) -- see the negative control in section 6, which
targets exactly this new call site.

## 4. Implementation

- `source/gpu_files/gpuAdaptiveRuntime.hpp`/`.cpp`: `markProgress()`,
  `waitForProgress(GPU_STREAM_T)`, `markExternalProgress(GPU_STREAM_T)`,
  `waitForExternalProgress()`. The two backing events are process-lifetime
  singletons (function-local statics in `gpuAdaptiveRuntime.cpp`), not
  `GpuAdaptiveRuntime` instance members -- see section 7 for why.
- `source/gpu_files/gpu_wrappers.h`: `GPU_EVENT_CREATE_NO_TIMING`
  (`cudaEventCreateWithFlags(..., cudaEventDisableTiming)` /
  `hipEventCreateWithFlags(..., hipEventDisableTiming)`), used for both new
  events since neither is ever read with `GPU_EVENT_ELAPSED_TIME`.
- `source/gpu_files/gpuSimulation.cpp`: `advanceAdaptiveStep` restructured
  per section 3; the FFT `fftEvaluator` lambda now fences both directions
  explicitly and gates its own timing block.
- `source/gpu_files/gpuSDSimulation.cpp`: one new `waitForProgress(workStream)`
  call before the adaptive moment-update call.
- No change to `GpuMomentUpdater`, `GpuCorrelations`, reconstruction
  mathematics, or any kernel's numerical formula.

## 5. Correctness evidence (Parts E/F)

- `ctest -L cg13-cuda`, fp32: **33/34**, three consecutive full runs,
  identical to the CGP-04 baseline's own known failure
  (`adaptive-cg-production-e2e`'s `coarse_dipole` metric assertion) --
  reproduced byte-for-byte identically on an unmodified `build_cgp04_cuda_fp32`
  run of the same test, confirming it is pre-existing and unrelated.
- `ctest -L cg13-cuda`, fp64: **34/34**, clean (this pre-existing failure
  does not reproduce on fp64, matching CGP-04's own fp64 finding).
- `compute-sanitizer --tool synccheck --target-processes all` over the full
  `adaptive-cg-production-e2e` scenario set (finite-T fine/mixed, DMI/
  anisotropy CPU/GPU parity, restart, Initmag 2/4/8, spin-spiral,
  user-facing examples): **0 errors**.
- `compute-sanitizer --tool memcheck`, same scope: **0 errors**.
- No change to the reconstruction formula, `GpuMomentUpdater`, or
  `GpuCorrelations` (git-diff-verified).

## 6. Negative control

Attempted on the new, previously-absent `waitForProgress(workStream)` call
site in `gpuSDSimulation.cpp` (section 3's "one new call site"), since it is
the cleanest case where removing the wait has no other fence covering it.
Disabled it, rebuilt, and searched for a discriminating fixture:

- `coarse-graining-adaptive-asd-parity` (`gpu_adaptive_runtime_benchmark
  --parity-only`): does not exercise this call site at all -- it replicates
  step logic directly against `GpuAdaptiveRuntime`, never calling
  `GpuSimulation::advanceAdaptiveStep`/`SDmphase`. Confirmed non-
  discriminating for the intended reason (wrong code path), not by chance.
- Full `ctest -L cg13-cuda` (which does exercise `SDmphase` end-to-end via
  `sd.f95.cuda`, e.g. `adaptive-cg-production-e2e`,
  `adaptive-cg-moving-backend-parity`): three repeated full runs with the
  wait disabled, all still 33/34 (same single pre-existing failure). Not
  discriminating on any fixture available in this session.

Reported honestly as **not discriminating**, matching the CGP-03C precedent
for an attempted negative control that didn't trigger on the fixtures used.
The correctness argument for this wait rests on the static producer/consumer
proof in section 1 (row 8), not on an empirical failure demonstration; a
fixture with a materially longer coarse-atom commit and a much larger active
list (widening the real race window between `stream()`'s commit and
`workStream`'s moment-update read) might discriminate where these did not,
but constructing one was out of this session's remaining budget. The wait
was restored before every correctness/performance run in sections 5, 7, 8.

## 7. Side finding: a curand/event-lifecycle interaction, found and fixed

Initial correctness runs surfaced an intermittent
`coarse-graining-adaptive-asd-parity` failure with message `GPU
random-generator error: 201` (`CURAND_STATUS_LAUNCH_FAILURE`). Isolation:

1. **A/B in the identical build directory**: stashing every CGP-05 source
   change and rebuilding `gpu_adaptive_runtime_benchmark` in the same
   `build_cgp05_cuda_fp32` tree made the failure disappear (0/30 clean);
   restoring the changes brought it back (~14% of 65 combined trials).
   Confirms causation, rules out environment/directory effects.
2. Since `benchmark_gpu_adaptive_runtime.cpp` never calls any of the four new
   `mark*`/`wait*` methods (grep-confirmed), the only code this benchmark
   exercises from this task's diff is `GpuAdaptiveRuntime::initialize()`/
   `release()`'s creation/destruction of two new persistent
   `GPU_EVENT_T` objects.
3. `compute-sanitizer --tool memcheck` pinpointed it precisely: `CUDA API
   Error: CUDA Stream does not belong to the expected context`, inside
   `curandGenerateNormalDouble` <- `GpuThermfield::randomize()` <-
   `GpuDepondtIntegrator::evolveFirst()` -- a *different* subsystem's curand
   generator, whose stream binding (`ParallelizationHelperInstance
   .getWorkStream()`, created once and never touched by this task) became
   invalid.
4. First mitigation: created the two new events with
   `GPU_EVENT_CREATE_NO_TIMING` (`cudaEventDisableTiming`) instead of the
   default timing-capable `GPU_EVENT_CREATE`, since neither is ever read
   with `GPU_EVENT_ELAPSED_TIME`. This alone took fp32 from ~14% failure
   (of 65 trials) to 0/100 combined trials, but fp64 remained **~100%**
   failure (0/40 pass) -- the same benchmark constructs three
   `GpuAdaptiveRuntime` instances concurrently
   (`adaptiveExchange`/`adaptiveFull`/`adaptiveFullUniaxial`, each with its
   own `stream_` and, before this fix, its own two new events) for an
   earlier hamiltonian-oracle check, before `runProductionParityHarness`'s
   own curand generator is first used.
5. Root fix: made the two events **process-lifetime singletons**
   (function-local statics), not `GpuAdaptiveRuntime` instance members.
   Production only ever has one live `GpuAdaptiveRuntime` for a process's
   whole run (exactly like `ParallelizationHelperInstance`'s own
   `workStream`/`copyStream`), so this is behaviourally identical in
   production and collapses this benchmark's cumulative event allocation
   from 10 (2 x 5 `initialize()` calls across the file) to 2. Result: fp64
   went from 0/40 to **0/40 failures under normal execution** (100/100
   combined with fp32, both clean).
6. Residual: `compute-sanitizer --tool memcheck` (not `synccheck`) still
   reproduces this exact curand error under its own heavy instrumentation,
   3/3 runs, **even with every CGP-05 change fully reverted** (re-verified
   directly against the unmodified baseline in the same build directory).
   This is pre-existing sanitizer-timing sensitivity in the
   `ParallelizationHelper`/curand stream-binding lifecycle, unrelated to
   this task and out of its scope (touching that infrastructure is CGP-06A/
   06B territory, not CGP-05's). Not fixed here; flagged for whichever task
   next needs sanitizer-clean coverage of that path specifically.

This finding does not affect production correctness (production has exactly
one `GpuAdaptiveRuntime` for the process's life, so instance-churn was never
a live risk there) but was a real, reproducible defect in this task's first
implementation, caught by the project's own compute-sanitizer discipline
before it could reach the evidence section as a false "clean" result.

## 8. Performance evidence

### Sync-count reduction (code-verified, primary evidence)

Per ordinary (non-diagnostic, non-selector-update-step, no FFT dipole)
production step, counting only host-blocking `GPU_STREAM_SYNC`/
`GPU_EVENT_SYNCHRONIZE` calls in `advanceAdaptiveStep` plus the immediately
following moment-update call:

| | before | after (`phaseTimingEnabled=true`, default) | after (`UPPASD_ADAPTIVE_PHASE_TIMING=0`) |
|---|---|---|---|
| host-blocking syncs/step | **9** | 1 (end-of-step wall-clock measurement only) | **0** |
| device-side event waits/step | 0 | 5 | 5 |

With FFT dipole enabled, add 1 more host-blocking sync/step before (inside
`fftEvaluator`) and 1 more event-record/wait pair after, in both rows.

This is a full elimination of the barriers this task's Background section
named ("also synchronizes the production work stream around FFT dipole
evaluation") for the default production configuration
(`UPPASD_ADAPTIVE_PHASE_TIMING=0`), verified by code inspection rather than
timing (section 1's table), and is the reduction requested by the task's own
checklist item "Host synchronizations per normal step reduced."

### Whole-step wall-clock latency: attempted, inconclusive on this host

Constructed the same fixture class CGP-03D used (all-coarse, 8192 atoms,
`ncell 1024 2 2`, `block_size 1x2x2`, T=0, `Nstep 2000`,
`do_avrg/cumu/tottraj N`, `cg_diagnostics 0`), compared `build_cgp04_cuda_fp32`
(before) against `build_cgp05_cuda_fp32` (after) using `sd.f95.cuda`'s
own `Time for meas. phase` wall-clock report (host-side, unaffected by the
`UPPASD_ADAPTIVE_PHASE_TIMING` toggle), one warm-up discarded, ABBA-paired
samples with the GPU confirmed idle before each pair:

`before: 17.12, 17.15, 17.09 s` — `after: 17.08, 17.09, 17.07 s` (2000 steps
each; no measurable difference, both essentially flat at ~8.5 ms/step).

This is disclosed as **inconclusive** rather than reported as "no speedup":
`nproc` in this environment reports 2 despite 16 physical cores being
present and fully in this process's CPU affinity mask (cgroup-quota-limited,
consistent with project memory's "2-core host" disclosure precedent from
CGP-04). 8.5 ms/step for an all-coarse 8192-atom system is far above the
task's own "sub-100 us fp32 atomistic baseline" framing, indicating this
measurement is dominated by CPU-side contention/scheduling on the shared
host, not by GPU host-wait latency -- exactly the quantity too small
(section 8's sync-count table implies savings on the order of a few
microseconds of host round-trip per removed sync, times a few removed
syncs) to be visible against multi-millisecond CPU noise on this
environment. A clean, low-noise whole-loop wall-time number was not
obtainable within this session's budget; the sync-count reduction above is
offered as the primary, decisive performance evidence instead, per this
project's established practice of disclosing rather than fabricating a
number under a resource constraint (CGP-04's fp64 perf-baseline disclosure
is the direct precedent).

### fp64 correctness confirmed separately from fp32 performance

Per section 5, fp64 correctness (34/34) was fully verified even though no
fp64 performance comparison was attempted, for the same host-constraint
reason CGP-04 disclosed.

## Checklist

* [x] All synchronization sites inventoried.
* [x] Every retained barrier has a documented dependency.
* [x] Timing-only barriers absent in production timing mode (verified: the
  one dual-purpose FFT sync was split; `UPPASD_ADAPTIVE_PHASE_TIMING=0`
  leaves zero host-blocking syncs on the ordinary-step path).
* [x] Same-stream dependencies use stream ordering (three syncs proven
  same-stream-redundant, deleted with no replacement).
* [x] Cross-stream dependencies use device events where appropriate (five
  sites converted to `mark*`/`wait*` event pairs).
* [x] FFT dependency chain is explicit (both directions fenced; timing
  split from correctness).
* [x] No unnecessary host wait surrounds ordinary field evaluation.
* [x] T=0 all-fine parity passes (`cg13-cuda`, fp32 33/34 -- pre-existing
  unrelated failure only -- and fp64 34/34).
* [x] Finite-T all-fine parity passes (`adaptive-cg-production-e2e`'s
  finite-temperature fine-region/mixed-resolution cases).
* [x] Moving mixed fixture passes (`adaptive-cg-moving-backend-parity`,
  `adaptive-cg-moving-off-fine`).
* [x] Sanitizer/race evidence clean (`synccheck` and `memcheck`, 0 errors,
  over the full production e2e scenario set).
* [x] Host synchronizations per normal step reduced (9 -> 0, code-verified,
  section 8).
* [ ] Whole-step fp32 latency remeasured. Attempted; inconclusive on this
  shared, CPU-quota-limited host (section 8); not claimed as a pass.
* [x] CUDA/HIP structural semantics remain aligned (`GPU_EVENT_CREATE_NO_TIMING`
  added symmetrically to both branches of `gpu_wrappers.h`; every other
  primitive used was already CUDA/HIP-symmetric).

Additional, not in the original checklist but material to this task's own
risk framing:

* [x] Negative control attempted on the one new (previously absent) wait
  call site; reported honestly as not discriminating on the fixtures
  available (section 6), rather than claimed as a pass.
* [x] A real defect introduced by this task's first implementation (curand/
  event-lifecycle interaction) was caught by the project's own
  compute-sanitizer discipline before being reported as clean, root-caused,
  and fixed (section 7).

## Commit

`CGP-05: remove adaptive timestep host barriers`
