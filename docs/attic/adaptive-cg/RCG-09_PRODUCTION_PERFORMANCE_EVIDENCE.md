# RCG-09: honest production performance result

**Task:** RCG-09, "Establish an honest production performance result"
**Parent:** `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` §9
**Date:** 2026-08-14
**Status:** Measurement complete for CUDA/fp64 on this host. **No crossover
exists and none is claimed.** HIP unavailable.

**Superseded by RCG-09A.3:** the adaptive atomistic integration described in
this historical performance report was the transitional normalized-Heun path.
It has been deleted. Current adaptive fine spins use the shared production
Depondt predictor/corrector and production thermal-field generation; the
timings below are retained only as historical pre-migration evidence.

---

## 1. Summary

RCG-09 asked whether adaptive coarse graining is faster than UppASD's real
atomistic GPU path under identical physics. Measured against
`GpuHamiltonianCalculations` and `GpuDepondtIntegrator` — the same objects the
feature-off timestep loop calls — the answer on CUDA/fp64 is **no, by a wide
margin, and the gap grows with system size.**

At 16 384 atoms, one adaptive step costs **11.9 s** against **0.35 ms** for the
production atomistic step: about **3.4 × 10⁴ times slower**, not faster. There
is no fine fraction at which the adaptive path becomes cheaper, because the
dominant cost does not scale with the fine fraction at all.

The cause is isolated to one line. `evaluateAdaptiveAtomisticBonds` accumulates
the bilinear energy with `atomicAddEnergyTerm` — a compare-and-swap loop on a
single global `double` — and every bond thread calls it unconditionally. With
~10^5 bond threads contending on one address, the kernel serializes, and its
cost grows as `O(bonds^1.8)` and reaches `O(bonds^2.0)` at the largest size
measured, instead of linearly. A negative control that disables only that
accumulation makes the same kernel **10 489× faster** at 114 688 bonds and
restores sub-linear scaling.

This supersedes the `1.30x` claim in the parent blueprint's CG-10 evidence.
That number compared the adaptive runtime at a low fine fraction against the
adaptive runtime at 100% fine (finding F-10); it was never a comparison with
UppASD. Under the blueprint's own §13, item 9 fails, and the result may be
accepted only as a correctness/reference implementation.

---

## 2. What was measured, and against what

Two independent harnesses, deliberately not sharing a code path.

### 2.1 `PERF-ATOMISTIC-PROD` and `PERF-CG-SWEEP`, in process

`tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp`, rewritten for this
task. It builds one geometry and drives two things on it:

* the production path: `GpuHamiltonianCalculations::initiate` on neighbour
  lists in the layout a real run uploads, then the exact feature-off step
  sequence of `gpuSDSimulation.cpp` — `heisge`, `evolveFirst`, `heisge`,
  `evolveSecond`;
* the adaptive path: `GpuAdaptiveRuntime::integrateHeun` on the same atoms,
  bonds, moments and exchange constant.

**The comparison is licensed by a measurement, not by an assertion.** Before any
timing, both paths evaluate the field on the same state with every block fine,
and the harness requires them to agree. They agree to `2.8e-16` relative at
fp64 — round-off, not tolerance. If that check fails the harness reports no
speedup and exits non-zero.

Two limitations of this harness, recorded rather than hidden. First, its
absolute timings are not resolved at these sizes: §3.5 shows a 13× spread across
three samples of identical work, and lean-mode and instrumented-mode
measurements of the same kernels differing by 8×. Use it for scaling and
attribution, not for absolute cost. Second, its *production* arm is launch-bound
at the sizes measured. The reported production step wall
time barely varies with system size (55 µs at 1 024 atoms, 371 µs at 65 536),
and the harness's own `rng_and_integration_residual_us` comes out **negative** —
twice the isolated `heisge` time exceeds the whole four-call step, which is only
possible if per-call launch latency dominates the isolated measurement and
pipelines away inside the step loop. So the in-process production figures are a
sanity check and a scaling probe, **not** the production absolute baseline. The
absolute baseline is §3.1's end-to-end two-point measurement, which has no such
problem because it differences two runs of 1000 real steps.

### 2.2 End-to-end, the ordinary executable

`tests/coarse_graining/run_rcg09_perf_e2e.py` with the fixtures in
`tests/coarse_graining/perf/`. It runs the normal `sd` binary on two inputs
that differ only in whether adaptive coarse graining is enabled. No internal
staging, no direct kernel invocation. Setup is separated from steady state by
running each arm at a short and a long `Nstep` and dividing the difference —
process start, context creation, input parsing, neighbour-list construction and
adaptive setup cancel exactly.

### 2.3 The synthetic control, renamed and separated

The fused-multiply-add loop that finding F-10 objected to now lives alone in
`tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp`, built as
`gpu_adaptive_inactive_overhead_microbenchmark`. Its own output declares
`not_a_production_baseline=true`. It answers only whether an inactive adaptive
runtime costs anything; it contains no UppASD physics and cannot be read as a
baseline by anyone reading the artifact.

---

## 3. The result

### 3.1 End to end, 16 384 atoms, spin spiral — the headline

Ordinary `sd.f95.cuda`, `ncell 32 16 16` bcc with a two-atom basis and four
exchange shells, 1024 blocks of 16 atoms, fp64, one ensemble, `temp 0`,
`SDEalgh 1`. Steady state by two-point subtraction; raw output in
`docs/rcg09/e2e_two_point.txt`.

Adaptive arm, from the runtime's own step wall clock at `Nstep` 1 and 3:

| Quantity | `Nstep` 1 | `Nstep` 3 | Per step |
| --- | --- | --- | --- |
| adaptive-CG step wall | 11 757.8 ms | 35 614.5 ms | **11 928 ms** |
| of which atomistic phase | 11 740.8 ms | 35 582.0 ms | 11 921 ms (99.9%) |

Feature-off arm, three interleaved rounds at `Nstep` 100 and 1100 — a 1000-step
difference, because at `Nstep` 3 the whole run is 0.44 s and the per-step cost
is below the resolution of process wall time:

| Round | `Nstep` 100 | `Nstep` 1100 | Per step |
| --- | --- | --- | --- |
| 1 | 0.53 s | 0.86 s | 0.33 ms |
| 2 | 0.55 s | 0.91 s | 0.36 ms |
| 3 | 0.52 s | 0.87 s | 0.35 ms |

Median **0.35 ms/step**, MAD 0.01 ms.

**11 928 ms against 0.35 ms: the adaptive path is ~34 000× slower.**

An earlier draft of this document quoted ~60 ms for the feature-off arm and a
170× ratio. That figure came from UppASD's own "Time for one meas. iter" at
`Nstep` 5, which amortizes ~0.3 s of setup over five steps and is therefore
almost entirely setup. The two-point figures above are the correct ones and the
gap is two orders of magnitude larger than that draft claimed.

Every other adaptive phase together — coarse field, interface, selector,
polarization, compaction, integration — is under 1.6 ms per step. Adaptive setup
is ~38.9 s, separately (see RCG-09-FU2).

The selector left all 1024 blocks fine on this texture
(`resolution_counts fine=1024 interface=0 coarse=0`), so this is the all-fine
case. That does not rescue the comparison: §3.2 and §3.4 show the dominant cost
is insensitive to the fine fraction.

### 3.2 Atomistic phase scaling, all blocks fine

`gpu_adaptive_runtime_benchmark` on `build_rcg09_cuda_fp64`, 4 atoms per block,
fp64, all blocks fine. Median per-phase device time for the atomistic phase of
one Heun step. Field parity passes at every size (2.83e-16).

| Blocks | Atoms | Bonds | Atomistic phase | Growth per 4× bonds |
| --- | --- | --- | --- | --- |
| 64 | 256 | 448 | 512 µs | — |
| 256 | 1 024 | 1 792 | 1 522 µs | 3.0× |
| 1 024 | 4 096 | 7 168 | 15 374 µs | 10.1× |
| 4 096 | 16 384 | 28 672 | 176 233 µs | 11.5× |
| 16 384 | 65 536 | 114 688 | 2 902 048 µs | 16.5× |

Linear work would give 4× per 4× bonds. The measured growth is 10–16×, i.e.
`O(bonds^1.76)` then `O(bonds^2.02)` over the last two intervals: the cost is
becoming quadratic in bond count. A single all-fine field evaluation at
114 688 bonds takes **2.9 seconds**.

### 3.3 Negative control: the cause

Same binary, same geometry, one line disabled —
`atomicAddEnergyTerm(&kernels.energyTerms[0], ...)` in
`evaluateAdaptiveAtomisticBonds`, guarded behind
`RCG09_DISABLE_BOND_ENERGY_ATOMIC`. Nothing else changed. The field parity check
still passes, because the field never depended on that accumulation.

Both builds are `CMAKE_BUILD_TYPE=Release` with identical
`CMAKE_CUDA_FLAGS_RELEASE=-O3 -DNDEBUG`; the control build adds only
`-DRCG09_DISABLE_BOND_ENERGY_ATOMIC` to `CMAKE_CUDA_FLAGS`.

| Bonds | With energy atomic | Without | Ratio |
| --- | --- | --- | --- |
| 7 168 | 15 374 µs | 63.5 µs | **242×** |
| 28 672 | 176 233 µs | 100.0 µs | **1 762×** |
| 114 688 | 2 902 048 µs | 276.7 µs | **10 489×** |

Without it, 4× bonds costs 1.6–2.8× time: sub-linear, launch- and
bandwidth-bound, as a correctly parallelized bond kernel should be. With it, the
kernel is a queue in front of one 8-byte address.

This is causation, not correlation: one line, toggled, four orders of magnitude,
with the scaling exponent falling from ~2.0 to ~0.7 alongside it.

The patch was reverted before any committed build. It exists in this document
as evidence, not in the tree as a feature — skipping the energy accumulation
silently would break the energy contract RCG-06B established.

### 3.4 Why the fine fraction cannot rescue it

The bond kernel launches over *all* bonds and returns early only when neither
endpoint is atomistically owned. Coarsening reduces the number of threads that
do work, but the surviving threads still serialize on the same single address,
and the launch itself is still sized by the total bond count (RCG-08 §7.4).
Measured live-bond fractions fall roughly with the fine fraction — 100%, 79%,
54%, 29%, 17% across the sweep — while the atomistic phase falls far more
slowly, because contention on one address is not proportional to the number of
contenders doing useful work.

The in-process harness reports `production-crossover result=NOT_OBSERVED` at
every system size measured, which is the machine-readable form of the same
conclusion.

### 3.5 What the negative-control build shows about the design itself

The control build of §3.3 is the closest available proxy for "after the
accumulator is fixed", so its fine-fraction sweep is the most informative data
in this document. At 65 536 atoms / 114 688 bonds.

**Every column below is a phase *inside the adaptive step*.** None of them is
the production baseline. In particular the "atomistic phase" column is the
adaptive runtime's own short-range Hamiltonian evaluation over the atoms that
are still atomistically owned — it varies with the fine fraction precisely
because that is what coarse graining is for. The production atomistic baseline
is a single constant with no fine-fraction dependence, shown separately beneath
the table.

Per-phase device time within one adaptive Heun step:

| Fine fraction | Live bonds | Atomistic phase | Coarse phase | Interface phase | Integration phase |
| --- | --- | --- | --- | --- | --- |
| 1.000 | 100.0% | 276.7 µs | 57.3 µs | 56.0 µs | 165.4 µs |
| 0.750 | 75.0% | 196.7 µs | 14 500 µs | 231.9 µs | 137.0 µs |
| 0.500 | 50.0% | 158.4 µs | 42 923 µs | 346.4 µs | 129.0 µs |
| 0.250 | 25.0% | 120.4 µs | 83 918 µs | 355.8 µs | 79.5 µs |
| 0.125 | 12.5% | 85.7 µs | 107 303 µs | 389.6 µs | 105.2 µs |
| 0.0625 | 6.3% | 81.0 µs | 116 339 µs | 404.8 µs | 97.6 µs |
| 0.000 | 0.0% | 45.9 µs | 131 401 µs | 472.1 µs | 66.5 µs |

Reference line, *not* a column above: the production feature-off step on the
same geometry costs **1 170.5 µs**, the same at every row because the feature-off
path has no fine fraction.

**No factor between that number and the all-fine adaptive step is claimed here,
because this harness cannot support one at this size.** The all-fine adaptive
step's three raw lean-mode samples are 5 907 / 4 703 / 453 µs — a 13× spread,
MAD 1 204 µs on a median of 4 703 µs. The same work with phase timing enabled
measures 587.8 µs with MAD 14.9 µs. Those two figures differ by 8× for identical
kernels in identical order, and the production baseline is no better: its
samples are 1 171 / 1 176 / 719 µs and its `rng_and_integration_residual_us`
is −450 µs, i.e. two isolated `heisge` calls appear to cost more than the whole
four-call step containing them. Three repetitions of three iterations on a
contended device is simply not enough to resolve sub-millisecond differences.

The trustworthy comparison is §3.1's end-to-end two-point measurement, which
differences two 1 000-step runs and has none of these problems.

### 3.5.1 All-fine adaptive is not "atomistic plus overhead"

A natural expectation is that at 100% fine the adaptive path should fall back to
the ordinary atomistic computation and therefore cost slightly *more* than
feature-off. It should not, and the reason matters for interpreting any future
clean-device number:

* **Pair arithmetic is halved.** `heisge` walks a per-atom neighbour list and
  evaluates every pair twice, once from each endpoint. The adaptive path walks a
  compact unique-bond list and evaluates each pair once, scattering to both
  endpoints with atomics. Same physics — this is exactly what the `2.8e-16`
  parity check confirms — but roughly half the multiply-adds, traded for atomic
  traffic.
* **The integrators differ.** Feature-off runs Depondt (`SDEalgh 1`), which
  builds an effective field, performs a Rodrigues rotation, and calls
  `thermfield.randomize()` every step — 196 608 random numbers per step at this
  size, even at `temp 0`. The adaptive path runs a deterministic Heun update and
  generates none.
* **Launch structure differs.** Feature-off issues a handful of large kernels;
  the adaptive step issues roughly thirty, several sized by upper bound rather
  than actual work (RCG-08 §7.4).

So all-fine adaptive is a genuinely different implementation of the same
physics, not a fallback path, and it is not structurally required to be slower.
Whether it is faster or slower on a clean device is an open question this
harness has not answered.

Both differences were noted here to explain a timing observation, but both are
really physics questions — unique-pair scatter is uncontroversial for symmetric
exchange and not obviously sound for antisymmetric DMI, and an atomistic region
at finite temperature should carry a stochastic field. They are recorded for an
independent reviewer in
`docs/RCG-09_ADAPTIVE_ATOMISTIC_KERNEL_QUESTIONS.md`. RCG-09 did not
investigate or act on either; note in particular that the `2.8e-16` parity check
above runs on an isotropic-exchange fixture with `do_dm` disabled, so it is not
evidence about the DMI case.

Two things follow, and they point in opposite directions.

**The coarse-graining mechanism works.** The adaptive step's atomistic phase
falls from 276.7 µs to 45.9 µs as blocks coarsen, tracking the live-bond
fraction. That is the
intended effect and it is real, measured, and monotone. It also exposes the
ceiling: the reduction saturates at **6.0×**, and the residual 45.9 µs at zero
fine blocks is the full-system `O(N)` per-step work RCG-08-FU4 catalogued —
state saves, `clearAdaptiveAtomistic`, ghost prolongation/restriction, work-scan
initialization — which no amount of coarsening removes. **6× is therefore the
measured upper bound on what coarsening can save in the atomistic phase**, not a
number to be extrapolated past.

**The coarse operator has the same defect, mirrored.** The coarse phase grows
from 57 µs to 131 401 µs — a factor of 2 300 — in the direction where it is
supposed to be getting cheaper. `evaluateAdaptiveCoarseTensor` and
`finalizeAdaptiveCoarseLocal` accumulate `energyTerms[2]` through
`energyTerms[5]` through the same single-address FP64 CAS loop, contended by
every active block instead of every bond.

So the current implementation has **no operating point that wins**: all-fine is
bounded by the bond accumulator, all-coarse by the coarse accumulator, and every
mixed point pays both. That is one defect class in two places, not two unrelated
problems, and it is why RCG-09-FU1 is scoped below as a pattern to eliminate
rather than a line to change.

---

## 4. Measurement conditions, and what they do and do not support

### 4.1 Device state

**Every number in §3 was taken on a contended, thermally throttled device.**
Both RTX A4000s were being driven by another user's process at 63–75%
utilization, 92–93 °C, with `clocks_event_reasons.active=0x20` (thermal
slowdown) and ~4 GB of foreign allocation on each. This is exactly the condition
RCG-08-FU3 recorded and required RCG-09 not to repeat.

What that permits and forbids:

* **Permitted:** the §3.1/§3.2/§3.3 conclusions. A 10 489× ratio, a scaling
  exponent measured across a 256× range of bond counts, and a 3.4 × 10⁴
  end-to-end gap are not products of a busy neighbour. Contention on this host has been observed at the
  tens-of-percent level, not the four-orders-of-magnitude level.
* **Forbidden:** treating any absolute microsecond figure here as the clean-device
  cost of anything. They are upper bounds of unknown tightness.

`run_rcg09_perf_e2e.py` refuses to run on a device in this state and exits 4
unless `--allow-dirty-device` is passed; the numbers above were taken with the
in-process harness and with that flag, and are labelled diagnostic throughout.
**A clean-device re-measurement is required before any of these figures is cited
as a final quantity** — but not before the qualitative conclusion, which no
plausible clean-device correction can reverse.

### 4.2 Environment

| Item | Value |
| --- | --- |
| Device | NVIDIA RTX A4000, driver 610.57.04 |
| SM clock during measurement | 1125–1320 MHz of 1560 MHz max, thermally capped |
| Backend / precision | CUDA, fp64 (`UPPASD_PRECISION` default) |
| Build | fresh out-of-tree `build_rcg09_cuda_fp64`, `CMAKE_BUILD_TYPE=Release` |
| CPU threads | OpenMP reports 2 usable threads of 16 present |
| HIP | **unavailable** — no toolchain and no device, as in RCG-02..RCG-08 |

No result in this document is extrapolated to HIP. The CUDA and HIP kernels
share launch geometry through the `ADAPTIVE_LAUNCH` macro, which makes the same
defect structurally likely on HIP, but that is an argument and not a
measurement.

### 4.3 The instrument was removed from the measurement

RCG-08-FU2 recorded that per-phase device timing costs ten blocking host
synchronizations per step, ~38% of step wall time. Timing the adaptive path with
that enabled measures the instrument. This task added
`GpuAdaptiveRuntime::setPhaseTimingEnabled` — default **on**, so RCG-06C's
accepted per-phase evidence and every existing caller are unaffected — and the
harnesses take headline wall times with it off and the phase breakdown with it
on, reporting both plus the difference. Production runs can select it with
`UPPASD_ADAPTIVE_PHASE_TIMING=0`, and the run prints which mode it used.

This changed nothing about the conclusion. At a 12 s step, a 38% instrumentation
overhead is not what stands between adaptive coarse graining and a crossover.

### 4.4 One asymmetry, disclosed rather than normalized

`GpuDepondtIntegrator::evolveFirst` calls `thermfield.randomize()` every step
even at `temp 0`, so the production baseline generates `3·N·M` random numbers
per step that the deterministic adaptive Heun path never generates. The accepted
comparison uses the complete production step, because that is what a user pays.
The harness additionally reports a Hamiltonian-only floor (two `heisge` calls,
the work adaptive actually replaces) so a speedup arising only from skipping the
thermal field could not be mistaken for a coarse-graining result. At the
observed 3.4 × 10⁴ deficit the distinction is academic, but it is recorded because it
would matter if the defect in §3.3 were fixed.

---

## 5. Correctness on the measured commit

`ctest -L cg13-cuda` on the fresh out-of-tree `build_rcg09_cuda_fp64` passes
**29/29** — identical to RCG-08's fp64 baseline. The parity check inside the
benchmark (§2.1) additionally passes at `2.8e-16` at every system size measured.
Full commands and log in §7.

The only production-source changes this task made are the phase-timing opt-out
and its production selector; both are measurement controls that leave kernel
order, launches, and results identical. The benchmark and harness changes are
test-only.

---

## 6. Checklist reconciliation

Reconciled honestly; items that are not met are left unchecked rather than
reworded to fit.

| Blueprint item | State |
| --- | --- |
| The baseline calls the real production atomistic path | **Met.** `GpuHamiltonianCalculations` + `GpuDepondtIntegrator` in process; the ordinary `sd` binary end to end. |
| Baseline and adaptive cases use identical supported physics | **Met.** Verified numerically at `2.8e-16`, not asserted; e2e inputs differ only in the adaptive block. |
| Synthetic inactive-overhead control renamed and separated | **Met.** Separate source, separate target, self-declaring output. |
| Setup, warm-up, and steady-state scopes explicit | **Met.** Two-point subtraction e2e; explicit warm-up and repetition loops in process; setup reported separately. |
| CPU and GPU use appropriate wall/device timing | **Met.** Host steady-clock wall for totals, device events for phases. |
| Active atom, block, interface, and bond fractions reported | **Met.** Including live-bond fraction, which is the quantity §8.9 of the parent actually asks about. |
| Every relevant phase and host wait reported | **Met.** Eight phases, launch counts, phase synchronizations, compaction host wait. |
| Median, MAD, repetitions, and raw samples exist | **Met.** All emitted per point. |
| Device, driver, compiler, flags, clocks, precision recorded | **Met.** §4.2 and the harness metadata lines. |
| CUDA and HIP results separate | **Met vacuously.** CUDA only; HIP unavailable and marked so. |
| CPU thread count and affinity recorded | **Met.** |
| Crossover computed against the production atomistic baseline | **Met.** And the answer is that there is none. |
| Speedup exceeds its uncertainty margin | **NOT MET — there is no speedup.** Left unchecked. |
| No result extrapolated to an unavailable backend | **Met.** |
| All physics/parity fixtures pass on the measured commit | **Met.** 29/29 `cg13-cuda` on a fresh out-of-tree build; §7. |
| Human approves any restored performance wording | **Open.** No performance wording is restored by this task; the recommendation is that the `1.30x` claim be withdrawn. |

---

## 7. Fresh-build and suite evidence

Commit measured: `c00c3097` plus this task's working tree.

```text
cmake -S . -B build_rcg09_cuda_fp64 -DUPPASD_GPU_BACKEND=CUDA -DCMAKE_BUILD_TYPE=Release
cmake --build build_rcg09_cuda_fp64 -j 6
ctest --test-dir build_rcg09_cuda_fp64 -L cg13-cuda --output-on-failure
```

Fresh out-of-tree directory, no reused modules or objects. Build succeeded;
`100% tests passed, 0 tests failed out of 29`, total 310.11 s. This matches
RCG-08's recorded 29/29 fp64 baseline exactly — the phase-timing opt-out
changed no result, which is the point of leaving it default-enabled.

Full log: `docs/rcg09/ctest_cg13_cuda.log`.

Working-tree changes on the measured commit:

| File | Kind |
| --- | --- |
| `source/gpu_files/gpuAdaptiveRuntime.{hpp,cpp}` | production: phase-timing opt-out (`beginPhase`/`finishPhase`/`updateBlockState`) |
| `source/gpu_files/gpuSimulation.cpp` | production: `UPPASD_ADAPTIVE_PHASE_TIMING` selector and its status line |
| `tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp` | test-only: rewritten as PERF-ATOMISTIC-PROD/PERF-CG-SWEEP |
| `tests/coarse_graining/benchmark_gpu_inactive_runtime_overhead.cpp` | test-only: separated synthetic control |
| `tests/coarse_graining/run_rcg09_perf_e2e.py`, `tests/coarse_graining/perf/` | test-only: e2e harness and fixtures |
| `CMakeLists.txt` | test-only: new microbenchmark target |

### Raw artifacts

| File | Contents |
| --- | --- |
| `docs/rcg09/ctest_cg13_cuda.log` | complete 29/29 suite log |
| `docs/rcg09/scaling_with_energy_atomic.txt` | §3.2 sweep, full harness output at five sizes |
| `docs/rcg09/scaling_without_energy_atomic.txt` | §3.3 negative control, three sizes |
| `docs/rcg09/e2e_two_point.txt` | §3.1 both arms, ordinary `sd` binary, two-point |
| `docs/rcg09/device_state.txt` | `nvidia-smi` device and foreign-process state |

---

## 8. Follow-ups opened by this task

* **RCG-09-FU1 — single-address FP64 energy accumulation, everywhere it
  appears. Promoted to task RCG-09B on 2026-08-14.** §3.3 and §3.5. `atomicAddEnergyTerm` is a compare-and-swap loop on
  one global `double`, and eight call sites contend on it:
  `evaluateAdaptiveAtomisticBonds` (`[0]`, one per bond),
  `evaluateAdaptiveAtomisticOnsite` (`[1]`, one per active atom),
  `evaluateAdaptiveCoarseTensor` (`[2]`, `[3]`) and
  `finalizeAdaptiveCoarseLocal` (`[4]`, `[5]`, one per active block), and the two
  dipole kernels (`[6]`). Fixing only the bond site would move the bottleneck to
  the coarse site and leave the feature no better off — §3.5 shows both ends of
  the sweep are already blocked. The fix is a per-thread-block reduction into
  scratch followed by a tree reduction, which is what the *field* accumulators
  already do per atom; RCG-06B's unconditional FP64 energy contract must be
  preserved exactly, so removal is not an option. This blocks any future speedup
  claim. It is deliberately **not** fixed here: RCG-09 is a measurement task and
  the blueprint requires one defect class per patch.
* **RCG-09-FU2 — quadratic unique-bond construction at setup.**
  `build_unique_bonds` in `source/CoarseGraining/adaptivecgproduction.f90`
  scans every bond found so far for every directed neighbour, which is
  `O(directed × bonds)`. Measured setup cost at 16 384 atoms is ~38.9 s, against
  a complete 1100-step feature-off run of 0.87 s. It is a one-off cost, so it does not affect
  the per-step result above, but it dominates short production runs. This is the
  same defect class as F-08, which RCG-07 fixed for dilation, in a routine
  RCG-07 did not touch.
* **RCG-09-FU3 — clean-device re-measurement.** §4.1. Inherits RCG-08-FU3.
  Required before any absolute figure here is quoted as final. Does not affect
  the qualitative conclusion.
* **RCG-09-FU4 — HIP measurement.** Inherits RCG-08-FU1 / RCG-04-FU1 /
  RCG-05-FU1 / RCG-06-FU1.
* **RCG-09-FU5 — e2e fine-fraction sweep.** The e2e harness currently compares
  the selector-driven adaptive arm against feature-off. A static-mask sweep at
  controlled coarse fractions would let the e2e arm reproduce the in-process
  sweep end to end. Not required for the present conclusion, which is dominated
  by a fraction-independent cost.

---

## 9. Recommendation

Under remediation blueprint §13, items 1–7 and 10–12 may still be satisfiable,
but **item 9 fails**: no speedup claim survives comparison with the real UppASD
atomistic production path. The prescribed consequence is explicit — the result
may be accepted only as a correctness/reference implementation, the parent
blueprint must leave the GPU production-performance gates open, and speedup
wording must be removed.

Concretely, for RCG-10:

1. Withdraw the `1.30x` active-DOF crossover from the parent CG-10 evidence and
   replace it with a reference to this document.
2. Leave parent CG-10's "A real active-DOF crossover is measured" box unchecked.
3. Describe the GPU backend as a correctness prototype, as the remediation
   blueprint's §2.4 claim discipline has required throughout.
4. Treat RCG-09-FU1 as the prerequisite for any future performance work. Until
   it lands, no crossover is reachable at any system size, and re-running this
   benchmark will keep saying so.
