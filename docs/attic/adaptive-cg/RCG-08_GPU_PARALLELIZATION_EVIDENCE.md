# RCG-08 evidence: parallelize the adaptive GPU production path

**Task:** RCG-08 (`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` §9)
**Findings addressed:** F-09 (six adaptive GPU kernels explicitly single-threaded),
F-12 (Hillis–Steele compaction scan performs `O(N log N)` element work over
`O(log N)` launches)
**Date:** 2026-08-13
**Status:** Evidence recorded; **not closed** — see §9.

---

## 1. Dependency status, stated up front

RCG-08's declared dependency is RCG-07. At the time this session began:

- **RCG-07 was uncommitted.** Its CPU changes
  (`source/CoarseGraining/*.f90`, `CMakeLists.txt`,
  `tests/coarse_graining/test_static_hybrid_dilation_scaling.f90`, and
  `docs/RCG-07_CPU_ALGORITHMIC_SCALING_EVIDENCE.md`) were present only as
  working-tree modifications on top of `a1cd5bd3`.
- **RCG-07 is explicitly "not yet closed"** in the remediation blueprint: no
  independent reviewer distinct from its implementer has examined its
  race-freedom reasoning.

This session did not commit, alter, or re-review RCG-07. RCG-08's changes are
confined to `source/gpu_files/` and `tests/coarse_graining/`, which are
disjoint from RCG-07's file set, so the two remain separable into distinct
commits as §5.2 of the blueprint requires. **RCG-08 cannot be accepted before
RCG-07 closes**, independently of anything below.

The blueprint's sequencing rule was raised with the Human owner before work
started; the instruction was to do RCG-08 rather than RCG-09, precisely
because RCG-09's measurement is meaningless while F-09 stands.

---

## 2. Finding reproduction on the pre-fix commit

Working tree at `a1cd5bd3` plus RCG-07's uncommitted CPU changes. Fresh
out-of-tree build, no reused module or object files:

```bash
cmake -S . -B build_rcg08_base -DUPPASD_GPU_BACKEND=CUDA \
      -DUPPASD_PRECISION=DOUBLE -DCMAKE_BUILD_TYPE=Release
cmake --build build_rcg08_base -j 4
```

### 2.1 F-09 reproduced by direct source inspection

All six kernels named in F-09 began with the same guard:

| Kernel | Pre-fix line | Guard |
| --- | ---: | --- |
| `restrictAdaptiveMoments` | 152 | `if(blockIdx.x != 0 \|\| threadIdx.x != 0) return;` |
| `publishAdaptiveState` | 433 | same |
| `evaluateAdaptiveAtomistic` | 671 | same |
| `finalizeAdaptiveEnergy` | 1099 | same |
| `predictorAdaptiveHeun` | 1121 | same |
| `correctorAdaptiveHeun` | 1165 | same |

Each was launched `<<<1, 1, 0, stream_>>>`, so the guard was not merely
redundant: the launch geometry itself was one thread.

### 2.2 F-09 quantified

`gpu_adaptive_runtime_benchmark --blocks 2048 --atoms-per-block 4
--repetitions 7`, fp64, RTX A4000. Raw artifact:
`docs/rcg08/bench_pre_fp64.txt`.

At an all-fine mask (`requested_fraction=1.000000`, 8192 atoms, 24576 bonds),
the single-threaded atomistic kernel accounted for **103 530 µs of a
109 382 µs field evaluation — 94.6%**. Over a complete Heun step, the
single-threaded atomistic and integration phases together accounted for
159 411 µs of 178 268 µs (**89.4%**).

The reported `1.30x` "crossover" in that same artifact compares the adaptive
runtime's all-fine mode against its own coarsened mode — not against UppASD's
production atomistic path. That is F-10's finding and RCG-09's problem, not
this task's; it is named here only so the pre-fix numbers are not mistaken for
a production baseline.

### 2.3 F-12 reproduced

`launchCompaction()` ran `ceil(log2(N))` `scanAdaptiveWorkStep` launches, each
reading and writing all `3N` flags. Measured directly by the new negative
control in §5.3: **19 launches** at `N = 76800`.

---

## 3. What changed

All changes are in `source/gpu_files/gpuAdaptiveRuntime.{cpp,hpp}` plus the two
test/benchmark files. No Fortran, no CPU operator, and no physics convention
was touched.

### 3.1 The six F-09 kernels

| Pre-fix kernel | Post-fix | Work decomposition |
| --- | --- | --- |
| `restrictAdaptiveMoments` | same name, rewritten | one thread per `(channel, block, ensemble)` scalar, gathering that block's CSR atom range |
| `publishAdaptiveState` | same name, rewritten | one thread per spatial block |
| `evaluateAdaptiveAtomistic` | `clearAdaptiveAtomistic`, `evaluateAdaptiveAtomisticBonds`, `evaluateAdaptiveAtomisticOnsite`, `writebackAdaptiveAtomistic` | one thread per field component; per `(bond, ensemble)`; per `(active-atom slot, ensemble)`; per `(active-atom slot, ensemble)` |
| `predictorAdaptiveHeun` | `predictorAdaptiveHeunAtoms`, `predictorAdaptiveHeunBlocks` + two stream-ordered device copies | one thread per `(active-atom slot, ensemble)` and per `(active-block slot, ensemble)` |
| `correctorAdaptiveHeun` | `correctorAdaptiveHeunAtoms`, `correctorAdaptiveHeunBlocks` | as above |
| `finalizeAdaptiveEnergy` | unchanged behaviour, guard rewritten | **deliberately left single-threaded** — see §7.1 |

`grep -n 'blockIdx.x != 0\|threadIdx.x != 0' source/gpu_files/gpuAdaptiveRuntime.cpp`
now returns nothing.

Two design choices are worth calling out because they are what make the
parallel restriction and integration *bitwise* equal to the serial reference
rather than merely within tolerance:

- **Restriction gathers instead of scattering.** Rather than have every atom
  atomically add into its block's accumulator, each thread owns one output
  scalar and walks that block's CSR slice. `BlockTopology` fills `block_atoms`
  by scanning atoms `1..Natom` with a per-block cursor, so the slice is
  strictly ascending in atom index — the identical summands in the identical
  order the serial global atom loop used. This also removed the separate clear
  pass, since each thread writes its output completely.
- **The Heun stages keep the accepted update and the compact-list work
  definition verbatim.** Only the loop-to-thread mapping changed.

### 3.2 The F-12 compaction scan

`launchCompaction()` now runs a three-phase hierarchical scan
(`scanAdaptiveTiles` → recursive tile-total scan → `addAdaptiveTileOffsets`)
instead of a global Hillis–Steele sweep. Tile-total storage
(`compactionScanLevels_`) is sized at `allocate()` time from
`compactionScanLevelItems()`, the same helper `estimateBytes()` uses, so the
scan's temporary storage is part of the ordinary GPU memory preflight and the
two cannot drift. Nothing is allocated on the hot path.

The tile-local scan is Kogge–Stone in shared memory: `T·log2(T)` operations
per `T = 256` elements. `T` is a compile-time constant, so this is a fixed
factor of 8 per element, not a growing one — total element work is `O(N)` and
the launch count is `O(log_T N)`. This constant is stated rather than hidden;
a work-efficient Blelloch tile scan would reduce it to ~2, but the scan is far
from the hot path and the simpler formulation is easier to keep identical
across backends.

### 3.3 Portability

No CUB or rocPRIM. The scan is plain shared memory so CUDA and HIP compile the
same algorithm rather than two vendor primitives that would need separate
validation — and there is no HIP device here to validate on. Every launch site
this task introduced or rewrote goes through a single `ADAPTIVE_LAUNCH` macro,
so the CUDA and HIP spellings cannot diverge in grid shape or argument order.
Pre-existing launch sites were left in their explicit `#if/#elif` form to keep
this patch's diff to one defect class.

---

## 4. Race-freedom and determinism argument

This is the argument an independent reviewer should attack first; §5.4 records
the sanitizer evidence that tests it empirically.

### 4.1 Kernels that need no atomic at all

`restrictAdaptiveMoments`, `publishAdaptiveState`, and all four Heun kernels
are write-disjoint by construction:

- **Restriction** — thread `s` writes only `channelMomentSum[s]`,
  `coarseMoment[3s..3s+2]`, `coarseDirection[3s..3s+2]`. Distinct threads have
  distinct `s`.
- **Publication** — thread `b` writes `pendingState[b]`, `blockState[b]`,
  `stateAge[b]`, `transitionEpoch[b]`, and `atomDirection` /
  `transitionBackup` / `ghostDirection` only at atoms in block `b`'s CSR
  range. `validate()` already proves CSR membership is a bijection (each atom
  in exactly one block), so no two threads can address the same atom.
- **Heun stages** — the compact lists hold each atom/block at most once, so
  `(atom, ensemble)` and `(block, ensemble)` are unique per thread.

### 4.2 Kernels that do need atomics

Only the atomistic Hamiltonian:

- `evaluateAdaptiveAtomisticBonds` scatters onto both endpoints of each bond,
  and two bonds can share an endpoint. Endpoint accumulation uses `atomicAdd`
  on the field scratch. This is the same trade RCG-07 accepted for the CPU
  analogue of this loop: per-component atomics rather than an atom-sized
  per-thread reduction array, which would multiply an `O(natoms)` array by the
  thread count and undo RCG-06A. Unique-pair ownership is unchanged — each
  bond is still visited exactly once.
- Energy terms `[0]` and `[1]` accumulate through `atomicAddEnergyTerm`, the
  RCG-06B FP64 CAS accumulator, so they stay FP64 regardless of build
  precision. Their summation order becomes atomic-arrival order — exactly what
  the coarse energy terms `[2..6]` have done since CG-10.

`evaluateAdaptiveAtomisticOnsite` needs no atomic for the *field* (each active
atom is owned by one thread) and accumulates its axes into a thread-local
`double` before a single atomic energy add.

### 4.3 Why transitions stay deterministic

The blueprint requires deterministic accepted state transitions. Atomic
arrival order is not deterministic, so it matters exactly which quantities a
transition decision reads:

1. **Selector scores** — `atomicMaxSelector`. Maximum is associative,
   commutative and exact in floating point, so the result is independent of
   arrival order.
2. **The polarization gate** — a pure function of `restrictAdaptiveMoments`'
   outputs, which §3.1 establishes are bitwise identical to the serial
   reference.
3. **The reconstruction RNG** — still seeded from the
   `(globalSeed, block, channel, ensemble, epoch)` tuple, so each block draws
   its own stream regardless of the order blocks now execute in.
4. **The energy-jump gate** — the production GPU path
   (`GpuSimulation::advanceAdaptiveStep`) passes `acceptedBlockMask = nullptr`,
   so no production transition decision reads an atomically-accumulated
   energy at all.

No transition decision depends on an atomically-accumulated quantity.
Energies and fields therefore carry atomic-order noise at the last ulp; state
decisions and trajectories do not.

### 4.4 Ordering between the split kernels

The four atomistic launches must run in order: clear zeroes the scratch that
the bond scatter accumulates into atomically; the on-site pass then subtracts
into the same scratch non-atomically at atoms it owns exclusively; the
writeback is last. Consecutive launches on one stream do not overlap, which
supplies exactly that ordering — the same property RCG-05E already relies on
for the dilation double buffer. The predictor's two save copies are issued as
stream-ordered device-to-device transfers before the predictor kernels, so
they are equally ordered.

---

## 5. Tests and negative controls

### 5.1 Environment

| | |
| --- | --- |
| Commit | `a1cd5bd3` + RCG-07's uncommitted CPU changes + this task's changes |
| Device | 2× NVIDIA RTX A4000 (sm_86), driver 610.57.04 |
| Toolkit | CUDA 13.3 (V13.3.73), `-arch=native` → sm_86 |
| Host compiler | GNU Fortran / GCC 13.3.0, CMake 3.28.3 |
| Flags | `-O3 -DNDEBUG -std=c++20`, `CMAKE_BUILD_TYPE=Release` |
| Precisions | fp64 (`UPPASD_PRECISION=DOUBLE`) and fp32 (`SINGLE`) |
| HIP | **unavailable** — no toolchain (`hipcc`, `rocminfo` absent) and no device |

### 5.2 Regression suites

Fresh out-of-tree builds, no reused module or object files.

| Suite | Pre-fix | Post-fix |
| --- | --- | --- |
| `ctest -L cg13-cuda`, fp64 | 29/29 pass | **29/29 pass** |
| `ctest -L cg13-cuda`, fp32 | not measured | **28/29** — see §5.5 |

The fp64 suite includes the `backend-parity` (3 tests), `moving-parity` (3),
`ownership` (3), `production` (1) and `rejection` (1) labels, so CPU/GPU field,
energy, trajectory and ownership parity, and the RCG-04 moving fixtures, are
all covered by that count.

A pre-fix fp32 suite run was **not** performed, so no pre/post fp32 comparison
is claimed. The fp32 evidence here stands on its own: the post-fix fp32 suite
result, plus the individual diagnosis of its one remaining failure (§5.5). The
two fp32 failures observed *before* the RCG-05-FU4 fix were on post-fix RCG-08
code, not on a pre-fix baseline, and are described as such in §5.5.

### 5.3 New test: `testHierarchicalCompaction` (F-12 negative control)

Added to `tests/coarse_graining/test_gpu_adaptive_runtime.cpp`. Every
pre-existing fixture in that file has 8 atoms and 4 blocks — a single
256-element tile — so none of them reaches the tile-total level or the offset
propagation between levels, and all of them would keep passing if either were
wrong.

The new fixture uses 300 blocks × 256 basis sites = **76 800 scan items**,
which forces two tile-total levels (76 800 → 300 → 2) and therefore the
complete down-propagation chain. Block states follow a coprime-strided
irregular pattern so runs of fine/coarse blocks straddle tile boundaries.

It asserts two independent things:

1. **Correctness** — all three compact work lists match a host reference
   computed directly from the block states.
2. **Complexity** — the compaction launch count stays near `log_256(N)`.

**Negative control, executed rather than argued.** `launchCompaction()` was
temporarily reverted to the pre-fix Hillis–Steele sweep (keeping the new
launch counter) and the suite rebuilt and rerun:

```text
terminate called after throwing an instance of 'std::runtime_error'
  what():  compaction still scales its launch count like log2(N) (F-12):
           got 19 launches for 76800 scan items
```

19 launches, matching the hand-derived `ceil(log2(76800)) + 2 = 17 + 2`. The
fixed build reports 7 and passes. The source was then restored and re-verified
(`grep -c scanAdaptiveWorkStepLegacy` → 0).

Two things are worth recording about that run. First, the *correctness*
assertions passed under **both** algorithms — the legacy sweep reached the
launch-count check, which means the new hierarchical scan and the old sweep
produce byte-identical compact lists on this 76 800-item irregular fixture.
That is a direct equivalence check between the two implementations, not just a
check of the new one against a host reference. Second, it confirms the
correctness assertions alone could not have distinguished the two: the old
scan was correct, only slow, so the complexity assertion is doing the real
discriminating work here.

### 5.4 Sanitizers (fp64, CUDA)

`compute-sanitizer` from CUDA 13.3 against the complete
`gpu_adaptive_runtime_tests` binary, which exercises restriction, both field
evaluations, selector, dilation, publication with accepted and rejected
transitions, both reconstruction schemes, a Heun step, the polarization gate,
and the new 76 800-item compaction fixture:

| Tool | Result | Raw artifact |
| --- | --- | --- |
| `memcheck --leak-check full` | **0 errors, 0 bytes leaked in 0 allocations** | `docs/rcg08/san_memcheck.txt` |
| `racecheck --racecheck-detect-level info` | **0 hazards displayed (0 errors, 0 warnings)** | `docs/rcg08/san_racecheck.txt` |
| `initcheck` | **0 errors** | `docs/rcg08/san_initcheck.txt` |
| `synccheck` | **0 errors** | `docs/rcg08/san_synccheck.txt` |

`initcheck` is the one that matters most for this patch: parallel restriction
*removed* the separate clear pass on the grounds that each thread now writes
its output scalar completely, and `initcheck` is what would catch that
reasoning being wrong. `racecheck` covers the shared-memory Kogge–Stone tile
scan, the only new shared-memory code.

### 5.5 fp32, and one pre-existing bug fixed to unmask it

The first fp32 run (post-fix RCG-08 code, before the fix described below)
showed two failures. Both were checked individually rather than attributed to
their tracked predecessors by assumption:

1. **`coarse-graining-gpu-adaptive-runtime` (subprocess aborted)** — run
   directly, the message was
   `descriptor layout: refineThreshold shifted across a copy`. That is
   **RCG-05-FU4** exactly: `testSelectorPolicyDescriptorLayout()` compared a
   `real`-valued field against a hardcoded `double` literal, and
   `float(0.123456)` does not round-trip to `double 0.123456`. It is a
   test-only defect, already confirmed by RCG-05G not to be a real
   descriptor-layout problem.

   Because that test runs first in `main()`, its abort **masked every later
   test in the binary — including all of RCG-08's own fp32 evidence**. This is
   the same masking RCG-06E disclosed for its GPU memory-preflight item and
   left in place.

   Rather than inherit the mask, RCG-05-FU4 was fixed here, exactly as its own
   task description specifies: compare `copied.*Threshold` against the pre-copy
   `policy.*Threshold` value rather than a literal. A sibling assertion that
   the two sentinels remain distinguishable was added so the new comparison
   cannot pass on an aliased field, and the reverse check on the
   default-constructed policy — which passed at fp32 only by accident — was
   corrected the same way. **This is a three-line test-only change in a
   file this task was already editing, carrying no production behaviour**, and
   it is called out separately here so it can be reviewed (or reverted)
   independently of RCG-08's own work.

2. **`adaptive-cg-production-e2e`** — this is the pre-existing **RCG-04-FU5**
   CUDA fp32 `gpu_fft_static_mixed` dipole-energy assertion, already
   reconfirmed at fp32 by both RCG-05G and RCG-06E. Not introduced here and
   not investigated here.

With RCG-05-FU4 fixed, `gpu_adaptive_runtime_tests` passes at fp32 —
including `testKernelParityAndWorkflow`, `testPolarizationGate` and the new
`testHierarchicalCompaction` — so RCG-08's parallel kernels have fp32 evidence
in their own right and not merely by fp64 extrapolation.

Resulting fp32 suite state: **`ctest -L cg13-cuda` 28/29**, the single failure
being RCG-04-FU5 above. Its log line is verbatim the assertion RCG-04-FU5
names:

```text
assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0
AssertionError
```

The same fp32 run also shows the RCG-04I frozen fp32 budgets still
discriminating correctly — all three moving-parity negative controls
(`freeze`, `perturb-one-component`, `+q vs -q chirality`) fail as designed
with `max_radians` far above the 0.0015 budget — so parallelization did not
quietly widen the precision envelope those fixtures police.

---

## 6. Measurements

### 6.1 Measurement environment — read this before the numbers

**This is a shared machine and it was not idle.** Throughout the measurement
window another user's processes held ~4 GB on each of the two RTX A4000s and
drove both to 84–91% utilization, with SM clocks throttled to 1230–1560 MHz
against a 2100 MHz maximum at 92–93 °C. Per-GPU state was sampled at each
round into `docs/rcg08/ab_env.txt`.

Consequences, stated plainly:

- **Absolute timings below are contended and are not a clean device
  measurement.** Device-event phase timers include time the context spends
  descheduled while time-slicing against the other workload.
- To keep the comparison meaningful, pre-fix and post-fix binaries were run
  **interleaved in alternating order across three rounds** (`docs/rcg08/ab.sh`
  ordering, raw output in `ab_pre.txt` / `ab_post.txt`), so both variants
  experienced comparable contention. The **ratios** are therefore the
  defensible result; the absolute microsecond figures are not.
- Every conclusion drawn here is a **ratio or a within-run fraction** — the
  pre/post speedups in §6.2 and the "94.6% of the field step" figure in §2.2,
  which normalizes a phase against its own run. No conclusion rests on an
  absolute microsecond value, and none should be quoted as one.
- §2.2's pre-fix baseline was taken earlier in the session and its contention
  state was not sampled, which is a second reason to treat §6.2's interleaved
  pairing, not §2.2, as the pre/post comparison of record.

Configuration: `--blocks 2048 --atoms-per-block 4 --repetitions 5
--iterations 10`, fp64, 8192 atoms, 24576 bonds.

### 6.2 Paired interleaved A/B, fp64

All-fine mask (`requested_fraction = 1.0`) — the case F-09 dominated:

| Round | Pre-fix wall (µs) | Post-fix wall (µs) | Pre atomistic (µs) | Post atomistic (µs) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 110 951 | 40 560 | 105 565 | 34 627 |
| 2 | 109 522 | 37 381 | 104 270 | 32 144 |
| 3 | 112 485 | 40 183 | 107 888 | 34 157 |
| **median** | **110 951** | **40 183** | **105 565** | **34 157** |

- Field evaluation: **2.76× faster**
- Atomistic phase alone: **3.09× faster**

Complete two-stage Heun step (`adaptive-step` line):

| Round | Pre-fix wall (µs) | Post-fix wall (µs) | Pre integration (µs) | Post integration (µs) |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 175 024 | 35 980 | 44 321 | 58.96 |
| 2 | 160 802 | 35 324 | 42 074 | 60.76 |
| 3 | 154 740 | 34 651 | 38 228 | 52.75 |
| **median** | **160 802** | **35 324** | **42 074** | **58.96** |

- Complete step: **4.55× faster**
- Integration phase (restriction + both Heun stages): **714× faster**

The integration figure is the clearest signature of F-09: those loops walked
8192 atoms on one lane and did essentially no arithmetic per atom, so
serialization was the entire cost. The atomistic phase improves by "only" 3.1×
because it is now bound by something else — see §6.5.

Ranges do not overlap between pre- and post-fix in any round for any of the
four quantities, so the improvement is far outside the round-to-round spread
even in this contended environment.

### 6.3 Launches and synchronizations per complete Heun step

New `adaptive-launches` benchmark line, identical across all three post-fix
rounds:

```text
adaptive-launches per_step(atomistic=8.00 coarse=10.00 interface=6.00
  selector=0.00 polarization=0.00 compaction=0.00 integration=4.00)
  total=28.00 phase_syncs=10.00
```

28 kernel launches and **10 blocking host synchronizations** per step.
`compaction=0` because the benchmark's `integrateHeun` loop performs no
selector update; the compaction path's launch count is measured separately by
`testHierarchicalCompaction` (§5.3): 7 launches at 76 800 items against the
pre-fix 19.

### 6.4 Theoretical occupancy (sm_86)

From `ptxas -v` on the production translation unit, raw output in
`docs/rcg08/ptxas.txt`. Occupancy is register-limited warps per SM against the
sm_86 maximum of 48, at the 256-thread launch geometry these kernels use:

| Kernel | Registers | Warps/SM | Occupancy |
| --- | ---: | ---: | ---: |
| `clearAdaptiveAtomistic` | 10 | 48 | 100% |
| `writebackAdaptiveAtomistic` | 16 | 48 | 100% |
| `scanAdaptiveTiles` | 16 (+1024 B smem) | 48 | 100% |
| `addAdaptiveTileOffsets` | 16 | 48 | 100% |
| `finalizeAdaptiveEnergy` | 20 | 48 | 100% |
| `restrictAdaptiveMoments` | 40 | 48 | 100% |
| `evaluateAdaptiveAtomisticBonds` | 40 | 48 | 100% |
| `predictorAdaptiveHeunAtoms` / `Blocks` | 44 | 42 | 87.5% |
| `correctorAdaptiveHeunAtoms` / `Blocks` | 45 | 42 | 87.5% |
| `evaluateAdaptiveAtomisticOnsite` | 56 | 36 | 75% |
| `publishAdaptiveState` | 110 | 18 | 37.5% |

`publishAdaptiveState` is the outlier: the constrained-cone reconstruction is
register-hungry. It runs only at selector-update intervals rather than every
step, and the aligned default path is much cheaper, so it was not tuned. This
is theoretical occupancy from register pressure, not achieved occupancy —
achieved occupancy was not measured because the device was never idle enough
for a profiler run to mean anything (§6.1).

### 6.5 Where the atomistic phase time now goes

The atomistic phase improved 3.09×, not by the ~3 orders of magnitude the
thread-count change might suggest. Two effects, both visible in the artifacts:

1. **Unaccounted (host) time barely moved**: 14.69 ms/step pre-fix,
   13.54 ms/step post-fix. It was 9.1% of the pre-fix step and is **38.3% of
   the post-fix step**. With 10 blocking phase synchronizations per step (§6.3)
   plus a synchronous 64-byte device-to-host energy copy inside every
   `evaluateHybrid`, this is host round-trip cost, not device work — and under
   90% contention from another process each round trip is expensive.
2. **Device contention** inflates the device-event phase timers themselves.

Both point the same way: after RCG-08, the adaptive GPU step is no longer
kernel-serialization-bound, and the next real target is the per-phase
synchronization structure (§7.2), not more thread-level parallelism.

---

## 7. Remaining serial sections and known costs

RCG-08's prompt asks for remaining serial sections to be reported, not merely
for the six kernels to be fixed. These are the ones that remain.

### 7.1 `finalizeAdaptiveEnergy` — single thread, deliberately

Seven FP64 loads, six adds, one store, of accumulators every other
energy-producing kernel has already reduced in parallel. It is `<<<1,1>>>`
because there is exactly one element of work, not because `O(N)` work was
serialized onto one lane — which is what F-09 objected to and what the other
five kernels actually did. Parallelizing a seven-term sum would add a launch
and a reduction buffer to save nanoseconds.

This means the RCG-08 checklist item "no production hot kernel is guarded to a
single device thread" is met **in substance but not literally**: one kernel
still runs on one thread. It is recorded here rather than papered over, and
the guard was rewritten to use the shared thread-index helper so it stays
correct if the launch geometry is ever widened.

### 7.2 Per-phase host synchronization — the largest remaining serial section

`finishPhase()` blocks the host on an event at every phase boundary. That is
how RCG-06C obtains per-phase device times, but it also means adjacent phases
can never overlap, and each boundary costs a full host round trip. The new
`phaseSynchronizations` counter measures this directly (§6.3).

This was invisible while the six kernels dominated. Now that they do not, it
is a leading cost. **Removing it is deliberately out of RCG-08's scope**: the
per-phase device timers are RCG-06C's accepted evidence, and making them
optional changes what RCG-06C measured. It is recorded as a follow-up rather
than changed here.

### 7.3 Full-system `O(N)` passes that survive

Independent of active fraction, each step still performs:

- the predictor's two state saves — `3·atoms·ensembles` and
  `3·channels·blocks·ensembles` reals, now stream-ordered device copies rather
  than a serial kernel loop, but still full-system;
- `clearAdaptiveAtomistic` over `3·atoms·ensembles` field components;
- `prolongateAdaptiveGhosts` / `restrictAdaptiveInterface` /
  `commitAdaptiveGhosts` over `atoms·ensembles`;
- `initializeAdaptiveWorkScan` and `scatterAdaptiveWork` over
  `max(atoms, blocks)`.

These are bandwidth-bound rather than serialized, but they set a floor on how
far the coarse fraction can reduce the step cost, and they are the natural
next target after RCG-08. The CPU analogue of this observation is RCG-07's own
documented remaining `O(n_blocks)` pass.

### 7.4 Compact lists sized by upper bound, not by count

Kernels that walk a compact list are launched with a grid sized by the
*maximum* possible work (`atoms·ensembles` or `blocks·ensembles`) and bound-check
each thread against `runtime.workCounts[·]`, which lives in device memory.
This avoids a host round trip to read the counts, at the cost of launching
threads that return immediately when the active fraction is small. The real
Hamiltonian and integration *work* is still reduced — only the launch
geometry is not. Device-side launch sizing would need a device-side launch or
a persistent-kernel formulation; neither is attempted here.

### 7.5 Divergence in cone reconstruction

`publishAdaptiveState` runs a 100-iteration bisection per transitioning block
in the constrained-cone path. One thread per block means neighbouring lanes in
a warp may take very different iteration counts. This only affects blocks
actually transitioning under a nonzero cone angle, which no accepted fixture
uses today (`cg_cone_angle` defaults to 0, and RCG-06-FU2 tracks nonzero-cone
acceptance). Recorded so it is not rediscovered as a surprise if that
follow-up ever lands.

### 7.6 FP64 `atomicAdd` compute-capability requirement

`evaluateAdaptiveAtomisticBonds` uses `atomicAdd` on `real*`, which in an fp64
build is `atomicAdd(double*, double)` and requires compute capability ≥ 6.0.
This constraint is **pre-existing, not introduced here** —
`restrictAdaptiveInterface` and `atomicAddDerivativeStencil` have used the same
construct since CG-10. It is noted because `gpuAtomicDouble.hpp` documents
that this codebase enforces no minimum compute capability, so the constraint
is real but unenforced; the energy accumulator deliberately uses a CAS loop
for exactly this reason while the field accumulators do not.

---

## 8. Checklist reconciliation

Each RCG-08 checklist item against direct evidence in this document.

| # | Item | State | Evidence |
| --- | --- | --- | --- |
| 1 | No production hot kernel is guarded to a single device thread | **not met as literally stated** — left unchecked in the blueprint | `grep` for the pre-fix guard returns nothing and all six F-09 kernels are parallel (§3.1), so the finding's substance is closed; but `finalizeAdaptiveEnergy` is a production kernel on the hot path that still runs on one thread, because it is an `O(1)` seven-term sum rather than serialized `O(N)` work (§7.1). Ticking this would require restating the item to fit |
| 2 | Restriction covers all active atoms/blocks in parallel | met | §3.1, one thread per `(channel, block, ensemble)`; parity fixture §5.2 |
| 3 | Atomistic unique-pair fields and energies are race-safe | met | §4.2 argument; racecheck/memcheck/initcheck/synccheck all 0 (§5.4) |
| 4 | Predictor and corrector integrate compact active work in parallel | met | §3.1; norm-preservation and trajectory fixtures in §5.2 |
| 5 | Publication preserves deterministic accepted state transitions | met | §4.1 write-disjointness, §4.3 determinism argument; the existing two-runtime bitwise cone-determinism fixture passes |
| 6 | Energy finalization uses the accepted FP64 reduction | met | unchanged RCG-06B FP64 accumulators; §7.1 |
| 7 | Compaction performs linear work and uses tracked temporary storage | met | §3.2; preflight via shared helper; negative control 19 → 7 launches (§5.3) |
| 8 | Launch and synchronization counts are reported per phase | met | §6.3; new `adaptive-launches` benchmark line |
| 9 | FP64 fields, energies, trajectories, and decisions match CPU references | met | 29/29 `cg13-cuda` fp64 including `backend-parity` and `moving-parity` (§5.2) |
| 10 | FP32 has a separate accepted error budget | met (budget inherited, not derived) | RCG-04I's frozen fp32 budgets are separate from fp64's, unchanged by this task, and re-exercised at fp32 — including all three moving-parity negative controls still failing as designed (§5.5). No new budget was derived because no precision contract changed |
| 11 | Multiple ensembles and anisotropic buffer descriptors pass | met | `Fixture` uses 3 ensembles; `GEO-ANISO-BUFFER`/`ownership` label passes in §5.2 |
| 12 | CUDA sanitizer reports no memory or race errors where available | met | four tools, all 0 (§5.4) |
| 13 | HIP execution and sanitizer evidence are recorded where available | **not available** | no HIP toolchain or device; §5.1. Not extrapolated from CUDA |
| 14 | The serial reference, if retained, is test-only and clearly labelled | met, vacuously | **no serial reference was retained** — the serial kernels were replaced, not kept behind a switch, so nothing can accidentally dispatch them |
| 15 | Feature-off device inventory and timing remain unchanged within noise | **met for inventory; timing weak** | `inventory_delta_bytes=0` and `result=PASS` in all three post-fix rounds (`docs/rcg08/ab_post.txt`). See note below |

**Note on item 15.** The inventory half is exact and deterministic: an
inactive `GpuAdaptiveRuntime` allocates nothing, and that is asserted directly
rather than timed. The timing half passed its own adaptive 3·MAD budget in all
three rounds, but the observed deltas were `+20.5%`, `+0.3%` and `+43.7%` —
the budget widened with the noise rather than the measurement being tight.
Under §6.1's contention this control has weak discriminating power and should
not be cited as a strong feature-off timing guarantee. The structural argument
is stronger than the measurement here: this task changed no code on the
feature-disabled path, which executes none of the kernels it touched.

---

## 9. Why this is not closure

RCG-08 is **not closed** by this document. Three separate reasons, none of
which the evidence above can remedy on its own:

1. **Its dependency is not closed.** RCG-07 is uncommitted and explicitly
   awaiting an independent review of its race-freedom reasoning (§1). The
   blueprint's sequencing rule is that no task is accepted before all of its
   declared dependencies. This alone is decisive.
2. **Independent review is mandatory and has not happened.** RCG-08's task
   definition requires review by "Sol or Opus/Terra not responsible for the
   kernel patch" plus "Luna/Sonnet for device parity harness". The
   race-freedom and determinism argument in §4 is written by the same session
   that wrote the kernels; that is exactly the separation the blueprint
   forbids. §4 is structured as a set of falsifiable claims specifically so a
   reviewer can attack them.
3. **HIP evidence does not exist.** Item 13 is unavailable, not passing. The
   `ADAPTIVE_LAUNCH` macro makes the CUDA and HIP launch geometry identical
   *by construction* (§3.3), which is the strongest claim available without
   hardware — but it is a structural argument, not execution evidence, and it
   is not offered as a substitute.

Item 1 is additionally left unchecked in the blueprint rather than reworded to
fit: `finalizeAdaptiveEnergy` still runs on one thread (§7.1). The substance of
F-09 — `O(N)` work serialized onto a single lane in the hot adaptive path — is
closed, and the measurements in §6 are what that closure looks like; but the
checklist sentence as written is false, so it stays unchecked.

### 9.1 Follow-ups this task opened

- **RCG-08-FU1 — HIP execution and sanitizer evidence.** Re-run §5.2's suites,
  §5.4's four sanitizers and §6's benchmark on real HIP hardware. Mirrors the
  still-open RCG-04-FU1 / RCG-05-FU1 / RCG-06-FU1; may be picked up alongside
  them rather than duplicated. Dependencies: HIP toolchain/hardware only.
- **RCG-08-FU2 — per-phase synchronization structure.** §7.2: ten blocking
  host synchronizations per step, now ~38% of step wall time (§6.5). Removing
  or making them optional conflicts with RCG-06C's accepted per-phase timing
  evidence, so it needs its own scope decision. This is the single largest
  remaining cost in the adaptive GPU step and the natural predecessor to any
  RCG-09 measurement.
- **RCG-08-FU3 — clean-device re-measurement.** Every timing in §6 was taken
  on a machine another user was driving to ~90% utilization with thermal
  throttling (§6.1). The ratios are paired and interleaved and the effect size
  is large, but the absolute numbers are not a clean device measurement and
  RCG-09 must not reuse them as one.
- **RCG-08-FU4 — full-system `O(N)` passes.** §7.3 lists the per-step passes
  that remain proportional to total system size regardless of coarse fraction.
  They set the floor on achievable crossover and are the natural next
  optimization target.
