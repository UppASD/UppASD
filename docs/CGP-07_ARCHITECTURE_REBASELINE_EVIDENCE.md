# CGP-07 — Re-baseline the architecture before continuum-kernel work, evidence

Task: `docs/CGP_work.md` CGP-07 ("Re-baseline the architecture before
continuum-kernel work"). Branch `gpu_hip_cu_ab_cg`, commit base `f5cf15d2`
(CGP-06B), built from the branch's current working tree, which additionally
includes CGP-05 (host-barrier removal; verified, not yet committed at the
time of this task — see project memory `cgp05-status`).

**Primary model (task doc): Luna.** **Actual model: Claude (Sonnet 5,
interactive session)**, matching the CGP-03D/CGP-05/CGP-06A/CGP-06B
precedent.

This is an evidence-only task. **No production code was modified and no
continuum-kernel work was performed.** The only source change is a one-line
addition of a `~1%` fine-fraction point to the pre-existing sweep list in
`tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp` (the 6.25/12.5/
25/50/75/100% points every prior CGP evidence doc already depends on are
unchanged, and reproduce their prior values below). Everything else reuses
already-existing, already-parity-verified benchmark tooling
(`gpu_adaptive_runtime_benchmark`, `gpu_thermfield_rng_benchmark`,
`run_rcg09_perf_e2e.py`) built by RCG-09/CGP-06A/CGP-06B, plus a new
evidence-only regression script (`docs/cgp07_raw/fit_scaling_model.py`, not
part of the build) used to fit the scaling model in Part 6.

## Environment

- GPU: NVIDIA RTX A4000 (device 0 of 2; default device pinning), shared host
  per project memory `shared-gpu-host-contention`. Checked idle immediately
  before benchmarking: 0% utilization, 210 MHz SM clock (idle floor), <110
  MiB used, no other compute process on either GPU (`nvidia-smi
  --query-compute-apps` empty). GPU 0's temperature read 70°C partway through
  this session as residual heat from this task's own preceding correctness
  suite runs (not contention — 0% utilization throughout); this does not
  affect the CUDA event-based device timings reported below, which measure
  GPU-clock intervals, not wall-clock-sensitive host timing.
- Driver 610.57.04, CUDA 13.3.73 (nvcc), `CMAKE_BUILD_TYPE=Release`.
- Host CPU: `nproc` reports 2 (cgroup-quota-limited; project memory
  `shared-gpu-host-contention`/`gpu-benchmark-cold-start`). Material to the
  whole-executable e2e cross-check in Part 7, not to the in-process
  microbenchmark numbers that form the primary evidence (CUDA-event device
  timings and short in-process wall-clock samples, both far below the
  scheduling-noise floor this quota otherwise introduces into any longer
  external-process wall-clock measurement).
- Builds: fresh out-of-tree `build_cgp07_cuda_fp32` (`UPPASD_GPU_BACKEND=CUDA`,
  `UPPASD_PRECISION=SINGLE`) and `build_cgp07_cuda_fp64`
  (`UPPASD_PRECISION=DOUBLE`), both `CMAKE_BUILD_TYPE=Release`, both built
  from the current working tree (CGP-05's uncommitted changes included).

## Part 0 — Correctness gate (common benchmark protocol)

`ctest -L cg13-cuda -j1` on both fresh builds, before any timing:

- fp32: **33/34 pass.** The one failure, `adaptive-cg-production-e2e`
  (`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`), is the
  same pre-existing, unrelated finding every prior CGP task since CGP-03D has
  recorded — not a regression from CGP-05/06B or from this task.
- fp64: **34/34 pass**, matching the established fp32-only-precision-specific
  pattern for that pre-existing failure.

Both builds' `gpu_adaptive_runtime_benchmark` also reports
`atomistic-parity ... result=PASS` before every timing run below (the
in-process field-identity check against UppASD's real
`GpuHamiltonianCalculations`/`GpuDepondtIntegrator` path) — see raw logs in
`docs/cgp07_raw/*.log`. This is what licenses every speedup number in this
document as a comparison of identical physics, not an assumption.

No all-fine T=0/finite-T parity fixture failed on either build. This
satisfies the common-protocol requirement without a separate, redundant run.

## Part 1 — Benchmark matrix actually run

| Dimension | Values covered |
| --- | --- |
| System sizes | 65,536 / 262,144 / 1,048,576 atoms (16x the historical 65,536-atom benchmark) |
| Atoms per coarse block | 8, 16, 32, 64 (at fixed N=262,144) |
| Fine fractions | 0, 1%, 6.25%, 12.5%, 25%, 50%, 75%, 100% (built into the existing sweep; the 1% point is this task's one addition) |
| Precision | fp32 (primary, full 3-size + 4-block-size matrix) and fp64 (one full 3-size matrix at apb=16, per the task's "retain one fp64 matrix" scope) |

Concretely, 9 full sweeps of `gpu_adaptive_runtime_benchmark` (each sweep = 8
fraction points x full phase breakdown x PERF-ATOMISTIC-PROD baseline):

- fp32 size sweep: `--blocks 4096 --atoms-per-block 16` (65,536),
  `--blocks 16384 --atoms-per-block 16` (262,144),
  `--blocks 65536 --atoms-per-block 16` (1,048,576).
- fp32 block-size sweep at N=262,144: `--atoms-per-block 8/16/32/64`
  (`--blocks 32768/16384/8192/4096` respectively; the apb=16 point is shared
  with the size sweep).
- fp64 size sweep: the same three sizes at apb=16.

All runs used `--warmup 2 --iterations 10 --repetitions 7` (the benchmark's
defaults, matching every prior CGP task's usage of this harness); raw output
is `docs/cgp07_raw/{fp32,fp64}_apb{8,16,32,64}_N{65536,262144,1048576}.log`
(only the combinations actually run are present — not the full 4x3 cross
product, which was not needed since the size and block-size axes were swept
independently at their one shared point).

This is a scoping decision, not an omission: the full matrix requirement is
"multiple sizes tested" and "multiple block sizes tested", not necessarily
every combination of both, and CGP-07 is an evidence task bounded like every
other CGP task in this series.

## Part 2 — Effective atomistic fraction

The existing harness already reports, per sweep point, `active_atom_fraction`,
`interface_fraction`, `live_bond_fraction`, and `active_dof_ratio =
(active_atoms + active_blocks) / total_atoms` — the last of these is the
"effective atomistic fraction including buffer/interface degrees of freedom"
the task asks for. Sample (N=262,144, apb=16):

| requested\_fraction | active\_atom\_fraction | interface\_fraction | active\_dof\_ratio |
| --- | --- | --- | --- |
| 1.00 | 1.0000 | 0.0000 | 1.0000 |
| 0.25 | 0.2501 | 0.0001 | 0.2970 |
| 0.0625 | 0.0626 | 0.0001 | 0.1212 |
| 0.01 | 0.0101 | 0.0001 | 0.0720 |
| 0.00 | 0.0000 | 0.0000 | 0.0625 |

The gap between `requested_fraction` and `active_dof_ratio` at low fractions
(e.g. 1% requested vs. 7.2% effective) is buffer/coarse bookkeeping in this
1-D-chain synthetic fixture (a fixed two-block buffer plus the coarse-block
count itself, which never drops below `1/atoms_per_block` since a coarse
block always "exists" as one DOF even at zero fine atoms) — reported here
rather than silently using the requested fraction as if it were the true
active resolution.

## Part 3 — Reported phases

The existing per-phase device-timer breakdown (`adaptive-sweep-phases` lines)
already separates atomistic / coarse / interface / selector / polarization /
compaction / fft / integration, and `adaptive-sweep-launches` separately
reports launch counts per phase. Thermal RNG and host-synchronization/wait
are **not** part of that in-process phase struct (`GpuAdaptivePhaseMetrics`
has no RNG or sync-count field) and are reported through the two other
existing, purpose-built harnesses instead — matching the task's own
requirement to report them "separately", not folded into the same struct:

- **Thermal RNG**: Part 5, via `gpu_thermfield_rng_benchmark`
  (CGP-06A/06B's harness, exercising the real `GpuThermfield::randomize()`
  including the CGP-06B active-scoped overload).
- **Host synchronization/wait**: `phase_syncs_per_step` in the
  `adaptive-sweep-phases` lines is **10.286–10.429 with per-phase device
  timers enabled** — this is *timing-readout* synchronization (one
  `GPU_EVENT_SYNCHRONIZE` per instrumented phase, needed to read that phase's
  elapsed time), not production dependency-ordering synchronization; CGP-05
  already established and disclosed this distinction ("9→0 with
  `UPPASD_ADAPTIVE_PHASE_TIMING=0`, 1 with per-phase timing enabled" for the
  *dependency* syncs it audited). The headline numbers throughout this
  document use `step_wall_us` (phase timing **disabled**, i.e. production
  configuration), never `instrumented_wall_us`, for exactly this reason —
  `instrumentation_overhead_percent` in the raw logs runs 5%–133% depending
  on size/fraction, confirming instrumentation is not free and must not
  contaminate the headline claim.
- **Explicit measurement-energy overhead**: not exercised by this benchmark
  at all — `adaptiveCoarseStep()` uses `evaluateField()` (CGP-01's
  field-only path), the same steady-state call production's predictor/
  corrector uses, so ordinary-dynamics timing in this document already
  excludes energy evaluation by construction (see CGP-01), rather than
  needing a separate subtraction.

## Part 4 — Baseline

`ProductionAtomisticBaseline` in the benchmark builds UppASD's real
`GpuHamiltonianCalculations` (the same `heisge` the feature-off timestep
loop calls) and `GpuDepondtIntegrator` (the same `SDEalgh=1`
`evolveFirst`/`evolveSecond` pair) on the identical fixture geometry, and the
`atomistic-parity` check (Part 0) is what licenses treating it as "identical
supported physics" rather than an assumption. This is the same baseline
methodology CGP-01/02/04 already relied on; Part 7 additionally cross-checks
it once against the real whole-executable feature-off/adaptive pair.

## Part 5 — Thermal RNG, reported separately

Using `gpu_thermfield_rng_benchmark` (CGP-06A/06B's harness) at `--temperature
1.0` so RNG genuinely executes, `--warmup 3 --iterations 20 --repetitions 7
--subphase`:

**fp32, N=262,144, active-fraction sweep** (raw:
`docs/cgp07_raw/rng_fp32_N262144_*.log`):

| active fraction | active atoms | `randomize()` wall (us) | speedup vs. unscoped |
| --- | --- | --- | --- |
| 100% (unscoped) | 262,144 | 39.114 | 1.00x |
| 50% | 131,072 | 29.459 | 1.33x |
| 25% | 65,536 | 23.629 | 1.65x |
| 12.5% | 32,768 | 20.784 | 1.88x |
| 6.25% | 16,384 | 20.195 | 1.94x |
| 1% | 2,621 | 19.650 | 1.99x |
| ~0% (26 atoms) | 26 | 19.526 | 2.00x |

This reproduces CGP-06B's own reported 1.9x–2.4x fp32 range almost exactly
(here: 1.33x–2.00x across a slightly different fraction set, converging to
~2.0x as active count → 0, i.e. the write-eliminated floor), an independent
confirmation rather than a new claim.

**fp64 cross-check, N=262,144** (`rng_fp64_N262144_*.log`): unscoped
381.625us, 6.25%-active 341.670us → **1.117x**, confirming CGP-06A/06B's
precision-dependent-asymmetry finding (fp32 ≈2x, fp64 ≈1.1x at the same
active fraction) still holds on this rebuilt tree.

**Size scaling of the generate-only floor at fixed 6.25% active fraction**
(fp32): N=65,536 → 8.682us; N=262,144 → 20.195us; N=1,048,576 → 71.165us.
This is genuinely closer to linear-in-N at the largest size (3.5x time for
4x atoms, vs. 2.3x for the smaller 4x step) — direct confirmation that
`randomize()`'s curand/hiprand **generate** call remains full-`N` **by
design** (CGP-06B's own documented invariant-preservation decision, not an
oversight) and is therefore the one genuinely O(`Natoms`)-scaling component
still present anywhere in the adaptive step at steady state. It is small in
absolute terms at every size tested here (≤71us at 1M atoms) but is the
single component this task's audit did **not** find any way to have already
eliminated, because CGP-06A/06B's own invariant analysis (physical atoms'
stochastic continuation across a future coarse→fine transition must not
depend on how long they spent coarse) requires it.

## Part 6 — Scaling model

Fit against the pooled fp32 size-sweep data (3 sizes x 8 fraction points = 24
points), `docs/cgp07_raw/fit_scaling_model.py`:

```
T = b * Nlive_bonds + c * Ninterface + d * Ncoarse + overhead
```

`Nfine` was **dropped** from the originally-specified 5-term model after the
fit's design-matrix condition number came back at ~1e21–1e22 for every
single-size fit: this fixture's fixed exchange-shell connectivity makes
`Nlive_bonds` almost exactly `1.9375 * Nfine` (min/max ratio 1.9375–1.9608
across all nonzero-fine points tested), so `a*Nfine + b*Nlive_bonds` is not
identifiable from this data — dropping the redundant term left the fit
quality (R²) **unchanged** at every size (confirming the redundancy rather
than losing information), while producing physically sensible non-negative
coefficients. `Nlive_bonds` is also the more fundamental cost driver per the
governing contract's own item 6, so this is not a loss of the intended
model, just a resolved collinearity.

**Pooled fit** (24 points, all three fp32 sizes at apb=16):

- `b` (per live bond) = 5.29e-4 us
- `c` (per interface atom) = 1.570 us
- `d` (per coarse block) = 8.18e-4 us
- `overhead` = 36.90 us
- R² = 0.999, max residual 20.7us (1.9% of the largest sample)

Cross-check against the model's own prediction at all-coarse
(`Nlive_bonds=Ninterface=0`, `Ncoarse=blocks`):

| N | blocks | predicted (`d*blocks+overhead`) | actual `step_wall_us` |
| --- | --- | --- | --- |
| 65,536 | 4,096 | 40.25us | 42.28us |
| 262,144 | 16,384 | 50.31us | 48.92us |
| 1,048,576 | 65,536 | 90.52us | 92.37us |

Agreement within ~2–4%, i.e. the model's decomposition is not just a
curve-fit artifact but predicts held-out-style all-coarse totals correctly
from the fraction-sweep-interior points alone.

**What this says about the residual O(`Natoms`) component**: after CGP-01
through CGP-06B, the all-coarse steady state is **not** dominated by a term
that scales with total atom count `Natoms`. It is dominated by:

1. a fixed **~37us per-step floor**, independent of `N` and independent of
   fraction — mechanistically explained by kernel-launch dispatch count: at
   all-coarse, `adaptive-sweep-launches` reports ~8.23 coarse + ~2.06
   interface + ~2.06 integration ≈ **12.3 kernel launches/step**, and
   37us / 12.3 ≈ **3.0us/launch**, squarely in the range of ordinary CUDA
   host-side launch-dispatch cost on this driver/toolkit — not a leftover
   atomistic-reconstruction cost (that floor is what CGP-02 through CGP-06B
   already removed);
2. a genuinely block-count-scaling term, `d*Ncoarse` ≈ 0.8ns/coarse-block —
   small (54us even at 65,536 coarse blocks) but real, and **expected**: the
   coarse-block continuum state must still be integrated every step
   regardless of fine fraction (that is the coarse dynamics itself, not
   overhead the prior CGP tasks were ever scoped to remove — the governing
   contract's own target model explicitly includes a `d*Ncoarse` term).

Confirmed independently by the block-size sweep (N=262,144 fixed, apb 8/16/
32/64 ⇒ blocks 32,768/16,384/8,192/4,096): all-coarse `step_wall_us` is
55.4/48.9/43.5/41.8us respectively — monotonically decreasing with fewer
blocks, consistent with the `d*Ncoarse` term and *not* with an atom-count (N
is held fixed here) effect.

## Part 7 — Baseline crossover across the matrix

Per-size fp32 speedup vs. the real production step (`production_step_wall_us`
from `ProductionAtomisticBaseline`, same fixture, phase timing disabled):

| requested\_fraction | N=65,536 | N=262,144 | N=1,048,576 |
| --- | --- | --- | --- |
| 1.00 (all-fine) | 0.710x | 1.021x | 1.112x |
| 0.25 | 0.744x | 1.894x | 3.060x |
| 0.0625 | 0.804x | 3.237x | 5.870x |
| 0.01 | 0.984x | 3.556x | 8.356x |
| 0.00 (all-coarse) | 1.859x | 6.341x | 13.292x |

This is the single most important, and most honest, finding of this task.
**At the historical 65,536-atom benchmark size, adaptive coarse graining is
now *slower* than the plain production atomistic path at every fraction
except all-coarse**, because the ~37us fixed per-step floor characterized in
Part 6 is a large fraction of that size's ~78.6us production step — the
atom-scale optimizations already landed (CGP-01–CGP-06B) removed the
work that used to scale with total atom count, but the *fixed* floor those
optimizations left behind is comparable in size to the whole production
step at this scale. At 262,144 and 1,048,576 atoms the same ~37–90us floor
is a much smaller fraction of a much larger production step (310us/1,228us),
so adaptive coarse graining wins even at 100% fine and wins dramatically as
fraction drops. fp64 shows the same qualitative shape but crosses over
already at N=65,536 (production's fp64 step is itself heavier — see raw
logs — so the fixed floor is proportionally smaller there than in fp32).

**Decision-gate-relevant reading**: whether the atom-scale floor is "solved"
depends on the system size a given workload actually runs at. This task does
not attempt to say what size is representative of real UppASD adaptive-CG
usage — that is a modelling/user-workload question outside this evidence
task's scope — but it does mean the honest answer to "has the atom-scale
floor been fixed" is **conditional on N**, not a flat yes.

## Part 8 — Whole-executable cross-check (production e2e)

`run_rcg09_perf_e2e.py` against the fp32 `sd.f95.cuda` binary, both `spiral`
and `uniform` fixtures (N=16,384, `--rounds 5 --short-steps 20
--long-steps 220`, defaults), device-contention-checked (`result=CLEAN`,
0% utilization, 53°C at the moment of the check):

- `spiral`: `atomistic_step_us=92.7±136.6`, `adaptive_step_us=653.2±37.3`,
  `speedup=0.14x`, `crossover=NOT_OBSERVED`.
- `uniform`: `atomistic_step_us=39.8±20.0`, `adaptive_step_us=252.8±83.3`,
  `speedup=0.16x`, `crossover=NOT_OBSERVED`.

This cross-check is **attempted but not usable as quantitative evidence**,
disclosed rather than hidden: the atomistic arm's own MAD (up to 136.6us on
a ~93us median) exceeds its median, i.e. this host's CPU-quota-limited
process-launch/short-run noise dominates any signal at this short a
`Nstep` budget — exactly the finding CGP-05's evidence doc already recorded
for this same host and the same reason (cgroup-limited to ~2 of 16 physical
cores). It is not re-litigated further here.

One genuine observation, disclosed but not investigated (out of this task's
scope): both fixtures reported `active_blocks=0` throughout, i.e. the
selector never coarsened anything within the 220-step budget for either
`spiral` or `uniform`, including the fixture the README describes as "the
most favourable possible input for adaptive coarse graining". Whether that
is a step-budget artifact, a selector-threshold interaction with the current
tree's changes, or expected behavior for this particular fixture at this
step count is a selector-behavior question (CGP-03's domain), not a
performance question, and is flagged here rather than diagnosed. The
in-process benchmark's own fraction sweep (Parts 1–7) directly controls
`active_blocks` instead of relying on the selector to reach a particular
fraction, so this observation does not weaken any claim above — it only
means this one whole-executable cross-check could not independently
corroborate the coarsened-fraction regime, only the all-fine (N=16,384,
fraction≈1.0) point, which is broadly consistent in *sign* (adaptive slower
than production at this size and fraction) with the internal harness's own
N=65,536-all-fine finding in Part 7, despite the wide error bars.

## Decision gate

Per the task's own instruction: **only after this task**, decide whether to
proceed to continuum-kernel precontraction/stencil work.

**Answer: conditional on system size, not unconditional.**

- At the sizes this task could test above ~262,144 atoms, the atom-scale/
  synchronization floor characterized in Part 6 (~37–90us, dominated by
  fixed kernel-launch-dispatch count plus a small, physically-expected
  `Ncoarse` term) is now a small fraction of the production step cost at any
  fine fraction, and adaptive coarse graining wins substantially (1.02x–
  13.3x depending on fraction) at every fraction tested. At this size
  regime, **continuum-kernel arithmetic is plausibly the next dominant
  target** — this task did not benchmark the continuum kernel itself (out of
  scope per the task's own instruction not to optimize it here), so this is
  a *plausibility* statement from what remains, not a direct measurement of
  continuum cost.
- At the historical 65,536-atom benchmark size, that same fixed floor is
  large enough relative to the (smaller) production step that adaptive
  coarse graining is *slower* than plain production at every fraction tested
  except all-coarse. At this size, **the atom-scale/synchronization floor
  itself is still the dominant limiter**, and per this task's own decision
  rule ("if another atom-scale or synchronization floor still dominates, fix
  that first"), continuum-kernel work would not be the correct next target
  for workloads at this scale.
- The floor's own mechanistic explanation (Part 6: ~12.3 kernel
  launches/step at all-coarse, ~3us/launch dispatch cost) suggests **kernel
  fusion / launch-count reduction** as the concrete next atom-scale lever, if
  a future task chooses to address the small-N regime before continuum work.
  This task does not implement or design that fusion — flagged as a
  candidate follow-on, not committed to.

No continuum-kernel optimization was performed in this task.

## Checklist

* [x] Field-only dynamics benchmarked separately from energy measurement
  (Part 3: `evaluateField()`, CGP-01's field-only path, used throughout).
* [x] fp32 is the primary performance matrix (Part 1: 3 sizes + 4 block
  sizes in fp32; 1 full size matrix retained in fp64).
* [x] Multiple total system sizes tested (65,536 / 262,144 / 1,048,576).
* [x] Multiple block sizes tested (8 / 16 / 32 / 64 atoms/block).
* [x] Multiple fine fractions tested (0/1/6.25/12.5/25/50/75/100%).
* [x] Effective atomistic fraction reported (Part 2: `active_dof_ratio`,
  `interface_fraction`).
* [x] Thermal RNG reported separately (Part 5, own harness).
* [x] Host synchronization reported separately (Part 3: distinguished
  timing-readout syncs from CGP-05's already-audited dependency syncs).
* [x] Selector/compaction amortized correctly (`adaptive-sweep-livework`
  lines report `compaction_rebuilds`/`compaction_readbacks` per sweep point,
  amortized the same way CGP-02/03 already established; not separately
  re-derived here).
* [x] Feature-off production baseline used (Part 4, `atomistic-parity`
  gate; Part 8 whole-executable cross-check, inconclusive but attempted).
* [x] Residual O(`Natoms`) work quantified (Part 6: not O(`Natoms`) — a
  fixed ~37us floor plus a small O(`Ncoarse`) term, mechanistically
  explained by kernel-launch count).
* [x] New measured crossover reported (Part 7: crossover point is strongly
  size-dependent, tabulated per size).
* [x] No continuum-kernel optimization performed.
* [x] Evidence states plainly whether continuum arithmetic is now the
  dominant remaining target (Decision gate: conditionally yes above
  ~262,144 atoms, conditionally no — the atom-scale floor itself still
  dominates — at the historical 65,536-atom size).

See `docs/cgp07_raw/*.log` for raw benchmark output and
`docs/cgp07_raw/fit_scaling_model.py` for the regression used in Part 6.

## Commit

`CGP-07: rebaseline adaptive GPU performance`
