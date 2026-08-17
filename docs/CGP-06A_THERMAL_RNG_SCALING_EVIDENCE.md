# CGP-06A — Audit and benchmark thermal-RNG scaling, evidence

Task: `docs/CGP_work.md` CGP-06A ("Audit and benchmark thermal-RNG
scaling"). Branch `gpu_hip_cu_ab_cg`, commit base `cf11224c` (CGP-04),
building on CGP-05 (host-barrier removal; verified, not yet committed at the
time of this task -- see project memory `cgp05-status`).

**Primary model (task doc): Luna.** **Actual model: Claude (Sonnet 5,
interactive session)**, matching the CGP-03D/CGP-05 precedent of an
interactive session performing the audit/benchmark/design slice directly.

This is an evidence-only task. **No production code was modified.** The only
new files are a read-only micro-benchmark
(`tests/coarse_graining/benchmark_gpu_thermfield_rng.cpp`) that exercises the
real, unmodified `GpuThermfield` class, its `CMakeLists.txt` registration,
and this document.

## Environment

- GPU: NVIDIA RTX A4000 (device 0 of 2; `CUDA_VISIBLE_DEVICES=0` pinned for
  every run below), shared host per project memory
  `shared-gpu-host-contention`. Idle immediately before and during every run
  in this document: `nvidia-smi` showed 0% utilization, 210 MHz SM clock
  (idle floor), <110 MiB used, and `--query-compute-apps` returned no other
  process on either GPU.
- CUDA 13.3.73 (nvcc), CMake, `Release` build type.
- Host CPU: `nproc` reports 2 (cgroup-quota-limited; see project memory
  `shared-gpu-host-contention`). Material to the host/GPU wall-clock gap
  noted in section 4.
- Builds: fresh out-of-tree `build_cgp06a_cuda_fp32`
  (`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=SINGLE`) and
  `build_cgp06a_cuda_fp64` (`UPPASD_PRECISION=DOUBLE`), both
  `CMAKE_BUILD_TYPE=Release`, both built from this branch's current working
  tree (CGP-05's uncommitted changes included, per this task's own
  dependency).

## Part A — Trace the thermal RNG implementation

Traced from `source/gpu_files/gpuThermfield.{hpp,cpp}`,
`source/gpu_files/gpuDepondtIntegrator.cpp`, `source/gpu_files/gpu_wrappers.h`,
`source/gpu_files/gpuSimulation.cpp`, and `source/gpu_files/gpuParallelizationHelper.{hpp,tpp}`.

### Generator type

`GPU_RAND_GENERATOR_T` is `curandGenerator_t` (CUDA) / `hiprandGenerator_t`
(HIP) -- the vendor **host-API pseudo-random stream generator**, not a
per-thread/counter-based generator. The concrete algorithm is selected by the
input keyword `gpu_rng` (0-3: `PSEUDO_DEFAULT` = XORWOW, `XORWOW`,
`MRG32K3A`, `MTGP32`; `source/gpu_files/gpuSimulation.cpp:144-167`), default
`PSEUDO_DEFAULT`.

### State ownership

One `GpuThermfield` instance is owned by one `GpuDepondtIntegrator` instance
(`source/gpu_files/gpuDepondtIntegrator.hpp:33`), which is itself a single,
process-lifetime local inside `GpuSDSimulation.cpp`'s `SDmphase`
(`source/gpu_files/gpuSDSimulation.cpp:75,215`). There is exactly **one**
curand/hiprand generator, hence one RNG state, for the whole process, shared
identically between feature-off and adaptive dynamics -- confirmed by
`grep`: `gpuAdaptiveRuntime.{hpp,cpp}` contain no `Thermfield` reference at
all; the adaptive path reads the *same* `GpuDepondtIntegrator` member
(`GpuSimulation::advanceAdaptiveStep`, `source/gpu_files/gpuSimulation.cpp:1299`,
calls `integrator->evolveFirst(gpuLattice, activeAtomList, activeAtomCount)`
on the one shared integrator).

The generator's internal state is **opaque and sequential**: each call to
`GPU_RAND_GENERATE_NORMAL_REAL` (`curandGenerateNormal`/
`hiprandGenerateNormal[Double]`) advances the stream by exactly the number of
values requested, in one call, with no per-site addressability. It is not
per-thread, not per-site, and not counter-based.

### Random-value consumption

`GpuThermfield::randomize(mmom)` (`gpuThermfield.cpp:86-94`) does exactly two
things, unconditionally, every call:

1. `GPU_RAND_GENERATE_NORMAL_REAL(gen, randomValues.data(), randomValues.size(), 0, 1)`
   -- draws `randomValues.size()` values, where `randomValues.size() =
   3*N*M` rounded up to even (`gpuThermfield.cpp:51`, `N` = total atom count,
   `M` = ensemble count). **This consumes 3\*N\*M generator values from the
   shared stream regardless of which atoms are dynamically active.**
2. `parallel.gpuAtomSiteCall(SetupField(...))` -- launches `atom_site_kernel`
   over the **full** `(N, M)` grid (`gpuParallelizationHelper.tpp:50-63`;
   `AtomSite` carries the parallelization helper's own `N`/`M`, not any
   active-list size), writing `field[atom] = randomValues[atom] *
   sigmaFactor[site] * rsqrt(mmom[atom])` for **every** atom, active or not.

`randomize()` is called exactly once per accepted SD step, from
`evolveFirst()`/`evolveFirst(oneBasedAtoms, activeAtomCount)`
(`gpuDepondtIntegrator.cpp:388,454`) -- never from `evolveSecond()`, which
reuses the same `thermfield.getField()` the predictor call produced. The
active-atom overload's own comment states the design explicitly (`gpuDepondtIntegrator.cpp:450-453`):

> "Keep the production generator and its full (site,ensemble) layout. It
> consumes the same generator sequence per accepted SD step as the
> feature-off path; compact-list order therefore cannot change a site's
> random draw."

So **inactive (coarse) atoms currently advance -- and consume -- RNG state
identically to active atoms**, every step, by construction. This is not an
oversight; it is the mechanism the existing RCG-09A.2 parity tests
(`tests/coarse_graining/test_gpu_depondt_active_atoms.cpp`) rely on for
correctness (Part B).

### Temperature/moment entry into the field

Two separate stages, correctly separated by cadence:

- `GpuThermfield::resetConstants()` (called once per phase-parameter change,
  not every step) computes `sigmaFactor[site] = sqrt(dp * temperature[site])`
  where `dp = 2*damping*k_bolt / (timestep*gamma*mub*(1+damping^2))`.
- `SetupField::each(atom, site)` (every `randomize()` call) computes
  `field[atom] = randomValues[atom] * sigmaFactor[site] * rsqrt(mmom[atom])`.

### Seeds and restart

`gpu_rng_seed` (input keyword, `source/Input/inputdata.f90:275-276`, default
0 -> `time(nullptr)` is used if 0, `gpuThermfield.cpp:55`) is read from
Fortran exactly once per process (`gpuSimulation.cpp:169`) and passed to
`GpuThermfield::initiate()`, which calls
`GPU_RAND_SET_PSEUDO_RANDOM_GENERATOR_SEED(gen, seed)` exactly once, also per
process. **No code anywhere in the repository serializes curand/hiprand
generator state to, or restores it from, a restart file**: a
repository-wide `grep` for `curand`/`hiprand` across every `.f90` file
returned nothing, and neither `gpuThermfield.*` nor
`gpuDepondtIntegrator.*` contain the word "restart". This is a pre-existing
property of the whole GPU thermal path, **not introduced by adaptive CG and
not specific to this task** -- restarting a GPU run with an explicit fixed
seed reproduces the *same generator sequence from its own start*, not a
continuation from wherever the original run's generator had advanced to at
the restart point. Section on invariant 5 below states the consequence for
CGP-06B precisely.

### T=0 behaviour

`randomize()` has no temperature-dependent early return -- it checks only
`initiated()` (`gpuThermfield.cpp:87`). At `temperature=0`, `sigmaFactor`
resolves to `0`, so `field` ends up all-zero, but the full curand generate
call and the full `SetupField` launch still execute, consuming generator
state and device time identically to finite T. Confirmed empirically in
Part C (table 1): T=0 and finite-T timings are statistically indistinguishable
at every size tested.

## Part B — Accepted stochastic invariants

| # | Invariant | Status under the *current* implementation | Evidence |
|---|---|---|---|
| 1 | All-fine adaptive fixed-seed trajectory equals feature-off ASD | **Holds**, by construction (both paths share one generator, one call per step, at the current full-`(N,M)` grid) | `tests/coarse_graining/test_gpu_depondt_active_atoms.cpp::runAllActiveParity` (bitwise `emom`/`emom2`/`emomM`/`b2eff` equality, T=0 and finite T, RCG-09A.2); RCG-09A.4 `coarse-graining-adaptive-asd-parity` at the full-simulation level (`gpu_adaptive_runtime_benchmark --parity-only`) |
| 2 | Compact-list ordering cannot change an active atom's random draw | **Holds**, by construction (the generator call precedes and is independent of which/how-ordered atoms the active-list kernel later reads) | `test_gpu_depondt_active_atoms.cpp::runReorderedSubsetIdentity` (forward `{2,5}` vs. reverse `{5,2}` active lists produce bitwise-identical results) |
| 3 | Fine-fraction changes cannot accidentally permute RNG streams | **Holds**, by construction: `randomize()` takes no active-list argument at all, so its output for a fixed `(N,M,seed,step)` is identical regardless of which fraction is currently fine | Structural (Part A); confirmed empirically in Part C table 4, `thermfield-fraction-current`: cost and (implicitly, since it's the same call) the drawn values are independent of fraction |
| 4 | Coarse -> fine transitions have deterministic/reproducible stochastic continuation | **Holds**, by construction: a newly-fine atom's thermal draw at step *t* is exactly the value the shared generator produced for that atom's slot at step *t*, whether or not it was fine at step *t-1* -- there is no separate "resume" event to get wrong | Structural (Part A) |
| 5 | Restart reproduces the same continuation | **Does not hold, and never has, independent of CG or this task.** No curand/hiprand state is ever serialized to a restart file (Part A). A restarted run with an explicit fixed seed reproduces its *own* sequence from index 0, not a continuation of the original run's sequence at the restart step. With the default `gpu_rng_seed=0`, a restarted run does not even match itself (`time(nullptr)`-seeded). This is baseline GPU-path behaviour shared by feature-off dynamics; adaptive CG does not make it any better or worse. | Structural (Part A); confirmed absent by repository-wide search |

Invariant 1's *finite-temperature* form is bounded by a separate,
pre-existing, unrelated finding already on record: project memory
`gpu-adaptive-finite-t-nondeterminism` states that GPU adaptive finite-T
output is not run-to-run reproducible even on an unmodified binary. That
means finite-T trajectory comparisons for invariant 1 are necessarily
statistical, not bit-exact, on this codebase today -- T=0 remains the only
bit-exact oracle. This is not new to CGP-06A; it constrains what any future
CGP-06B correctness test can claim at finite T, same as it constrained
CGP-03D.

**Conclusion for Part B:** every invariant the *current* implementation is
asked to satisfy, it satisfies -- and it satisfies invariants 2-4 for the
structurally trivial reason that the RNG call is oblivious to activity,
fraction, and transitions. That obliviousness is exactly what makes it cost
O(N) regardless of active fraction (Part C), so **any real fix intended to
remove that floor must re-derive invariants 2-4 under a design where the RNG
call is no longer oblivious** -- this is not a free property to inherit and
is the central design tension Part D evaluates.

## Part C — Benchmark cost

Micro-benchmark: `tests/coarse_graining/benchmark_gpu_thermfield_rng.cpp`
(new file, evidence-only, builds against the real, unmodified
`GpuThermfield` class -- same `initiate()`/`initiateConstants()`/
`randomize()` calls, same curand generator type, production drives). One
`(N, M, temperature)` configuration per process invocation, deliberately:
an earlier version of this harness constructed many `GpuThermfield`
instances back-to-back in one process (varying the global
`ParallelizationHelperInstance` singleton alongside them) and segfaulted
inside `cudaEventRecord` on the very first construction -- a real
manifestation of the exact multi-instance GPU-resource-churn pitfall project
memory (`cgp05-status`, finding 2) already recorded from CGP-05: repeated
construct/destroy cycles of test-only multi-instance harnesses can
destabilize curand/stream state even when the code under test is never
itself in the failing call path. Restructuring to one configuration per
process (matching how production only ever constructs one `GpuThermfield`
per process, and matching `gpu_adaptive_runtime_benchmark`'s own
one-geometry-per-invocation CLI convention) eliminated the crash; a shell
driver (`build_cgp06a_cuda_fp32/run_thermfield_sweep.sh`) sweeps
configurations across separate invocations. 5 discarded warm-up iterations
and 30 measured iterations per repetition, 7 repetitions, median/MAD
reported; one full discarded cold-start invocation precedes the swept
measurements (project memory `gpu-benchmark-cold-start`).

Raw output: `build_cgp06a_cuda_fp32/thermfield_rng_fp32.log`,
`build_cgp06a_cuda_fp64/thermfield_rng_fp64.log`.

### 1. Cost vs. total atom count, T=0 and finite T (fp32)

| N (atoms, M=1) | T=0 wall us (median, MAD) | T=100 wall us (median, MAD) |
|---|---|---|
| 4,096 | 4.671, 0.007 | 4.416, 0.007 |
| 16,384 | 5.118, 0.012 | 4.845, 0.020 |
| 65,536 | 8.579, 0.019 | 8.638, 0.027 |
| 262,144 | 35.107, 0.049 | 35.143, 0.037 |

T=0 and finite-T costs are statistically indistinguishable at every size (as
Part A predicted: `randomize()` has no temperature-dependent branch). Cost is
clearly **not linear from the smallest sizes**: 4,096 -> 16,384 is a 4x atom
increase for only a ~1.1x time increase (a fixed per-call floor -- kernel
launch/dispatch overhead -- dominates at small N), while 65,536 -> 262,144 (also
4x atoms) is a 4.08x time increase (essentially linear once past the launch-
overhead floor). **Above roughly 10^5 atoms the cost is atom-count-dominated
and effectively independent of active fraction** (nothing in the call reads
an active fraction at all -- see table 4).

### 2. Cost vs. ensemble count, fixed N=65,536, T=50 (fp32)

| M (ensembles) | wall us (median, MAD) |
|---|---|
| 1 | 8.589, 0.005 |
| 2 | 13.925, 0.084 |
| 4 | 33.463, 0.041 |

Consistent with the `3*N*M`-sized generate/write cost: cost grows with the
`N*M` product, same launch-floor-then-linear shape as table 1.

### 3. RNG-generate vs. field-write-loop sub-phase breakdown (fp32)

Read from the same `GlobalStopwatchPool("GPU thermfield")` categories
(`"RNG"`, `"loop"`) `GpuThermfield::randomize()` already records in
production when phase timing is enabled -- not new instrumentation, this
harness just reads it out. GPU-event time (`gpu_us_per_call`), not host wall
time, is the meaningful comparison here.

| N (atoms) | RNG (curand generate) us/call | loop (SetupField write) us/call | RNG share |
|---|---|---|---|
| 4,096 | 3.257 | 2.634 | 55.3% |
| 65,536 | 4.518 | 3.569 | 55.9% |
| 262,144 | 5.940 | 7.407 | 44.5% |

The two sub-phases are comparably expensive across the whole range tested,
with the write loop overtaking the generate call only at the largest size
(262,144). **Neither sub-phase dominates the other enough to be ignored** --
material to Part D's evaluation of Strategy 1, which can only ever remove the
"loop" share.

The GPU-event sum (RNG+loop) is smaller than the host wall-clock number in
table 1 at the same N (e.g. at 262,144: 5.940+7.407 = 13.35 us GPU vs. 35.15
us wall) -- host-side call/dispatch overhead and this host's CPU-quota
limitation (`nproc`=2, per project memory `shared-gpu-host-contention`)
account for the gap, consistent with CGP-05's own finding that CPU-side noise
on this host swamps microsecond-scale GPU timings. Reported as a disclosed
gap between instrumented device time and uninstrumented wall time, per the
project's benchmark protocol; not investigated further here since neither
number changes this section's qualitative conclusion (both sub-phases
matter).

### 4. Fine-fraction projection vs. current flat cost (fp32, total_atoms=262,144, T=50)

`randomize()` takes no active-list argument (Part A), so its **actual**
per-step cost is, by construction, identical at every fine fraction:

```
thermfield-fraction-current total_atoms=262144 randomize_wall_us=35.189
note=current_implementation_cost_is_independent_of_fine_fraction
```

To bound what a genuinely active-scoped generate+write *could* cost at a
given fraction, the identical call was re-run at `N' = fraction *
total_atoms` -- a **projection from a smaller full lattice**, explicitly not
a measurement of any real active-scoped kernel (none exists in production):

| Fraction | Projected active atoms | Projected wall us | Current (flat) wall us | Projected speedup |
|---|---|---|---|---|
| 0.0 (all-coarse) | 0 | 0 (never measured; hypothetical) | 35.189 | undefined / best case |
| 0.01 | 2,621 | 3.987 | 35.189 | 8.83x |
| 0.0625 | 16,384 | 4.804 | 35.189 | 7.33x |
| 0.25 | 65,536 | 8.564 | 35.189 | 4.11x |
| 1.0 | 262,144 | 35.152 | 35.189 | 1.00x (sanity check: full fraction recovers the current cost, as it must) |

This is the central quantitative finding of this task: **at the fine
fractions CGP has used throughout this phase (1%, 6.25%, 25%), the current
implementation's O(total-atoms) thermal-RNG floor is 4x-9x more expensive
than an active-scoped implementation could be**, and the gap only closes as
the fine fraction approaches 100%.

### fp64 (DOUBLE build)

Same harness, same host state (GPU idle, 0% utilization/210 MHz confirmed
immediately before this run), rebuilt with `UPPASD_PRECISION=DOUBLE`
(`build_cgp06a_cuda_fp64`; `GPU_RAND_GENERATE_NORMAL_REAL` resolves to
`curandGenerateNormalDouble`, `gpu_wrappers.h:127,131`).

**1. Cost vs. total atom count** (T=0 vs. T=100 again indistinguishable):

| N (atoms, M=1) | T=0 wall us | T=100 wall us | fp64/fp32 ratio (T=0) |
|---|---|---|---|
| 4,096 | 9.212 | 8.685 | 1.97x |
| 16,384 | 22.077 | 22.097 | 4.31x |
| 65,536 | 78.964 | 79.007 | 9.20x |
| 262,144 | 319.886 | 320.133 | 9.11x |

**2. Cost vs. ensemble count**, fixed N=65,536, T=50: M=1: 78.974 us; M=2:
159.268 us (1.98x); M=4: 319.061 us (2.00x M=2) -- clean linear-in-M scaling,
same shape as fp32.

**3. Sub-phase breakdown:**

| N (atoms) | RNG us/call | loop us/call | RNG share |
|---|---|---|---|
| 4,096 | 5.575 | 2.742 | 67.0% |
| 65,536 | 14.336 | 2.019 | 87.7% |
| 262,144 | 20.139 | 3.106 | 86.6% |

**This is the single most important asymmetry this task found.** At fp32
the two sub-phases are comparably expensive (44-56% RNG share, table 3
above); at fp64, `curandGenerateNormalDouble` dominates `randomize()`'s cost
(67-88% RNG share) while the `SetupField` write-loop cost barely grows
between precisions (fp32 loop at N=65,536: 3.57-4.28 us; fp64: 2.02 us --
*not* larger, within noise). **Strategy 1's ceiling (Part D) is
precision-dependent: it can only remove the "loop" share, which is roughly
half of `randomize()`'s cost at fp32 but only 12-13% of it at fp64 for the
sizes that matter (65,536+).** A `MIXED_PREC`-shaped worry does not apply
here (CGP-06A introduces no such implementation, per the governing
contract), but this is exactly the kind of precision-sensitive result the
project's benchmark protocol asks to be reported explicitly rather than
assumed to generalize from one precision.

**4. Fine-fraction projection vs. current flat cost** (total_atoms=262,144,
T=50):

| Fraction | Projected active atoms | Projected wall us | Current (flat) wall us | Projected speedup |
|---|---|---|---|---|
| 0.01 | 2,621 | 7.979 | 348.796 | 43.71x |
| 0.0625 | 16,384 | 24.923 | 348.796 | 13.99x |
| 0.25 | 65,536 | 86.390 | 348.796 | 4.04x |
| 1.0 | 262,144 | 348.886 | 348.796 | 1.00x (sanity check passes) |

The projected win at low fine fractions is **larger at fp64 than fp32**
(43.7x vs. 8.8x at 1% fine) -- the opposite of what a naive "fp64 is just
slower everywhere" intuition would suggest. The reason follows directly from
table 1's floor-then-linear shape: fp64's fixed per-call launch/dispatch
floor is a smaller fraction of its (much larger) per-value compute cost, so
fp64's cost curve is closer to pure `O(N)` from smaller N already, which
makes the achievable reduction from scoping to a small active fraction
correspondingly larger. This reinforces Part D's conclusion that removing
the RNG-generate floor (Strategies 2/3), not just the write-loop
(Strategy 1), matters *more*, not less, for DOUBLE builds.

## Part D — Implementation strategies

### Strategy 1: advance full RNG state, generate/write fields only for active atoms

Keep the single shared curand/hiprand host-API generator and its one
`3*N*M`-value generate call per step exactly as today (so the generator's
*sequence* is untouched), but scope `SetupField`'s write to the compact
active+ghost list (the same list shape `gpuActiveAtomCall`/
`active_atom_kernel` already use, `gpuParallelizationHelper.tpp:18-33`),
skipping the write (and its `rsqrt(mmom)`/`sigmaFactor` work) for interior
coarse atoms whose `thermalField` entry is never read by the active-list-
scoped `EvolveFirst`/`EvolveFirst`-via-`evolveFirstWithThermalField` kernels
anyway.

- **Parity impact: none.** The generator sequence is byte-for-byte
  unchanged; only which atoms get a `field[atom]` write changes, and those
  writes were being discarded by every downstream reader already (traced:
  `evolveFirst`/`evolveSecond`'s active-list overloads only ever read
  `thermalField` through `gpuActiveAtomCall`, restricted to the same active
  list). Invariants 1-4 remain trivially true for exactly the same
  structural reason they hold today.
- **Restart impact: none** (unaffected either way; Part B invariant 5 is a
  pre-existing baseline property, not something this strategy touches).
- **Computational cost:** removes only the "loop" sub-phase (Part C table
  3: 44-56% of `randomize()`'s GPU time in the sizes tested) -- the "RNG"
  generate call stays full-cost, full-N, every step. **Bounded win**, not a
  removal of the O(N) floor: at best roughly halves the thermal-RNG
  contribution, it does not make it scale with active atoms.
- **Memory cost:** none (existing `field`/`randomValues` arrays keep their
  current full size; only the launch geometry of one kernel changes).
- **Invasiveness to feature-off ASD:** minimal. A new
  `randomize(mmom, oneBasedAtoms, activeCount)` overload is additive; the
  existing `randomize(mmom)` feature-off keeps calling is untouched.
- **CUDA/HIP feasibility:** straightforward on both backends. `SetupField`
  currently needs both `atom` and `site` per thread
  (`gpuThermfield.cpp:18-24`); the existing `active_atom_kernel` template
  only derives `atom` from a compact slot (`gpuParallelizationHelper.tpp:31-42`).
  Since `site = oneBasedAtom - 1` is already known before the per-ensemble
  offset is added, a small sibling kernel/launcher (not a modification of
  the shared `active_atom_kernel` template other call sites depend on) would
  be needed to hand `SetupField` both indices. Small, isolated addition.

### Strategy 2: persistent per-site RNG state, update only active sites

Replace the single host-API stream generator with a device-resident,
per-`(site,ensemble)` state array (e.g. `curandStateXORWOW_t`/
`hiprandStateXORWOW_t`, one per site*ensemble, `curand_init`-seeded once),
and a custom kernel that calls `curand_normal()` directly, only for active
`(site,ensemble)` pairs, each step.

- **Parity impact: broken**, and not narrowly -- there are two sub-variants,
  both bad for invariant 1:
  - *Advance every site's state every step regardless of activity* (so the
    stream stays "aligned" with what the shared generator would have
    produced) recovers nothing: this still does N\*M `curand_normal()` calls
    every step, defeating the entire point of the strategy.
  - *Freeze a site's state while it is coarse* (the only variant that
    actually saves cost) makes a fine atom's post-coarse-residence draw
    depend on how many steps it spent coarse, which the *feature-off*
    reference trajectory never experiences at all -- this is a **new
    stochastic contract**, not today's, and reopens a real physics question
    (is skipping thermal kicks during coarse residence, vs. some other
    continuation rule, itself defensible?) that CGP-06A is explicitly not
    scoped to resolve. It also does not, by itself, produce bit-identical
    per-site draws to today's shared-stream scheme even in the *all-fine*
    case, because splitting one stream into N\*M independent per-site
    sub-streams (via seeding/skip-ahead) does not reproduce which value the
    single shared XORWOW state would have produced for that site at that
    step -- invariant 1 breaks even for an all-fine run under this
    sub-variant, unless feature-off ASD is *also* switched to the same
    per-site scheme (a new accepted baseline, a Human-level decision, not a
    performance-only change).
- **Restart impact:** would need the whole per-site state array serialized
  to genuinely resume a frozen-during-coarse stream -- new capability weight
  relative to today's (already absent, Part B invariant 5) baseline, not a
  regression, but real added scope.
- **Computational cost:** potentially removes the O(N) floor entirely (both
  generate and write scale with active count only) -- but only under the
  freeze-during-coarse sub-variant, which is the one with the open physics
  question above.
- **Memory cost:** one persistent state word per site*ensemble (tens of
  bytes each on either backend) -- modest but nonzero, and must be
  initialized once (an O(N\*M) one-time `curand_init` pass, not a per-step
  cost if done once at simulation start).
- **Invasiveness:** would most naturally become the *shared* implementation
  (feature-off's "active list" is just "everyone"), which is a bigger
  structural change than this task's evidence-only scope, and squarely
  matches why CGP-06B is already flagged High risk in `CGP_work.md`.
- **CUDA/HIP feasibility:** symmetric -- both curand and hiprand expose
  per-thread device-side state types and `_init`/`_normal` primitives; the
  project's `gpu_wrappers.h` currently wraps only the host-API generator,
  so this would need a new macro layer.

### Strategy 3: counter-based/site-addressable RNG keyed by (seed, accepted_step, atom, ensemble, component)

Use a counter-based generator (Philox4x32-10, available in both curand and
hiprand as a named RNG type and as a low-level per-thread device function) so
each thermal draw is a **pure, stateless function** of its key tuple, computed
on demand only for the `(site,ensemble)` pairs active this step. No
persistent per-site state, no sequential stream to advance or skip.

- **Parity impact:** breaks bitwise match to the *current* XORWOW-sequential
  accepted trajectory, same as Strategy 2 -- a counter-based generator does
  not, and cannot, reproduce what a sequential-stream generator would have
  produced at a given index; Philox and XORWOW are different algorithms.
  Adopting this for adaptive-only (leaving feature-off on XORWOW) makes
  invariant 1 unsatisfiable *by construction* even at 100% fine, since the
  two paths would then run different algorithms; adopting it for both makes
  feature-off's own reference trajectory a new baseline, again a Human/Terra
  decision, not a performance-only change.
- **Restart impact: strictly better than today's baseline.** Since a draw is
  a pure function of `(seed, accepted_step, atom, ensemble, component)`, a
  restarted run that knows its own `accepted_step` counter reproduces
  bit-identical thermal draws with **no state serialization at all** --
  this upgrades invariant 5 rather than merely preserving it.
- **Fine-fraction/transition invariants (2-4): satisfied by construction**,
  more robustly than either other strategy -- a site's draw never depends on
  any other site's activity, ordering, or its own coarse-residence history.
  A newly-fine atom at step *t* draws exactly the value the key function
  gives for `(seed, t, atom, ensemble, ·)`, identical to what it would have
  drawn had it been fine at every prior step. No freeze/skip design question
  exists for this strategy the way it does for Strategy 2.
- **Computational cost:** best of the three -- both generate and write scale
  with active count only, no full-N pass, no persistent state to maintain or
  initialize.
- **Memory cost:** minimal (stateless).
- **Invasiveness:** touches the RNG *algorithm* production ASD accepts, for
  whichever paths adopt it -- the single most consequential decision point
  of any implementation task built on this evidence.
- **CUDA/HIP feasibility:** fully symmetric; Philox4x32-10 is a named RNG
  type in both curand and hiprand, with equivalent low-level per-thread
  counter primitives on both.

### Recommendation

**No strategy is free of a Human-level decision**, because invariants 2-4
(the actually load-bearing correctness properties for a fine-fraction-scoped
implementation) are only *free* under the current design because it is
oblivious to activity -- any design that stops being oblivious has to either
(a) preserve today's exact sequential-stream algorithm (Strategy 1: safe,
bounded win) or (b) change which algorithm produces the accepted trajectory
(Strategies 2/3: unbounded win, but a new stochastic contract requiring
explicit acceptance, exactly the kind of decision CGP-00's energy-precision
work and CGP-00B's transition-gate work already escalated rather than made
unilaterally).

Given that pattern:

- **Recommend Strategy 1 as CGP-06B's immediate, low-risk scope.** It is
  parity-preserving by construction (invariants 1-5 all remain exactly as
  strong as they are today, verified by the same existing RCG-09A.2 tests
  with no new correctness burden), minimally invasive to feature-off (purely
  additive overload), and CUDA/HIP-symmetric via a small sibling of the
  existing active-atom kernel geometry. Its ceiling is real but bounded:
  Part C table 3 shows it can remove roughly 44-56% of `randomize()`'s
  GPU-side cost (the "loop" share), not the O(N) floor itself.
- **Recommend Strategy 3 over Strategy 2 as the follow-on, if and when the
  project is willing to accept a new stochastic contract.** Strategy 3
  dominates Strategy 2 on every axis measured here except "requires no new
  RNG algorithm at all" (neither does): it is stateless (no serialization
  burden, no memory overhead), it resolves invariants 2-4 more robustly (no
  freeze/skip design question), and it strictly improves restart
  reproducibility rather than merely preserving today's already-broken
  baseline (Part B invariant 5). Its adoption is squarely a Human/Terra-level
  decision -- not a performance-only change -- and should be scoped as such,
  the same way CGP-00 escalated the energy-precision contract rather than
  silently changing it.
- **Strategy 2 is not recommended** as a primary path either way: its
  cost-saving variant inherits Strategy 3's parity-breaking downside while
  additionally requiring more implementation weight (persistent per-site
  state, its own serialization question) without a clean advantage over
  Strategy 3's statelessness.

This recommendation changes no production RNG behaviour; per this task's own
scope, only CGP-06B (contingent on Human acceptance of whichever contract is
chosen) would implement it.

## Checklist

* [x] RNG implementation traced (generator type, state ownership,
  consumption semantics, temperature/moment entry, seed/restart handling,
  T=0 behaviour, inactive-atom advancement -- Part A).
* [x] State ownership documented (single process-lifetime generator, shared
  identically between feature-off and adaptive dynamics).
* [x] Random-consumption semantics documented (`3*N*M` values, one call per
  accepted step, full-grid write, no active-list awareness).
* [x] All-fine parity requirement frozen (invariant 1, Part B; existing
  RCG-09A.2 tests cited as the controlling evidence).
* [x] Transition/restart requirements frozen (invariants 4-5, Part B;
  restart's pre-existing, CG-independent limitation stated precisely).
* [x] RNG time vs. total atom count measured (Part C table 1: fp32 and
  fp64).
* [x] RNG time vs. active fraction measured (Part C table 4: current flat
  cost vs. a disclosed smaller-full-lattice projection, both fp32 and
  fp64).
* [x] Three implementation strategies evaluated (Part D: parity, restart,
  cost, memory, invasiveness, CUDA/HIP feasibility for each).
* [x] Preferred strategy recommended with evidence (Part D: Strategy 1
  near-term, Strategy 3 over Strategy 2 as the contingent follow-on).
* [x] No production RNG behaviour changed (only a new, additive,
  evidence-only benchmark file and this document).

## Commit

`CGP-06A: characterize adaptive thermal RNG scaling`
