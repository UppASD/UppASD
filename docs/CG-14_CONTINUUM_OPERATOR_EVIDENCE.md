# CG-14 evidence: reduce the GPU continuum operator constant

Task: `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` §11, "Task CG-14".
Branch: `gpu_hip_cu_ab_cg`. Baseline tip: `abda6534` plus the RCG-10
reconciliation working tree.
Date: 2026-08-15.

**Headline.** The continuum operator kernel is 3.01x faster and the coarse
phase 2.54x faster at the RCG-09C reference size. The 10x target in the task's
acceptance table was **not** reached, and section 4 explains why it was never
reachable by common-subexpression elimination: the kernel was already running
at 68% of this device's FP64 pipeline peak, so the 5.88x of redundant *total*
instructions this task removed could only buy the 2.78x of redundant *FP64*
instructions inside it.

**A crossover now exists, and it is marginal and geometry-dependent.** Against
the production atomistic GPU path at 65 536 atoms, in the all-coarse limit,
with 16 or more atoms per block, the benchmark's own conservative crossover
test passes on 17 of 20 measured runs at ~1.05x. It fails on 0 of 20 baseline
runs. At four atoms per block there is no crossover at any size, in either arm.
It exists only at a fine fraction of exactly zero, and **it does not survive
fp32** (§6.6): in single precision the atomistic reference gains 3.61x against
the adaptive path's 2.04x and the crossover disappears, which shows it is partly
an artifact of this device's 1:64 FP64 ratio. Section 6 states all of this
plainly and section 9 records that **no speedup wording has been restored
anywhere**, and that §14 item 9's amendment is left standing pending Human
decision.

---

## 1. What changed

One kernel, `evaluateAdaptiveCoarseTensor` in
`source/gpu_files/gpuAdaptiveRuntime.cpp`, plus the four device helpers it
calls. No physics, sign, prefactor, tensor-index contract, ownership rule,
launch geometry, or public interface changed. `finalizeAdaptiveCoarseLocal`,
`prolongateAdaptiveGhosts`, the CPU operator and the Fortran layer are
untouched.

| # | Change | Kind |
|---|---|---|
| 1 | `plusBlock` becomes `plusBlockTriple`: the three forward-neighbour indices are computed once per block-thread instead of ~189 times | exact CSE |
| 2 | The four coarse-direction vectors are read once into `delta[3][3]`, and the three physical gradients are computed once each instead of 27 times | exact CSE |
| 3 | The derivative stencil accumulates into a local `contribution[4][3]` and is committed by at most 12 atomics, replacing 324 | **re-association** |
| 4 | `blockVolume[block]`, `coarseBlockMask` for self and the three neighbours, and the constant `inverseBlockTranspose` are hoisted into registers | exact CSE |
| 5 | `crossComponent` evaluates only the component of `m x grad` the spiralization energy reads, instead of all three | exact CSE |
| 6 | `basisCross` evaluates `e_k x v` without the four products that multiply an exact zero | exact for finite `v` |
| 7 | The derivative is summed over `q` (exchange) and over `k` (spiralization) before one stencil call per `p`, instead of one call per term | **re-association** |

Changes 1, 2, 4 and 5 preserve every arithmetic operation and its order.
Change 6 preserves the value for every finite input; it differs only where
`crossDevice` would have formed `0 * inf`, which the runtime's finiteness
validation already excludes. Changes 3 and 7 are re-associations of the field
sum and are treated as such throughout — see section 5.

The spiralization *energy* deliberately keeps its own `(k, p)` loop rather than
being folded into the p-major field loop, so that its accumulation order, and
therefore its value, is unchanged. This costs about 27 FP64 instructions per
block-thread and is what makes the bitwise energy result in section 5.1
possible.

## 2. The grounding claims, re-verified rather than trusted

The task's §4.2 redundancy table was checked by source inspection **and**
confirmed against the hardware counters. Per block-thread at 16 384 blocks:

| Quantity | Claimed | Measured (baseline) |
|---|---|---|
| stencil atomics | 324 (+27 direct = 351) | **351.0** — `smsp__inst_executed_op_global_red.sum` 179 712 / 512 warps |
| `plusBlock` calls | ~189 | 189 by inspection (12 per exchange term x 9, 9 per spiralization term x 9) |
| `physicalGradient` calls | 27 | 27 by inspection |
| `blockVolume` loads | 18 | 18 by inspection |

The atomic count matches the source-inspection prediction exactly, which is
the check that the reading of the kernel was correct.

## 3. Kernel-level result

Nsight Compute, `--kernel-name regex:evaluateAdaptiveCoarseTensor
--launch-skip 200 --launch-count 1`, at 16 384 blocks / 65 536 atoms,
`--texture spiral`, fine fraction 0. Baseline binary built out of tree from a
pristine copy of the pre-change worktree.

| Metric | baseline | CG-14 | ratio |
|---|---|---|---|
| kernel duration | 252 896 ns | **83 936 ns** | **3.01x** |
| kernel constant | 15.44 ns/block | **5.12 ns/block** | 3.01x |
| instructions executed (warp) | 3 547 328 | 602 816 | 5.88x |
| FP64 pipe instructions | 598 016 | 215 040 | 2.78x |
| global reduction (atomic) instructions | 179 712 | 3 072 | 58.5x |
| FP64 pipe, % of peak sustained | 67.86% | **81.95%** | — |
| registers per thread | 80 | 104 | — |
| local memory spilling | 0 | 0 | — |

Per block-thread: 6 928 -> 1 177 instructions, 1 168 -> 420 FP64 instructions,
351 -> 6 atomics. (Six rather than twelve because `commitDerivativeStencil`
skips exactly-zero contributions; section 5.3 shows why that is exact rather
than approximate.)

## 4. Why 10x was not reachable, stated as a measurement

This is the substantive finding of the task, and it contradicts the projection
in the prompt's §4.4 and the parent blueprint's acceptance table.

The redundancy table counted neighbour-index arithmetic, global loads and
atomics. Those are real and this task removed 5.88x of the total instruction
count. But the kernel's binding resource is neither of them: on the RTX A4000
(GA104, compute capability 8.6) FP64 runs at 2 operations per SM per clock,
and Nsight reports the *baseline* kernel already at **67.86% of the FP64
pipeline's sustained peak**. Integer index arithmetic and L2 atomics were
never the limiter; they were hidden behind an FP64 pipe that was two-thirds
saturated.

The consequence is arithmetic:

- FP64 instructions per block-thread fell 1 168 -> 420, a factor 2.78.
- Kernel time fell by 3.01x — the extra 8% is the pipe utilisation rising from
  67.9% to 82.0% as the non-FP64 work that used to interleave with it went
  away.
- At 82.0% of peak, **at most a further 1.22x is available from scheduling,
  occupancy or memory work.** Everything beyond that requires removing FP64
  instructions.

So the operator is now within ~20% of what this device can execute for the
arithmetic it performs. Reaching 10x would have required removing ~90% of the
FP64 operations, and the FP64 operations are the tensor contraction, the
gradients and the derivative stencil — that is the method, not overhead. The
§4.4 projection was sound arithmetic applied to the wrong quantity.

A successor with an appetite for it has one identified route, recorded here
rather than attempted: the exchange scatter
`Sum_p Sum_q 2 V B[p,d] C[p,q] grad_q` can be pre-contracted on the host into a
single constant 3x3 matrix `M[d][d'] = Sum_pq B[p,d] C[p,q] B[q,d']`, reducing
the exchange half of the kernel from ~324 to ~27 FP64 operations. It is not a
common-subexpression elimination — it changes which numbers are multiplied,
and it is only valid when every `(p, q)` term of a block is owned, so it needs
an ownership-gated fast path with a general fallback. That is a separate task
with its own correctness argument, not something to fold into CG-14.

## 5. Numerical effect, identified and bounded

### 5.1 The coarse energies are bitwise identical

The exchange energy keeps its `(p, q)` accumulation order and the
spiralization energy keeps its `(k, p)` order, and the gradients feeding them
are bit-identical (change 2 preserves the operation order of
`physicalGradient` exactly). The prediction is that both coarse energy terms
come out bit-for-bit unchanged on any fixture whose input state is unchanged.

Measured across RCG-04I's 19 `moving_*` fixtures, comparing the baseline CUDA
fp64 binary against the CG-14 CUDA fp64 binary (same CPU reference, same
device, same session):

- **`coarse_exchange`: bitwise identical on 12 of 13** coarse-active fixtures.
- **`coarse_spiralization`: bitwise identical on all 5** fixtures that have it.
- **`coarse_anisotropy`: bitwise identical on all 5** fixtures that have it.

The single exception is `moving_wall_adaptive`'s `coarse_exchange`, at 1.4e-15
relative. That fixture is the only *adaptive* one in the set: its block state
evolves in response to the field, so an ulp-level field difference feeds back
into which blocks are coarse and the input state itself diverges. The
bitwise-energy property is a statement about a fixed input state, and
`moving_wall_adaptive` does not hold its input state fixed.

### 5.2 The field re-association, measured against the run-to-run floor

Changes 3 and 7 re-associate the derivative sum. Bitwise identity is not
available and is not claimed — and could not be, because the pre-change field
was itself accumulated by `atomicAdd` in scheduler order and was never
reproducible run to run (the parent blueprint's trap list says so).

The right comparison is therefore against that pre-existing nondeterminism.
Using the `coarse_norm2` field checksum, the most sensitive published quantity:

| Comparison | worst relative difference |
|---|---|
| baseline vs CG-14, 13 coarse-active fixtures | **2.65e-15** (`moving_dmi_chiral_bs1_minus`) |
| CG-14 vs CG-14, same binary run twice | **1.24e-15** (`moving_dmi_chiral_bs1_minus`) |
| baseline vs CG-14, `coarse_exchange` energy | 1.4e-15 (adaptive fixture only) |
| CG-14 vs CG-14, `coarse_exchange` energy | 3.5e-16 (adaptive fixture only) |

The change-induced difference is within a factor of ~2 to 4 of the noise the
unmodified code already produces between two runs of the same binary. Four of
the thirteen fixtures — every `moving_all_coarse_bs*` — show **exactly zero**
difference.

Against RCG-06B's frozen production budget this is
`2.65e-15 / 1.0e-3 = 2.6e-12` of the allowance, and `3.1e-11` of the observed
CPU/GPU summation-order floor of `8.6e-5` that dominates that budget. Every
headline backend-parity metric — `angular_max_rad`, `component_max`,
`restart_max` — is unchanged to all printed digits on all 19 fixtures, and
`adaptive-cg-moving-backend-parity` passes with its frozen budgets untouched.

### 5.3 The zero-skip in the stencil commit is exact, not approximate

`commitDerivativeStencil` skips an atomic whose contribution is exactly zero.
`clearAdaptiveCoarse` initialises the field to `+0.0` and every later update is
an addition, so the accumulator can never hold `-0.0` (`a + (-a)` is `+0.0`
under round-to-nearest, and `+0.0 + (-0.0)` is `+0.0`). Adding `+/-0.0` to any
reachable accumulator value is therefore the identity, and skipping it is
exactly equivalent.

## 6. Performance: measurement and result

### 6.1 Protocol

Both arms are fresh out-of-tree builds. The baseline arm is built from an
rsync'd pristine copy of the pre-change worktree at
`/tmp/cg14_baseline`, so no build directory is shared and the stale-`.mod`
trap cannot apply. Per RCG-09C section 2 and RCG-08-FU3: at each size the first
invocation is a discarded warm-up, then the arms alternate **ABBA**, four runs
per arm (eight for the confirmation run), medians reported. `nvidia-smi` was
sampled before and after every size.

Device state, recorded rather than assumed. Across every measurement in this
section: **no compute processes other than the benchmark** (`nvidia-smi
--query-compute-apps` empty throughout); utilisation 0–47%, all of it the
benchmark; SM clocks 1725–1935 MHz while active;
`clocks_throttle_reasons.active` `0x0` while active (`0x1` = GpuIdle when
idle); temperature 50–60 °C, no thermal throttle observed. GPU 0 was allowed
to cool below 52 °C before the fine-fraction sweep and before the confirmation
run. This host's two RTX A4000s are shared with another user; they were idle
for the duration, which is recorded above rather than assumed, and the ABBA
interleaving is designed so that any drift hits both arms equally. The
production-oracle row is identical code in both binaries and agrees between
arms to within 1.6% at every size, which is the check that both arms saw the
same device.

Command, both arms:

```sh
gpu_adaptive_runtime_benchmark --blocks B --atoms-per-block A \
    --warmup 1 --iterations 3 --repetitions 3 --texture spiral
```

Raw JSON: `docs/cg14/cg14_abba_fine.json`, `docs/cg14/cg14_abba_fp32.json`,
`docs/cg14/cg14_abba_blocksize.json`, `docs/cg14/cg14_abba_confirm.json`.

### 6.2 Fine-fraction sweep, four atoms per block

All times microseconds per adaptive step; `coarse` is the device-event phase
time, `step wall` the untimed step. Medians of four runs per arm.

**blocks = 16 384, atoms = 65 536.** Production oracle: baseline 284.0,
live 279.4.

| fine frac | coarse blocks | coarse base | coarse live | step base | step live |
|---|---|---|---|---|---|
| 1.0000 | 0 | 16.6 | 16.4 | 324.9 | 316.6 |
| 0.7500 | 4 094 | 283.6 | 119.5 | 573.7 | 422.0 |
| 0.5000 | 8 190 | 367.8 | 122.4 | 616.2 | 395.1 |
| 0.2500 | 12 286 | 372.6 | 125.9 | 591.9 | 364.1 |
| 0.1250 | 14 334 | 528.4 | 203.6 | 720.1 | 421.5 |
| 0.0625 | 15 358 | 529.4 | 203.8 | 713.4 | 418.3 |
| 0.0000 | 16 384 | **523.9** | **206.4** | **656.3** | **365.0** |

Coarse-phase constant at all-coarse: **31.98 -> 12.60 ns per block (2.54x)**.
Whole step at all-coarse: 1.80x. Crossover: **not observed**, 0/4 runs in
either arm. Best adaptive step across the sweep is all-fine at 316.6 µs
against 279.4 µs production, i.e. 1.13x slower — the all-fine point does not
launch this kernel and is unchanged, as expected.

**blocks = 4 096, atoms = 16 384.** Production oracle: baseline 98.2, live 97.7.
All-coarse coarse phase 286.1 -> 121.9 (2.35x), constant 69.85 -> 29.77
ns/block, step 343.7 -> 196.6 (1.75x). Crossover not observed.

**blocks = 1 024, atoms = 4 096.** Production oracle: baseline 50.8, live 50.8.
All-coarse coarse phase 291.3 -> 123.8 (2.35x), constant 284.50 -> 120.93
ns/block, step 334.5 -> 182.3 (1.83x). Crossover not observed.

The coarse phase does not fall below ~120 µs at any size, including 1 024
blocks where the tensor kernel itself is a few microseconds. That floor is the
phase's other launches plus the event synchronisation `finishPhase` performs at
every phase boundary; it is the same instrument RCG-09C used and is reported
here for comparability, not as a property of the operator. Whole-step figures
use the untimed `step_wall_us`.

### 6.3 Block-size sweep at fixed 65 536 atoms

Four atoms per block is the least favourable geometry for the method, so the
sweep holds the atom count fixed and trades blocks for block size. Medians of
four runs per arm.

| blocks | atoms/block | coarse base | coarse live | ns/block base -> live | all-coarse step base | step live | production | crossover (live) |
|---|---|---|---|---|---|---|---|---|
| 16 384 | 4 | 523.9 | 206.4 | 31.98 -> 12.60 | 656.3 | 365.0 | 279.4 | 0/4 |
| 8 192 | 8 | 345.4 | 119.6 | 42.16 -> 14.60 | 477.4 | 276.7 | 274.4 | 0/4 |
| 4 096 | 16 | 262.9 | 114.9 | 64.19 -> 28.06 | 399.0 | **261.8** | 275.2 | **4/4** |
| 2 048 | 32 | 261.6 | 113.0 | 127.75 -> 55.19 | 392.8 | **261.1** | 274.5 | **3/4** |
| 1 024 | 64 | 261.5 | 112.1 | 255.37 -> 109.50 | 391.7 | **261.7** | 275.7 | **3/4** |

Confirmation run at 4 096 blocks / 16 atoms per block, eight runs per arm on a
cooled device: **live 7/8 PASS** (speedups 1.031, 1.046, 1.046, 1.048, 1.048,
1.050, 1.052), **baseline 0/8**. The one live miss is the first run of the
series, still warming — its production oracle read 323.5 µs against 274–275 µs
once settled, and its zero-fine step 309.4 µs against 261–262 µs.

Aggregating every run at 16 or more atoms per block: **live 17/20 PASS,
baseline 0/20.**

### 6.4 Does a crossover exist? Plainly

**Yes, marginally, and only in a specific corner.** The benchmark's own
conservative test — the adaptive step plus `3 x (MAD_production +
MAD_adaptive)` must fall below production times `(1 - 0.02)` — passes for the
CG-14 arm at the all-coarse point when there are 16 or more atoms per block, at
a measured **1.03x to 1.06x**. It does not pass for the baseline arm anywhere.

Everything that qualifies this:

- It exists **only at fine fraction 0.0**, the all-coarse limit. At every
  mixed resolution the adaptive step is still 1.3x to 1.5x slower than
  production.
- It **depends on block size**, decisively. At four atoms per block there is no
  crossover at 1 024, 4 096 or 16 384 blocks; at eight atoms per block there is
  none either. The answer to the task's question is that the crossover does
  depend on the four-atoms-per-block geometry, and that geometry does not have
  one.
- It is **~5%**, close enough to the 2% acceptance margin that the test flips
  between PASS and NOT_OBSERVED across repeats at two of the three qualifying
  sizes, and is only unanimous at 4 096 blocks / 16 atoms per block.
- It compares an all-coarse adaptive step against a full-resolution atomistic
  step. That is what coarse graining is for, and blueprint item 6 validates
  all-coarse long-wavelength accuracy, but it is a comparison between different
  resolutions and should not be read as a like-for-like speedup.

For context against the pre-CG-14 state: at 16 384 blocks the all-coarse step
was 2.02x the all-fine step before this change (656.3 / 324.9) and is 1.15x
after (365.0 / 316.6). The all-coarse penalty for coarsening, which is what
made the method a pessimisation per degree of freedom, is largely gone.

### 6.5 The coarse phase still dominates, but barely

At 16 384 blocks, all-coarse, after CG-14:

| Phase | before | after |
|---|---|---|
| coarse | 523.9 | **206.4** |
| interface (`prolongateAdaptiveGhosts` + `clearAdaptiveInterface`) | 124.3 | 123.4 |
| integration | 38.9 | 38.2 |
| atomistic | 25.9 | 25.9 |

The coarse phase is still the largest term, but its margin over the interface
phase fell from 4.2x to 1.67x. Two successors are named:

1. **`prolongateAdaptiveGhosts`, 123.4 µs**, unchanged by this task and the
   expected successor the prompt anticipated. RCG-09C section 6 records why it
   cannot simply be restricted to the ghost shell.
2. **The coarse phase's own non-tensor remainder.** The tensor kernel is now
   83.9 µs of a 206.4 µs phase, so more than half of the coarse phase is the
   four small coarse kernels plus launch and synchronisation latency. At
   1 024 blocks the phase does not go below ~120 µs at all. If the coarse phase
   is attacked again, that fixed remainder — not the operator — is where the
   time is.

### 6.6 The crossover does not survive fp32, and probably not a datacenter GPU

Added 2026-08-15 after the sections above, in answer to "can the FP64 parts go
away?". It qualifies §6.4 materially and is recorded here rather than left for
a reader to discover.

The same ABBA protocol, both arms rebuilt with `UPPASD_PRECISION=SINGLE`
(the production oracle is fp32 in both arms, so the comparison is internally
fair):

| | fp64 | fp32 | fp32 gain |
|---|---|---|---|
| coarse kernel (16 384 blocks) | 83.9 µs | **16.0 µs** | 5.25x |
| coarse-phase constant | 12.60 ns/block | **5.18 ns/block** | 2.43x |
| adaptive all-coarse step, 4 096 blocks / 16 atoms per block | 261.8 µs | **128.4 µs** | 2.04x |
| **production atomistic step, same geometry** | 275.2 µs | **76.2 µs** | **3.61x** |
| crossover | **PASS (1.05x)** | **NOT_OBSERVED (0.59x)** | — |

The kernel gets 5.25x faster in fp32 and the crossover **disappears**, because
the production atomistic path gains more from single precision (3.61x) than the
adaptive path does (2.04x). The atomistic reference does FP64 work proportional
to 65 536 atoms and is punished hard by this device's 1:64 FP64 ratio; the
adaptive all-coarse path does far less FP64 work but carries a fixed residual of
launch latency and atom-proportional passes (§6.5) that no precision change
touches.

Two consequences, stated plainly because they cut against §6.4:

1. **The fp64 crossover is partly an artifact of the measuring device.** It
   exists because an RTX A4000 penalises the atomistic reference more than it
   penalises the coarse operator. On a datacenter GPU with a 1:2 FP64 ratio
   (A100, H100) the atomistic reference would speed up by roughly the factor
   fp32 gives here while the adaptive path's fixed residual would not, so the
   crossover should be expected to narrow or vanish. **It has not been measured
   on such a device and must not be assumed to transfer.**
2. **Removing FP64 is not the lever it looks like.** It makes everything
   faster and the adaptive path relatively worse. FP64 is worth attacking for
   absolute throughput, not for the crossover.

One targeted FP64 reduction does survive this argument. In the fp32 build the
kernel still issues 26 624 FP64 instructions — 52 per block-thread, the
RCG-06B-mandated `double` energy accumulation — and they still occupy **48.7%
of the FP64 pipe**, against 36.3% overall SM throughput. That accumulator is
the largest single utilisation number in an otherwise-fp32 kernel. Accumulating
a block-thread's nine energy terms in `real` and widening to `double` only at
the block reduction would keep the global sum in FP64, where RCG-06B's
microbenchmark actually demonstrated the need (N up to 3e6), while removing the
per-term FP64 traffic. Estimated up to ~2x on the fp32 kernel. It is a
summation-precision change and would need its own bound; it is not attempted
here.

### 6.7 The ceiling on any further continuum-operator work

At the qualifying geometry (4 096 blocks, 16 atoms per block, all-coarse), the
step is 261.8 µs of which the coarse phase is 114.9 µs. **A continuum operator
that cost nothing at all** would give 146.9 µs against 275.2 µs production, i.e.
**1.87x — that is the hard ceiling on every remaining optimisation in this
kernel and its phase.** In fp32 the same arithmetic gives 70.2 µs against
76.2 µs, a ceiling of 1.09x.

The remaining 146.9 µs is interface (116.6), integration (24.2) and atomistic
(26.3) — all proportional to atoms, none of them falling with block size
(§6.5, and the flat 261.8 / 261.1 / 261.7 µs across 16, 32 and 64 atoms per
block in §6.3). Further work on the continuum operator is therefore
subject to sharply diminishing returns, and `prolongateAdaptiveGhosts` is not
merely the next target but the one that sets the ceiling for everything else.

## 7. Correctness

All results on **fresh out-of-tree builds** (`build_cg14_cuda`,
`build_cg14_cpu`, `build_cg14_cuda_fp32`, each `rm -rf`'d and reconfigured
after the last source edit).

| Gate | Result |
|---|---|
| `ctest -L cg13-cuda` (fp64) | **32/32 passed** |
| `ctest -L cg13-cpu` | **29/29 passed** |
| `coarse-graining-dispersion` (RCG-04-FU4), CPU and CUDA trees | passed |
| `gpu_adaptive_runtime_benchmark --parity-only --require-acceptance` | **exit 0, 205 PASS, 0 FAIL** |
| `adaptive-cg-moving-backend-parity` | passed, frozen budgets untouched |
| `gpu_adaptive_runtime_tests` incl. the new fixture | passed |
| `ctest -L cg13-cuda` (fp32) | 31/32; see section 7.3 |

### 7.1 The adjoint pair

`physical_forward_gradient` / `add_physical_gradient_transpose` are Fortran
routines in the CPU operator and were not touched. Their exactness is held by
`coarse-graining-tensor-operator` and `coarse-graining-dispersion`, both of
which pass unchanged in both trees. On the GPU side the same adjoint structure
is what the new host-reference fixture (section 7.4) checks directly: the
fixture compares the device field against a host derivative built by the
literal transpose scatter, so a broken adjoint shows up as a field mismatch.

### 7.2 The RCG-09A.4 parity gate

All 205 stages PASS with `--require-acceptance`, including bitwise
thermal-field identity and bitwise `depondt-predictor` /
`depondt-corrector-final-spin` on the T=0 cases.

The residual `1e-16`-level differences between the baseline and CG-14 parity
logs are **not** attributable to this change. Running the *same* binary twice
produces 25–34 differing lines; running baseline against CG-14 produces 34. The
all-fine parity path launches no coarse blocks and therefore never enters this
kernel; the wobble is the pre-existing `atomicAdd` nondeterminism of the
atomistic bond kernel.

### 7.3 fp32, reported separately, inheriting no fp64 claim

`ctest -L cg13-cuda` on a fresh fp32 build: 31/32, with
`adaptive-cg-production-e2e` failing at
`assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0`.

This is **pre-existing**: the identical assertion fails identically on a fresh
fp32 build of the *pristine baseline worktree* (`/tmp/cg14_baseline`,
31/32, same test, same assertion). It is RCG-04-FU5, already on record. Nothing
in the fp32 arm is introduced by CG-14, and no fp64 result in this document is
claimed for fp32.

### 7.4 Negative controls, and an oracle gap this task had to close

Five defects were injected into the rewritten kernel one at a time, rebuilt,
and run against the oracles cheapest-first. The first result was that **three
of the five were invisible to every existing fixture**:

| Defect | dispersion | existing unit tests | production-e2e | backend-parity |
|---|---|---|---|---|
| NC-1a: stencil neighbour along direction 2 dropped | pass | pass | pass | pass |
| NC-1b: stencil neighbour along direction 0 dropped | pass | pass | **FAIL** | **FAIL** (36 fixtures) |
| NC-2: `inverseBlockTranspose` index transposed | pass | pass | pass | pass |
| NC-3: basis-cross chirality inverted | pass | pass | **FAIL** | **FAIL** (15 fixtures) |
| NC-4: exchange derivative index `p`/`q` swapped | pass | pass | pass | pass |

The three misses are structural, not luck:

- Every geometry that runs end to end today is **orthogonal** — RCG-05-FU3:
  `neighbourmap.f90` rejects the skew fixture at setup — so
  `inverseBlockTranspose` is diagonal and transposing it is the identity
  (NC-2).
- The tracked materials are **cubic**, so `exchangeStiffness` is diagonal and
  `Sum_q C_pq grad_q` collapses to a multiple of `grad_p` (NC-4).
- Every tracked texture varies along **one axis**, and with a diagonal metric
  that makes the transverse stencil branches carry exact zero (NC-1a). The
  fixtures' block grid is 24x2x2, so direction 2 is a genuine distinct
  neighbour; the branch is live, it just carries nothing.

Also worth recording: `coarse-graining-dispersion` did not catch any of the
five. It exercises the Fortran `CoarseTensorOperator` on the CPU
(`tests/coarse_graining/test_coarse_dispersion.f90` uses
`evaluate_coarse_tensor_operator`, no GPU path), so it is not an oracle for a
GPU-only kernel edit. It is a sharp oracle for the operator's algebra, and it
passes, but the task prompt's expectation that it would be the most
discriminating gate for this change is not borne out.

**The response was to add the missing oracle, not to report the gap and move
on.** `testContinuumOperatorAgainstHostReference` in
`tests/coarse_graining/test_gpu_adaptive_runtime.cpp` builds a 3x3x3 periodic
all-coarse fixture with a dense non-symmetric `inverseBlockTranspose`, a dense
symmetric `exchangeStiffness` (`validate()` requires symmetry;
symmetric-but-not-diagonal is what discriminates the `p`/`q` contract), a dense
`spiralization`, and a direction field varying along all three block
directions. It compares the device energies and published coarse field against
a host reference written from the operator's definition in the *pre-CG-14*
shape — per-term gradients, per-term scatter — so agreement is evidence about
the rewrite rather than a restatement of it. The fixture asserts its own
non-triviality before comparing.

Re-running all five controls against it:

| Defect | first oracle to fail | reported |
|---|---|---|
| NC-1a | **new host-reference fixture** | field deviation 6.12 |
| NC-1b | **new host-reference fixture** | field deviation 6.14 |
| NC-2 | **new host-reference fixture** | coarse exchange energy mismatch |
| NC-3 | **new host-reference fixture** | field deviation 1.22 |
| NC-4 | **new host-reference fixture** | field deviation 0.81 |

**5 of 5 caught, and the new fixture is the first to fail in every case** — in
under a second, against the 3–8 seconds the whole-run fixtures need for the two
they can see. All five injections were reverted; the acceptance builds in
section 7 were made after restoring the source from a pre-injection copy and
deleting and reconfiguring every build tree.

## 8. The two options the task asked to be decided explicitly

### 8.1 Stencil scatter to gather over minus-neighbours: **not adopted**

The case for it is real: it removes the atomics entirely and would make the
coarse field deterministic by construction rather than by scheduling. The case
against it is now measured rather than argued.

A gather has each block-thread evaluate not only its own terms but those of its
three minus-neighbours, since those are the threads that contribute the `+`
half of the stencil to it. That means computing four sets of three gradients
instead of one, and running the term loop four times: the FP64 instruction
count, which section 4 shows is the *only* thing that sets this kernel's
runtime, would rise by roughly 4x. It also widens the read footprint from 4
distinct blocks to about 13 (self, its three plus-neighbours, the three
minus-neighbours, and each minus-neighbour's other two plus-neighbours).

The resource it would spend is the binding one; the resource it would save is
one that CG-14 has already reduced from 351 atomics per block-thread to 6, and
that Nsight shows contributes nothing to the FP64 pipe occupancy. A ~4x FP64
increase against a 3.01x-faster kernel would land the operator slower than the
pre-CG-14 baseline. It is not adopted.

Determinism by construction therefore **remains unclaimed** for the coarse
field, as it was before this task. If it is wanted, the cheap route is not a
gather: it is a fixed-order reduction over the four threads that contribute to
each address, which costs nothing in FP64 and is a scheduling change rather
than an algorithmic one. That is a separate task.

### 8.2 Orthogonal-cell fast path: **not adopted here, deferred with an estimate**

When `inverseBlockTranspose` is diagonal — which, per RCG-05-FU3, is every
geometry that currently runs end to end — `B[p + 3d]` vanishes unless `p == d`.
Then each physical gradient has one nonzero direction instead of three (27 FMA
-> 3 multiplies) and the stencil for physical `p` touches only neighbour `p`
(162 -> 54 FP64 operations). Estimated saving ~130 of the 420 FP64 instructions
per block-thread, i.e. **~30% of the kernel, ~25 µs of a 365 µs all-coarse
step, ~7%**.

Not adopted, for three reasons stated rather than implied:

1. It **cannot change the verdict in section 6.4.** The crossover already
   exists at 16+ atoms per block and already does not exist at four; a 7% step
   improvement moves neither boundary.
2. It doubles the number of code paths through an operator whose correctness
   rests on index contracts that, until this task added
   `testContinuumOperatorAgainstHostReference`, no fixture could discriminate
   at all. Two paths would need two sets of controls.
3. Section 6.5 shows the remaining coarse-phase time is more than half
   non-tensor. Optimising the tensor kernel further is now the smaller half of
   the smaller problem.

It is deferred rather than rejected. The new host-reference fixture uses a
dense non-symmetric `inverseBlockTranspose`, so it is exactly the guard that
would make a two-path implementation safe, and a successor should run the same
five negative controls against both paths.

## 9. Claim discipline

**No speedup wording has been restored in any document.** The parent
blueprint's §14 item 9 amendment and §14.2, and `docs/CG-13_RELEASE_VALIDATION.md`,
are untouched by this task.

Section 6.4 demonstrates a crossover for the first time — marginal, all-coarse
only, and at block sizes of 16 atoms and above. Per the task's instruction,
this document **proposes that §14 item 9's amendment be revisited** and does
not revisit it. That decision requires Human approval, and the following should
be weighed with it:

- the crossover is ~5% at the benchmark's own 2% acceptance margin, and is not
  unanimous across repeats at two of three qualifying sizes;
- it does not exist at four or eight atoms per block;
- it exists **only at a fine fraction of exactly zero**. At 6.25% fine atoms
  the adaptive step is already 0.86x production, i.e. slower. An adaptive run
  with any fine region at all does not benefit (§6.3);
- it compares an all-coarse step against a full-resolution atomistic one;
- **it does not survive fp32, and is partly an artifact of this device's 1:64
  FP64 ratio** (§6.6). It has not been measured on a GPU with datacenter FP64
  and should not be assumed to transfer.

On that evidence a defensible reading is that the amendment should stand as
written, with §6.4 recorded as a bounded, geometry- and device-specific
observation rather than as a speedup claim. Until Human approval is recorded
the amendment stands and the branch remains a correctness/reference
implementation.

## 10. Environment

```
backend CUDA; precision fp64 (fp32 reported separately in 7.3)
CUDA compiler /usr/local/cuda-13.3/bin/nvcc, release 13.3 V13.3.73
Release flags -O3 -DNDEBUG; CMAKE_CUDA_ARCHITECTURES=native (sm_86)
Fortran GNU Fortran (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
Linux 6.8.0-137-generic x86_64
CPU 11th Gen Intel(R) Core(TM) i9-11900 @ 2.50GHz
GPU 0: NVIDIA RTX A4000 (GA104, compute capability 8.6), 16376 MiB,
       driver 610.57.04, FP64 at 2 ops/SM/clock
Nsight Compute and Nsight Systems 2026.1.3 from the same toolkit
```

**HIP: not measured.** No HIP toolchain and no AMD device exists in this
environment, as in every RCG-0x session. Every launch this task touched goes
through the existing `ADAPTIVE_LAUNCH_SHARED` macro and no launch site changed,
so the CUDA and HIP spellings remain structurally identical by construction —
but that is a compile-time property and is not evidence of HIP behaviour.
RCG-09-FU4 stands. **Unavailable is not passing.**

## 11. Limitations

1. The 10x acceptance target was not met; 3.01x on the kernel and 2.54x on the
   coarse phase were. Section 4 gives the measured reason and the remaining
   headroom (~1.22x from scheduling, the rest only from fewer FP64
   instructions).
2. The crossover in section 6.4 is marginal, geometry-dependent, confined to a
   fine fraction of exactly zero, and does not survive fp32 (section 6.6). It
   is reported with every qualification rather than as a result.
2b. Only one device was measured. The FP64-roofline finding in section 4 and
   the crossover in section 6.4 are both properties of an RTX A4000's 1:64
   FP64 ratio. Neither has been measured on a datacenter GPU, where the
   atomistic reference would gain far more than the adaptive path.
3. The field re-association is bounded empirically (section 5.2), not proven.
   The bound is tight — within 2–4x of the code's own run-to-run floor — but it
   is a measurement over 19 fixtures, not a theorem.
4. `basisCross` differs from `crossDevice` on non-finite input (`0 * inf`).
   The runtime's `validate()` rejects non-finite topology and kernel inputs, so
   this is unreachable, but it is a behavioural difference and is recorded as
   one rather than described as an identity.
5. No skew-cell geometry runs end to end (RCG-05-FU3), so the non-orthogonal
   branch of this kernel is exercised only by the new unit fixture, never by a
   whole-run test.
6. HIP is unmeasured (section 10).
7. The coarse-phase floor of ~120 µs (section 6.2) was characterised but not
   attacked; it is now a larger share of the coarse phase than the operator is.
