# CG-10 GPU adaptive kernel contract

## Accepted scope

`GpuAdaptiveRuntime` now owns the backend-neutral CUDA/HIP implementation of
the accepted single-ferromagnetic-channel adaptive short-range workflow.  Its
optional `KernelInput` is inactive when `atomMoment == nullptr`; in that mode
the CG-09 device tensor inventory and execution path are unchanged.

When enabled, initialization validates and stages:

- atom moments and the eight-entry normalized-prolongation stencil;
- unique atomistic bond pairs and Fortran-order `(3,3,bond)` tensors;
- selector neighbour pairs;
- the physical inverse block transpose, exchange stiffness, and
  spiralization tensors in Fortran `(row,column)` order; and
- zero-, one-, or two-axis UppASD-compatible anisotropy.

All immutable arrays and runtime scratch are included in
`GpuAdaptiveRuntime::estimateBytes()`, allocated through `GpuTensor`, and
released through the existing tracked cleanup.  No kernel performs a runtime
allocation.

The first implementation intentionally uses stable serial reference kernels
inside each phase.  Compact active-atom and active-block lists control field
writeback and integration, while the accepted interface adjoint retains its
complete immutable projection stencil.  This establishes fp64 correctness
before scan/reduction and launch-shape tuning.

After fp64 acceptance, the hot prolongation and interface work was split
across atom/ensemble work items, and coarse tensor, local, dipole, and
writeback work was split across the compact active-block list.  Coarse
derivative and energy collisions use the backend's existing real-valued
atomic additions.  Ordered launches separate zeroing, derivative assembly,
local field conversion, all-grid dipole addition, and writeback so no
cross-grid synchronization is assumed.

Selector scoring now clears scores in parallel and evaluates
edge/ensemble work items with an exact nonnegative atomic maximum.  State
proposal and target-centric buffer dilation use separate ordered launches,
preserving the serial proposal/dilation semantics without write races.
Compaction uses a stable three-channel inclusive scan for active atoms,
active blocks, and interface atoms followed by ordered one-based scatter.
Both scan buffers are tracked allocations included in memory preflight.

## Physics and state contracts

The device implementation retains the CPU definitions:

```text
E_bond = -s_i^T K_ij s_j
B_i    = -(mu_B mu_i)^-1 dE/ds_i

E_A = V sum_pq A_pq (d_p m).(d_q m)
E_D = V sum_kp D_kp e_k.(m x d_p m)

B_coarse = -(mu_B mu_block)^-1 dE/dm
```

Physical forward gradients and their transpose use the same
`inverse_block_transpose(physical,lattice)` indexing as
`CoarseTensorOperator`.  A gradient density is included only when its full
forward stencil is coarse.  Atomistic bonds are included exactly once when
at least one endpoint is atomistic.  Ghost reactions use the exact
moment-weighted adjoint of normalized trilinear prolongation.

Selector scores are `max(0, 1-s_i.s_j)` and contribute symmetrically to both
endpoint blocks.  Proposals do not mutate published state.  Publication:

- is rejected unless called at a complete integration-step boundary;
- accepts an optional per-block energy-gate mask;
- leaves rejected directions, state, and transition epoch unchanged;
- uses `(global seed, block, channel, ensemble, prospective epoch)` for cone
  reconstruction; and
- rebuilds compact work only after publication.

Aligned reconstruction rejects a reduced or zero resultant.  Cone
reconstruction rejects zero, oversized, or unrepresentable resultants and
rolls back the block direction on failure.

The FFT dipole field is an optional already-dispatched all-grid input.  Its
energy is never masked by adaptive short-range state, and this work introduces
no adaptive dipole grid or resolution.

## Phase accounting

`GpuAdaptivePhaseMetrics` reports independent accumulated device-event times
for:

- atomistic bonds;
- coarse tensor terms;
- prolongation/interface restriction;
- selector scoring/proposals;
- compaction;
- the externally dispatched FFT phase; and
- restriction, reconstruction, and Heun integration.

`recordFftMilliseconds()` is the explicit seam from the existing FFT dipole
dispatch.  It records timing only and cannot change adaptive state.

## Parity fixtures

`coarse-graining-gpu-adaptive-runtime` now contains a mixed
coarse/buffer/fine fixture with two atoms per block.  It checks:

- exact preflight/tracker agreement and cleanup;
- moment-weighted restriction;
- atomistic unique-pair sign and field;
- coarse anisotropy energy and tensor-index convention;
- an analytic external-field energy derivative;
- independent all-grid dipole energy ownership;
- non-tie selector decisions and complete-step-only publication;
- accepted transition plus rejected energy-gate rollback;
- exact aligned reconstruction;
- repeatable constrained-cone reconstruction and requested resultant;
- compact list rebuilding; and
- normalized deterministic Heun updates without an in-stage state change.

The same source compiles under the CUDA/HIP backend selected by CMake.

## Precision budgets

FP64 is the acceptance precision.  For the dimensionless `O(1)` parity
fixture:

| Quantity | FP64 budget |
|---|---:|
| Per-term and total energy | `5e-12` absolute |
| Field components | `1e-11` absolute |
| Directional energy derivative | `2e-12` absolute |
| Restriction/resultant | `2e-10` absolute |
| Unit-direction norm squared | `2e-12` absolute |
| Selector score | `2e-12` absolute |

For physical-SI fixtures, use the larger of that fixture's scaled absolute
budget and `2e-11` relative error for nonzero values.  State parity fixtures
must keep every score at least `1e-10` away from both thresholds; exact
threshold ties are policy-dependent and are not parity fixtures.

FP32 is a separate, post-fp64 budget:

| Quantity | FP32 budget |
|---|---:|
| Per-term and total energy | `2e-5` absolute, `5e-5` relative |
| Field components | `2e-5` absolute, `5e-5` relative |
| Restriction/resultant | `2e-4` absolute |
| Unit-direction norm squared | `2e-5` absolute |
| Selector score | `2e-5` absolute |

FP32 state fixtures must remain at least `2e-4` away from either threshold.
HIP fp32 remains an independent execution-acceptance gate even though the
common source compiles through the HIP path.

## CUDA execution acceptance

The CUDA hardware run accepts:

- the fp64 `gpu_adaptive_runtime_tests` parity fixture;
- the fp32 `gpu_adaptive_runtime_tests` budget fixture;
- the existing fp64 `dipole_gpu_fft_tests` suite; and
- Compute Sanitizer memcheck of the fp64 adaptive fixture with zero errors.

The FFT run reports a maximum periodic two-basis field error of
`1.7764e-14`, maximum energy error of `1.1369e-13`, and direct-oracle field
and energy maxima of `3.8858e-16` and `1.8041e-16`, respectively.  The
projected-block field and energy errors are both `4.8573e-17`.

## Hardware gates and performance protocol

CUDA feature-off noise, phase overhead, and active-DOF crossover measurements
are recorded below.  HIP execution remains a separate hardware gate.  The
implementation does not substitute host wall time or a modelled crossover
for hardware measurements.

### Performance acceptance harness

Build and run the non-CTest hardware benchmark:

```bash
cmake --build build_gpu --target gpu_adaptive_runtime_benchmark -j
./build_gpu/bin/gpu_adaptive_runtime_benchmark \
  --blocks 2048 --atoms-per-block 4 \
  --warmup 2 --iterations 10 --repetitions 7 \
  --require-acceptance
```

Increase `--blocks`, `--iterations`, and `--repetitions` if the reported
median absolute deviations are comparable with the requested acceptance
margin.  Keep the GPU otherwise idle and record its model, driver, clock
policy, precision, and full benchmark output.

The three remaining checklist items have distinct evidence requirements:

1. **Feature-off performance:** the benchmark alternates a baseline
   atomistic kernel batch with the identical batch while a default,
   uninitialized `GpuAdaptiveRuntime` is in scope.  Its inventory delta must
   be zero.  The absolute median difference must be no larger than the greater
   of 3% of baseline or three times the combined median absolute deviations.
   The output line must say `feature-off ... result=PASS`.
2. **Selector/compaction overhead:** record the complete
   `adaptive-overhead` line.  It reports selector device and wall time,
   compaction device time, localized host wait, complete wall time, bytes per
   block-mask update, and overhead relative to a mixed-resolution field step.
   This gate requires a numeric report, not a pass/fail threshold.
3. **Active-DOF crossover:** the sweep compares each real compact mask with
   the all-atomistic median.  A crossover is accepted only when the candidate
   median plus three combined MADs is at least 2% below the atomistic median.
   The output must say `active-dof-crossover result=PASS` and records the
   active-DOF ratio, requested fine fraction, time, and speedup.

`--require-acceptance` returns exit status 2 if either the paired feature-off
gate or crossover is not observed.  A `NOT_OBSERVED` crossover is a valid
measurement but does not tick the blueprint box; it indicates that kernel
optimization or a larger problem is still needed.  The benchmark is compiled
from the same source for CUDA and HIP and supports the same command-line
protocol in fp64 and fp32 builds.

### CUDA performance measurement

The fp64 acceptance command above was run on an NVIDIA RTX A4000 with driver
610.43.02 using 2048 blocks, four atoms per block, two warmups, ten measured
iterations, and seven repetitions.  Before tuning, the robust sweep returned
`NOT_OBSERVED`: its all-atomistic median was 58568.20 us and medians increased
to 166315.81 us at a 0.250 active-DOF ratio.  This negative result triggered
the compact parallel coarse/interface work described above.

On the post-optimization acceptance run, paired feature-off medians were
38.349 us baseline and 38.357 us with the inactive runtime present, a +0.023%
delta with zero inventory change.  This passes the 3%/three-MAD feature-off
budget.

Before control-phase tuning, at the 50% requested fine fraction selector wall
time was 5990.19 us and compaction wall time was 725.77 us per update,
together 31.00% of the mixed field step.  After parallel selector and stable
scan compaction, selector wall time is 25.05 us and compaction wall time is
41.58 us: approximately 239x and 17.5x faster, respectively.  Their combined
overhead is 0.308% of the 21668.22 us mixed field step.  Selector device time
is 15.96 us; compaction device time is 35.56 us, including 14.56 us of
localized host wait.  Compaction still transfers exactly 8192 mask bytes per
update.

The post-optimization crossover is accepted at a 0.813232 active-DOF ratio.
The all-atomistic median was 40705.61 us with a 2.81 us MAD; the crossover
median was 31274.23 us with a 2.93 us MAD.  This is a 1.3016x speedup and is
well separated from the required 2% plus three-combined-MAD margin.  The
zero-fine median was 2368.32 us.  The optimized fp64 and fp32 parity fixtures
pass, Compute Sanitizer reports zero errors, and the fp64 FFT dipole suite
continues to pass.
