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

## Remaining hardware gates

HIP execution, feature-off noise measurements, phase overhead, and the
active-DOF crossover still require suitable hardware runs.  The
implementation does not substitute host wall time or a modelled crossover
for those measurements.
