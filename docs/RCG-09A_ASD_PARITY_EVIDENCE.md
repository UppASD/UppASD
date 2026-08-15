# RCG-09A.4 adaptive all-fine ASD parity

Status: harness implemented; local CPU/build validation is complete. CUDA
runtime evidence is recorded separately for the available-device handoff
host; the current validation container has no usable CUDA device. HIP remains
unavailable and is not inferred to pass.

## Contract

Fine atomistic adaptive dynamics is production UppASD Depondt ASD.

The feature-off arm is the real production GPU path: `GpuHamiltonianCalculations::heisge`
followed by `GpuDepondtIntegrator::evolveFirst/evolveSecond`. The all-fine
adaptive arm uses `GpuAdaptiveRuntime` only for adaptive field assembly and
the same `GpuDepondtIntegrator` calls for fine atoms. Its inputs differ only
by enabling adaptive configuration and setting every block to fine.

## Staged harness

`gpu_adaptive_runtime_benchmark --parity-only --require-acceptance` runs the
following hierarchy on exchange, DMI, uniaxial, combined, and finite-temperature
fixtures:

1. initial Hamiltonian field;
2. bilinear, on-site, and total term-resolved energies;
3. nonzero initial torque and displacement assertions;
4. exact thermal-field identity after the predictor draw;
5. Depondt predictor state;
6. Hamiltonian field at the predicted state;
7. exact second thermal-field reuse (Depondt does not draw again);
8. Depondt corrector/final spin;
9. six-step finite-temperature fixed-seed trajectories using one and two
   ensembles. The adaptive arm copies the feature-off thermal field, so the
   two-ensemble case validates indexing and state parity without requiring
   independently constructed cuRAND generators to produce identical streams.

Hamiltonian and trajectory comparisons use the established fp64 roundoff
budget; thermal fields require bitwise identity. Each stage reports bitwise
status, maximum absolute/relative error, and the first failing stage.

## Negative controls

The CMake fault-injection options are intentionally disabled by default:

| Mutation | Expected first failure |
|---|---|
| `RCG09A_NEGATIVE_DMI_SIGN` | initial Hamiltonian field/term energy |
| `RCG09A_NEGATIVE_NO_TRANSPOSE` | initial Hamiltonian field/term energy |
| `RCG09A_NEGATIVE_THERMAL_AMPLITUDE` | thermal field |
| `RCG09A_NEGATIVE_DAMPING` | T=0 Depondt predictor |
| `RCG09A_NEGATIVE_RNG_DISPLACEMENT` | finite-T thermal field |

These are build-time mutations only; no mutation is enabled in the production
source or acceptance target. They must be run as separate fault-injection
builds and must return nonzero after reporting the named stage.

## Local execution

The CUDA/fp64 target compiled successfully with:

```text
cmake --build build_rcg09a_cuda_fp64 --target gpu_adaptive_runtime_benchmark -j2
```

The direct parity command returned the project’s unavailable-backend code 77
on the local host at one point:

```text
ADAPTIVE-ASD-PARITY unavailable: no backend device
```

The CTest acceptance target did subsequently run with a visible CUDA device
and passed after the staged harness reused the feature-off production thermal
field. The DMI-sign mutation was also executed and failed at the initial
Hamiltonian field/energy stage as intended. The other four mutation builds
compiled, but their runtime executions were skipped when the device became
unavailable; they are not reported as passing by inference. CPU and HIP are
not reported as passing by inference.

The staged implementation now includes an explicit two-ensemble acceptance
fixture. It uses the feature-off thermal field as the sole stochastic oracle
and copies that field into the adaptive integrator, which is the parity
contract and avoids treating independent cuRAND stream identity as a physics
requirement. Runtime execution of this added case still requires an available
CUDA or HIP device; the current container has neither.

## Validation record (2026-08-15)

- CPU Release build: passed.
- CPU coarse-graining suite: 29/29 tests passed, including production ASD,
  setup rejection, moving parity, DMI, and trajectory tests.
- CUDA/fp64 build: `gpu_adaptive_runtime_benchmark` compiled successfully.
- CUDA runtime: direct parity execution returned 77 with
  `ADAPTIVE-ASD-PARITY unavailable: no backend device`; CTest consequently
  skipped the acceptance test through its documented `SKIP_RETURN_CODE`.
- HIP: not installed/available on this host and not claimed as passing.
