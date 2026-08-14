# RCG-09A.4 adaptive all-fine ASD parity

Status: harness implemented; accelerator execution is pending on a host with
an available CUDA or HIP device. The local validation host built the CUDA/fp64
target but reported no usable backend device, so this record does not claim
the acceptance checklist is closed.

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
9. six-step finite-temperature fixed-seed trajectory, using one ensemble.

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

The staged implementation also has an explicit two-ensemble fixture path, but
the current CUDA/fp64 host shows non-reproducible cuRAND normal-field reuse
between the independently constructed production oracles for `M=2`. The
single-ensemble finite-temperature path is therefore the acceptance fixture;
multi-ensemble parity remains an identified follow-up rather than being
silently relaxed.
