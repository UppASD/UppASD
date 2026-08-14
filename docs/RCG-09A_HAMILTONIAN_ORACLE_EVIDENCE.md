# RCG-09A.1 adaptive all-fine Hamiltonian oracle

Status: complete for the RCG-09A.1 Hamiltonian scope. Exchange, DMI, and deterministic uniaxial parity executed on CUDA/fp64; required DMI negative controls executed on CUDA/fp64. HIP is unavailable on the validation host and is explicitly not reported as passing.

## Production oracle and fixture

`tests/coarse_graining/benchmark_gpu_adaptive_runtime.cpp` now executes a Hamiltonian oracle before its performance sweep. The reference constructs `GpuHamiltonianCalculations` and calls `heisge(..., measure=true)`, the feature-off GPU ASD Hamiltonian entry point. The adaptive arm constructs `GpuAdaptiveRuntime` with every block fine. Neither arm uses an analytical field oracle.

The fixture has eight periodic blocks, two atoms per block, one ensemble, unit moments, deliberately non-collinear atom directions, scalar exchange `J=0.1`, and reciprocal DMI `D_ij=(0.37,-0.29,0.61)`, `D_ji=-D_ij`. It requires nonzero DMI field at both endpoints of its first canonical bond. The uniaxial evaluation uses a distinct normalized axis and nonzero `(K1,K2)` on every atom. It reports max absolute/relative field errors and exchange, DMI, uniaxial, and total energies. Production energy columns are per atom and are multiplied by `N`; adaptive energy is total J. Difference-of-evaluations isolates the DMI and uniaxial field terms, while `energyM(:,2)` is the production uniaxial term diagnostic.

The existing target is CUDA/fp64 (`real=double`); adaptive energy accumulation and production `energyM` are fp64. Acceptance is `1e-12` relative in fp64 (`5e-5` only in an explicitly fp32 build).

### CUDA/fp64 execution evidence

On `alcazar`, the normal 16-atom (`8x2`) run passed the production oracle:

| Quantity | Adaptive | Production `heisge` | Difference |
| --- | ---: | ---: | ---: |
| Exchange energy | -2.14250853461587454 | -2.14250853461587365 | 8.9e-16 |
| DMI energy | 6.45635050280112477e-1 | 6.45635050280111922e-1 | 5.6e-16 |
| Total energy | -1.49687348433576206 | -1.49687348433576162 | 4.4e-16 |
| Exchange field max relative | — | — | 1.87325218057106195e-16 |
| DMI field max relative | — | — | 2.13735599136465424e-16 |
| Total field max absolute / relative | — | — | 2.49800180540660222e-16 / 2.11180000688015182e-16 |

The fixture reported `nonzero_endpoint_DMI=true`; it is therefore not a zero-field or one-endpoint false pass. The exchange-only performance harness independently retained its prior all-fine parity (`1.17e-16` relative).

### CUDA/fp64 uniaxial execution

The expanded normal fixture passed on `alcazar` with `negative_control=off result=PASS`:

| Quantity | Adaptive | Production `heisge` | Difference |
| --- | ---: | ---: | ---: |
| Uniaxial energy | -1.31668151511199882 | -1.31668151511199882 | 0.0 |
| Total energy including uniaxial | -2.81355499944776088 | -2.81355499944776044 | 4.4e-16 |
| Exchange field max absolute / relative | — | — | 5.55111512312578270e-17 / 1.87325218057106195e-16 |
| DMI field max absolute / relative | — | — | 2.49800180540660222e-16 / 2.13735599136465424e-16 |
| Uniaxial field max absolute / relative | — | — | 1.11022302462515654e-16 / 4.27942472688887573e-16 |
| Total exchange+DMI field max absolute / relative | — | — | 2.49800180540660222e-16 / 2.11180000688015182e-16 |
| Total field including uniaxial max absolute / relative | — | — | 2.49800180540660222e-16 / 2.20415714914710008e-16 |

## Directed-list folding contract

Production `GpuHamiltonianCalculations::Heisge::each` and `ApplyHamiltonian:heisge` consume a directed entry as `D_ij x M_j`. Adaptive stores the matrix for `D_ij x s`: its first endpoint gets `K s_j` and its second gets `K^T s_i`. With `D_ji=-D_ij`, `K^T=-K`, so the second scatter is exactly `D_ji x s_i`. Scalar exchange needs `J_ji=J_ij`.

`build_unique_bonds` now validates these reciprocal conditions before folding. A missing reverse entry, unequal exchange, or DMI vector that is not its reverse partner is rejected with a named reciprocal-contract diagnostic. `coarse-graining-adaptive-hamiltonian-contract` calls those exact setup predicates and proves that a deliberately nonreciprocal DMI list fails; it also checks unequal `J` and an aliased DMI pair.

The sparse production list carries target atom indices but no periodic-image translation. Two images mapping to one `(i,j)` cannot be given a correct distinct key or displacement for adaptive buffer ownership. Adaptive setup therefore rejects more than one directed entry in either direction for a canonical exchange/DMI pair with `periodic-image alias`; this is an explicit capability boundary rather than unsafe merging.

## Ensemble moments

UppASD permits `mmom(atom,ensemble)`, while adaptive currently stores one `atom_moment_mub(:)` and one folded bond matrix sourced from ensemble one. Setup now rejects an ensemble that differs from ensemble one by more than 64 scaled fp64 eps, naming atom, ensemble, and both values. The same contract test exercises both the accepting identical-moment case and the rejecting differing-moment case. It therefore no longer silently uses ensemble-one moments for a distinct-moment ensemble.

## Hamiltonian capability table

No row is unknown.

| Production contribution | Adaptive classification | Basis |
| --- | --- | --- |
| Isotropic Heisenberg exchange | SUPPORTED_AND_PROVEN | CUDA/fp64 production-oracle energy and field parity above. |
| Dzyaloshinskii--Moriya interaction | SUPPORTED_AND_PROVEN | CUDA/fp64 term-resolved production-oracle parity above; both directed-list and transpose fault controls fail clearly. |
| Deterministic uniaxial on-site anisotropy | SUPPORTED_AND_PROVEN | CUDA/fp64 `heisge` anisotropy diagnostic and difference-field parity above. |
| Cubic or other non-type-1 on-site anisotropy | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `build_production_anisotropy` accepts only `taniso=1` uniaxial sites. |
| Tensorial exchange | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `ham_inp%do_jtensor`. |
| Symmetric anisotropic exchange | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `ham_inp%do_sa`. |
| Pseudo-dipolar interaction | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `ham_inp%do_pd`. |
| Biquadratic exchange | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `ham_inp%do_bq`. |
| Biquadratic DMI | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `ham_inp%do_biqdm`. |
| Ring, chiral, and four-spin terms | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `do_ring`, `do_chir`, `do_fourx`. |
| External, pulse, and demagnetizing fields | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `hfield`, `do_bpulse`, `demag`. |
| Legacy atomistic/macrocell dipole | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `do_dip`. |
| GPU periodic FFT dipole | SUPPORTED_AND_PROVEN | This is not an atomistic `heisge` short-range contribution. Adaptive execution delegates its all-grid field to the production `GpuHamiltonianCalculations::evaluateAdaptiveFftDipole` / `dipoleConvolution` owner; no legacy `do_dip` comparison is relevant. |
| Finite-temperature stochastic field | EXPLICITLY_REJECTED_BY_ADAPTIVE_SETUP | `Temp /= 0` is rejected before runtime setup. |

## Negative controls

The checked-in source is defect-free. The following fault-injection builds ran on CUDA/fp64 and both emitted the harness's expected `negative_control=enabled result=PASS`, where PASS means that parity failed clearly.

| Defect injected | DMI energy adaptive / production | DMI field max relative | Total field max relative | Result |
| --- | ---: | ---: | ---: | --- |
| `RCG09A_NEGATIVE_DMI_SIGN` | -6.45635050280110701e-1 / 6.45635050280111922e-1 | 2.0 | 1.97608635661279308 | Clear failure |
| `RCG09A_NEGATIVE_NO_TRANSPOSE` | 6.45635050280112477e-1 / 6.45635050280111922e-1 | 1.09084236991678885 | 1.03812398009028106 | Clear failure |

The missing-transpose build subsequently reports failure in the exchange-only performance benchmark. That is expected collateral: the injected defect removes the second endpoint scatter for every bilinear matrix, including the exchange part. It is not evidence against the defect-free normal run.

## Completion procedure

On the CUDA host, rebuild and run the expanded normal fixture:

```bash
cmake -S . -B build_rcg09a_cuda_fp64 -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=DOUBLE
cmake --build build_rcg09a_cuda_fp64 --target gpu_adaptive_runtime_benchmark -j2
./build_rcg09a_cuda_fp64/bin/gpu_adaptive_runtime_benchmark --blocks 8 --atoms-per-block 2 --warmup 1 --iterations 1 --repetitions 3 --require-acceptance
```

The recorded run reports `negative_control=off result=PASS` and roundoff-scale uniaxial field/energy differences. It is intentionally CUDA-only: this host has no HIP toolchain/runtime. The previously recorded DMI-sign and missing-transpose CUDA fault builds remain the required negative-control evidence.
