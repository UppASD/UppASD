# Adaptive exchange-CG purge postmortem

## Scope

The retired feature was the experimental **adaptive exchange coarse-graining**
implementation. This purge did not remove the successful coarse-block dipole
projection/integration work; `coarse` in that dipole path describes the
surviving macrocell/block projection.

## Branch lineage and purge history

The historical implementation remains at
`gpu_hip_cu_ab_cg` @ `438d14d721a2aac82745413000a84501f50b0f8e`.
The cleaned forward branch is `gpu_hip_cu_ab_dip`.

The preceding purge commits were:

- PURGE-01 `ae144521e5acf27f6287ef3c6444b2a8bbb2000a`
- PURGE-02 `9cb691120a8558d4a508114bc3203266f1a82408`
- PURGE-03 `d204400d9a98fef3e82d5d0274a7d01aa363cb31`
- PURGE-04 `9941156448f8df8d9f12aca613e07f5205f607d2`
- PURGE-05 `a2ed89dcd1135bfabbb5d76daeeb020a4e37757e`

PURGE-06 is the separate closure commit containing this postmortem and the
final historical-document and CI disposition.

## Removed at subsystem level

- adaptive Fortran operators and runtime;
- adaptive input configuration and its 18 retired keywords;
- the Fortran-to-C++ adaptive bridge;
- GPU adaptive runtime and reconstruction support;
- adaptive-only active-subset integration and thermal helpers; and
- adaptive tests and examples.

## Retained

The cleaned branch deliberately retains ordinary atomistic ASD, ordinary
production Depondt, CUDA/HIP infrastructure, CPU-Hamiltonian/convolution,
OPEN_FFT, EWALD3D, macrocell/coarse-block dipoles, `BlockTopology`,
`macrocells`, and `hamiltonianmacroblocks`. It also retains the independent
DMI correctness fixes and the DMI/BlockTopology regression oracles. The
[DMI evidence](RCG-02_DMI_HANDEDNESS_EVIDENCE.md) remains active, while the
adaptive-CG design/evidence set is preserved in
[`docs/attic/adaptive-cg/`](attic/adaptive-cg/).

The active WP10 dipole records remain in the docs tree:
`luna_wp10_7c_open_coarse_acceptance.md`,
`sol_wp10_7a_open_coarse_projection.md`, and
`terra_wp10_7b_open_coarse_integration.md`.

## CI disposition

The obsolete adaptive-CG clean workflow was repurposed as
`.github/workflows/clean-cpu-regression.yml`. It performs a clean CPU
configure/build and runs the surviving CPU-HAM, DMI dimer, BlockTopology, and
host-side dipole tests. It does not require CUDA/HIP hardware or refer to
removed adaptive tests.

## Validation record

The final PURGE-06 clean CPU production build used GNU Fortran 16.2.0 and
AppleClang 21.0.0.21000101 with `BUILD_TESTING=OFF` and
`UPPASD_GPU_BACKEND=OFF`; configuration and build completed at 100%.

A separate FFTW-enabled `BUILD_TESTING=ON`, `UPPASD_GPU_BACKEND=OFF`
configuration and build also completed at 100%. `ctest -N` reported 16 local
registrations. Fourteen passed, and the only two failures were the YAML-driven
`regression-test` and `asd-tests`, both blocked by the local Python
`ModuleNotFoundError: No module named 'yaml'`. No purge-induced code/test
failure appeared. The four available host dipole tests all passed.

Deterministic DMI, BlockTopology, and CPU-HAM-11 stdout hashes remained:

```text
DMI           13d26598b2c29b9809250943583218e4dc48c9bcda2ffcf2a4d2b74a1d0fc80b
BlockTopology 24bc68909ac0b8170d1196a3961fd1a7457e11627b9140ecef24b9663ac3772f
CPU-HAM-11    b33b000f51467d181b624951528466bbcacf67bb27313fc69096eaffc6c32646
```

GPU compilation was not performed: neither `nvcc` nor `hipcc` was available.
This is a validation limitation, not a code failure.

The final active source/build/test sweep found no adaptive implementation
symbols, retired input keywords, deleted-source references, or adaptive
CMake/CI dependency. It did find eight retained `BlockTopology`
comments/diagnostic strings that
still say “Adaptive coarse graining”; these are non-executable diagnostics in
surviving topology validation and were left unchanged because PURGE-06 is
restricted to documentation and CI.

The pre-purge ownership and evidence audit is committed as
[`PURGE-00_AUDIT.md`](PURGE-00_AUDIT.md). Its baseline data is intentionally
unchanged; it records the state audited at `438d14d7`.
