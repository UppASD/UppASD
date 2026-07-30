# CG-13 release validation

CG-13 is the release gate for the accepted adaptive coarse-graining model. The
suite runs the ordinary test executables and, for production cases, the normal
UppASD executable with an ordinary `inpsd.dat`. It does not use test-only
staged moments or bypass the production capability preflight.

## Running the suites

Configure a CPU, CUDA, or HIP build with `BUILD_TESTING=ON`, then build the
test targets. The CTest labels are the stable suite interface:

```sh
ctest --test-dir build_cpu -L cg13-cpu --output-on-failure
ctest --test-dir build_cuda -L cg13-cuda --output-on-failure
ctest --test-dir build_hip  -L cg13-hip  --output-on-failure
```

`coarse-graining` is a broader label for the component tests. `cg13-cpu` is
always available. The CUDA and HIP labels are added only by the matching
backend build, so a CUDA build cannot accidentally report HIP validation (or
vice versa). Hardware-dependent benchmark binaries are built by the backend
build but are not CTest pass/fail tests; record their stdout as an artifact.

The production executable suite includes feature-off, all-fine, all-coarse,
static mixed, adaptive transitions, DMI plus uniaxial anisotropy, initial
texture preparation, and CPU/GPU parity. The setup-rejection matrix runs each
input mutation in a fresh temporary directory and requires a nonzero exit
before the measurement phase. This keeps rejected combinations from being
silently accepted by a later fallback. If the ordinary input reader rejects an
inconsistent combination before adaptive preflight, that setup error is also
the tested negative result; cases reaching adaptive preflight retain its
named diagnostic.

## Reference cases and provenance

| Case | Executable evidence | Reference/provenance |
| --- | --- | --- |
| Spin spiral | `coarse-graining-smooth-projected-operator`, `coarse-graining-tensor-operator`, `initmag_spin_spiral` | Analytic periodic phase textures and the normal `Initmag 8` path |
| Domain wall | `coarse-graining-tensor-operator`, `coarse-graining-static-hybrid-operator`, adaptive solver | `sqrt(A/K)` width and translated wall fixtures |
| Skyrmion | static hybrid operator and adaptive solver | Translated radial texture and adaptive transition fixture |
| Static interface | static hybrid operator | Uniform patch, energy ownership, tangent derivative, and refinement checks |
| Adaptive transition | adaptive solver and production `adaptive_mixed` | Complete-step publication, hysteresis/dwell, rollback, and compact-list rebuild |
| Multi-channel | `coarse-graining-multichannel-tensor-operator` and channel restriction tests | Two-sublattice compensated AFM/Ferri reference; production single-channel gate remains explicit |
| Dipole | GPU `gpu_fft_static_mixed` when enabled | Periodic EWALD3D FFT field/energy coupling; short-range mask is independent |

The analytic fixtures live in `tests/coarse_graining/test_*.f90` and are
versioned with the operator implementation. Production fixtures live under
`tests/coarse_graining/e2e`; their `posfile`, interaction files, masks, and
expected diagnostic predicates are the reproducible reference inputs. The
suite compares CPU/GPU decisions, per-term energies, field checksums, and
restart trajectories where a backend is enabled.

## Precision and diagnostics

The default release reference is fp64. Analytic CPU tolerances are listed in
the test sources and are normally `1e-12` relative for energies and tensor
coefficients, with tighter zero-torque checks where appropriate. GPU fp64
parity permits the documented accumulated-kernel tolerance (`5e-4` relative
for energies and `8e-4` relative for field checksums). fp32 is a separate
backend/precision result and must use its own recorded tolerances; it is not
allowed to inherit fp64 claims.

Set `cg_diagnostics` to 1, 2, or 3 to report progressively more of the
resolved operator, ownership states, selector transitions, energy terms,
field/trajectory checksums, phase timings, FFT timing, and device memory.
Level 0 suppresses adaptive reports. Diagnostics are evidence, not a second
physics model.

## User model and input boundary

The production path is an optional low-temperature, fixed-length, two-scale
spin-dynamics model. It accepts a regular replicated periodic cell with exact
positive block divisibility, deterministic fixed-length Heun LLG at `T=0`,
one ferromagnetic dynamical channel, uniform Landé factor and damping, scalar
exchange, DMI, and deterministic cell-periodic type-1 uniaxial anisotropy.
The accepted initial phase can be ordinary spin dynamics, Metropolis, heat
bath, single-Q, two-Q/cube, three-Q, or VPO preparation; those stages finish
atomistically before adaptive ownership is constructed. GPU periodic dipole
coupling is additionally available as `gpu_dipole_mode EWALD3D_FFT` with
`do_gpu=Y` and `do_gpu_llg=Y`.

`alat` is mandatory and must be an explicit finite positive SI length in
metres. It is the characteristic length used to convert neighbour
displacements to material tensors. Exchange stiffness scales as `alat^-1`,
spiralization as `alat^-2`, and anisotropy density as `alat^-3`; omitting it
would silently change the physical model, so setup rejects the historical
implicit default. `block_size_x/y/z` partitions repetition-cell coordinates,
not arbitrary atom indices, and every repetition count must divide exactly.

Tensor exchange, higher-order pair/multisite interactions, random or
unsupported anisotropy, legacy dipoles, nonperiodic boundaries, explicit
device geometry, heterogeneous channel/damping data, finite-temperature
measurement dynamics, external/time-dependent fields, sparse/reduced/fixed
moment paths, legacy samplers, and unsupported runner combinations are
rejected during setup. `Initmag 4` restart input is rejected because adaptive
ownership, selector age/epoch, and reconstruction state are not serialized.
There is no implicit fallback to an atomistic run under the adaptive flag.

Low-temperature cone reconstruction is a deterministic geometric
reconstruction that preserves a requested coarse resultant. It is not
thermodynamic equilibration, and no finite-temperature equilibration claim is
made for it.

## Error budget and known limitations

The continuum tensor is a long-wavelength approximation. Error grows when a
texture varies on the scale of a block, at sharp interfaces, for high-q
spirals, and for atomistic modes not represented by the accepted channel
descriptor. Interface error comes from projected coarse ghosts, finite
buffer/halo width, block-centre quadrature, and discretized ownership; the
refinement and translated wall/skyrmion tests are regression checks, not a
promise of zero error for arbitrary textures.

`EWALD3D_FFT` uses a uniform basis-resolved periodic macro grid. Its dipole
error includes the finite grid/padding and reciprocal truncation/tolerance,
basis deposition/interpolation, periodic-image convention, and the macro
field approximation. The adaptive mask does not refine that FFT grid.
Dipole diagnostics must therefore be reported separately from short-range
interface error.

## Performance reporting

Run the backend benchmark with several repetitions and retain raw stdout, for
example:

```sh
build_cuda/bin/gpu_adaptive_runtime_benchmark --repetitions 7 \
  --require-acceptance > cg13-cuda-benchmark.txt
```

Each sweep row reports requested fraction, active atom count and fraction,
active block count and fraction, interface count and fraction, wall median and
MAD, and atomistic/coarse/interface/integration/FFT phase times. The overhead
line reports selector device time and host wall time, compaction device time,
localized host wait and wall time, uploaded mask bytes, device allocation, and
the combined overhead relative to the mixed-resolution step. The separate
`adaptive-step` line measures a complete two-stage Heun call and reports the
same phase decomposition plus unaccounted wall time. A crossover or speedup
without these active fractions and overheads is not a CG-13 performance
result.
