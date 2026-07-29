# CG-10.5 production integration

## Production entrypoint inventory

The normal executable enters `setup_simulation` in `uppasd.f90`. Complete
input, geometry, moments, damping, and Hamiltonian setup precede
`setup_adaptive_cg_production`; the adaptive capability check therefore sees
the resolved solver and all source arrays. Setup happens before legacy
Hamiltonian input storage is released and before GPU preflight.

CPU measurement-phase spin dynamics runs through `sd_mphase` in
`sd_driver.f90`. An enabled run dispatches each complete Heun step through
`adaptive_cg_cpu_step`; a disabled run reaches the unchanged legacy field and
integrator statements. Measurements consume `emom`/`emomM` before each step,
and reconstruction refreshes both before the next measurement.

GPU spin dynamics enters `FortranData_Initiate`, `GpuSimulation::initiateMatrices`,
and `GpuSDSimulation::SDmphase`. The Fortran setup seam stages production
topology, mutable state, projection, bond, selector, and tensor arrays before
memory preflight. `GpuAdaptiveRuntime` validates and copies them, after which
the staging pointers are cleared. The enabled measurement loop skips the
ordinary all-atom short-range Hamiltonian and integrator and calls the
adaptive Heun, reconstruction, selector, publication, and compaction path.

Normal cleanup prints the adaptive summary before moment/topology owners are
destroyed. Assigning the production state to its default value recursively
releases all owned Fortran allocatables; the GPU owner releases all tracked
device tensors. Setup failures unwind the Fortran owner, and GPU staging is
cleared on preflight/allocation exceptions.

## Capability boundary

The initial production model accepts regular replicated cells with exact
positive block divisibility, `P P P` boundaries, mode `S`, no initial phase,
deterministic Heun (`SDEalgh=1`), zero temperature, one ensemble-compatible
ferromagnetic dynamical channel, scalar exchange, uniform Landé factor and
damping, and no restart. CPU and GPU call the same Fortran validation before
either runtime is allocated.

The setup diagnostic names the offending keyword or capability when it
rejects Monte Carlo, GNEB, spin-lattice or other modes; initial-phase or
restart input; stochastic/finite-temperature dynamics; DMI, tensor exchange,
anisotropy, dipole, external/time-dependent fields, sparse/reduced/fixed
moments, energy-output paths not yet connected to hybrid accounting,
nonperiodic or explicit-device geometry, multiple channels, or heterogeneous
Landé/damping data. FFT dipole coupling is intentionally outside this first
production boundary rather than silently omitted.

The ordinary text `inpsd.dat` reader is the runtime input frontend. It maps
every CG keyword to `adaptive_cg_config_t`. JSON/YAML files in the repository
are conversion and test artifacts rather than independent runtime readers;
they must generate the documented text keywords before invoking UppASD.

## Ownership and synchronization

`AdaptiveCGProduction` owns topology, material metadata and diagnostics,
tensor/projection descriptors, unique production bonds, masks, selector and
hybrid runtime, atom moments, coarse state, and GPU staging storage for the
whole simulation phase. No test fixture supplies production coefficients or
maps. Auxiliary files are read during enabled setup and closed immediately.
With `do_adaptive_cg=N`, setup returns before opening a file or constructing
any owner.

Both CPU and GPU publish selector changes only after a complete corrector
stage. Coarse reconstruction is committed to the atomic direction array used
by ordinary measurements and output. GPU arrays remain alive through
preflight and upload; the C staging sentinel is cleared immediately after the
GPU owner has copied them.

## Executable evidence

`tests/coarse_graining/run_production_e2e.py` launches the normal binary with
real `inpsd.dat` files for feature-off, all-fine, all-coarse, static mixed,
adaptive transition, invalid-block, CUDA/HIP static-mixed, and CUDA/HIP
adaptive cases. CPU assertions compare active updates and prove that mixed
and coarse modes reduce short-range/integration work. GPU cases assert
compact active counts and are skipped only when the compiled backend reports
that no device is present.
