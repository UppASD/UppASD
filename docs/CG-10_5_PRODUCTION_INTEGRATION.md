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
positive block divisibility, `P P P` boundaries, mode `S`, deterministic Heun
(`SDEalgh=1`), zero measurement temperature, one ensemble-compatible
ferromagnetic dynamical channel, scalar exchange, uniform Landé factor and
damping, and no restart. `Initmag` modes 1 (random), 2 (cone), 3 (momfile), 5
(random Ising seed), and 8 (spin spiral) are accepted. `ip_mode N` starts CG
from that state; `ip_mode S` first runs the ordinary atomistic SD initial
phase and hands its completed texture to CG. CPU and GPU use the same Fortran
capability preflight before the initial phase can allocate a device.

The setup diagnostic names the offending keyword or capability when it
rejects Monte Carlo, GNEB, spin-lattice or other measurement modes;
non-SD initial-phase modes or restart input; stochastic/finite-temperature
measurement dynamics; DMI, tensor exchange, anisotropy, legacy `do_dip`,
external/time-dependent measurement fields,
sparse/reduced/fixed moments, energy-output paths not yet connected to hybrid
accounting, nonperiodic or explicit-device geometry, multiple channels, or
heterogeneous Landé/damping data. The GPU periodic path accepts
`gpu_dipole_mode EWALD3D_FFT`; CPU adaptive runs and `OPEN_FFT` remain outside
this periodic production boundary and are rejected explicitly.

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

Adaptive setup has two lifecycle points. A read-only preflight validates the
resolved inputs, geometry, Hamiltonian, and selected initial/measurement
solvers before `run_initial_phase`. Actual topology, material, selector, and
runtime construction occurs at the start of `run_measurement_phase`, after an
optional atomistic `ip_mode S` has updated `emom`, `emomM`, and `mmom`.
Hamiltonian input storage is retained until this handoff construction and
released immediately afterward. No coarse ownership or adaptive selector is
active during the initial phase.

Both CPU and GPU publish selector changes only after a complete corrector
stage. Coarse reconstruction is committed to the atomic direction array used
by ordinary measurements and output. GPU arrays remain alive through
preflight and upload; the C staging sentinel is cleared immediately after the
GPU owner has copied them.

## Field and energy ownership

Short-range scalar exchange bonds are owned once: a bond touching an
atomistic/interface block is evaluated by the atomistic operator with
projected coarse ghosts, while wholly coarse stencil terms are evaluated by
the coarse tensor operator. The reported atomistic-bilinear and coarse
exchange/spiralization/anisotropy/external terms therefore form one hybrid
total rather than overlapping atomistic and continuum totals.

`EWALD3D_FFT` is deliberately independent of that short-range ownership mask.
At each predictor and corrector evaluation, the production FFT owner reduces
the current atomic moment vectors to the uniform basis-resolved macro grid
and convolves it. For the accepted single dynamical FM channel:

- atomistic and interface atoms receive their basis channel's field once;
- a coarse block receives the magnetic-moment-weighted basis-field average
  in its one dynamical equation once;
- dipole energy uses `-1/2 mu_B sum_i mu_i m_i dot B_i`, with actual atomic
  directions in atomistic/interface blocks and the coarse channel direction
  as the source in coarse blocks.

The field prefactor is applied exactly once when the padded convolution field
enters adaptive equations. The adaptive mask never changes FFT grid
resolution or drops a requested dipole interaction.

## Diagnostics

`cg_diagnostics` levels 1 through 3 report the resolved operator, mask,
selector, reconstruction, thresholds/cadence, block topology, basis-to-FFT
and single-channel mapping, initial/final ownership counts, interface and
active work, and owned CPU/device bytes. Step/summary output additionally
reports accepted and rejected transition counts, transition selector reason,
before/after/jump energies on the CPU transition path, named hybrid energy
terms, atom/coarse field checksums, trajectory checksums, and accumulated
field, integration, reconstruction, selector, compaction, and FFT timings.
Level 0 suppresses these adaptive reports.

## Executable evidence

`tests/coarse_graining/run_production_e2e.py` launches the normal binary with
real `inpsd.dat` files for feature-off, all-fine, all-coarse, static mixed,
adaptive transition, invalid-block, CUDA/HIP static-mixed, and CUDA/HIP
adaptive cases. Nonuniform seeded CPU/GPU fixtures compare state decisions,
per-term energies, atom/coarse field checksums, transition counts, and final
restart trajectories. A GPU `EWALD3D_FFT` mixed-state fixture proves that the
production convolution contributes a nonzero adaptive dipole energy and FFT
timing. CPU assertions compare active updates and prove that mixed and coarse
modes reduce short-range/integration work. GPU cases assert compact active
counts and are skipped only when the compiled backend reports that no device
is present. Additional executable cases cover `Initmag=2` plus CPU/GPU
atomistic-SD handoff and the deterministic `Initmag=8` spin-spiral texture.

The ordinary inputs in `examples/AdaptiveCoarseGraining` cover a static mixed
mask and adaptive MAX_ANGLE selection. They contain no test-only hooks or
manually staged arrays and are executed by the same e2e test.
