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
ferromagnetic dynamical channel, scalar exchange, DMI, deterministic type-1
uniaxial anisotropy, uniform Landé factor and damping, and no restart.
`Initmag` modes 1 (random), 2 (cone), 3 (momfile), 5
(random Ising seed), and 8 (spin spiral) are accepted. `ip_mode N` starts CG
from that state. Spin-only atomistic preparation accepts `ip_mode S`
(spin dynamics), `M` (Metropolis), `H` (heat bath), `Q` (single-Q search),
`Y` (three-Q search), `Z` (two-Q/cube search), and `G` (VPO energy
minimization), then hands the completed texture to CG. CPU and GPU use the
same Fortran capability preflight before the initial phase can allocate a
device. `X` and `SX` replica-exchange runners remain outside this boundary.

The setup diagnostic names the offending keyword or capability when it
rejects Monte Carlo, GNEB, spin-lattice or other measurement modes;
unsupported initial-phase modes or restart input; stochastic/finite-temperature
measurement dynamics; tensor exchange, unsupported anisotropy kinds, legacy `do_dip`,
external/time-dependent measurement fields,
sparse/reduced/fixed moments, energy-output paths not yet connected to hybrid
accounting, nonperiodic or explicit-device geometry, multiple channels, or
heterogeneous Landé/damping data. The GPU periodic path accepts
`gpu_dipole_mode EWALD3D_FFT`; CPU adaptive runs and `OPEN_FFT` remain outside
this periodic production boundary and are rejected explicitly.

The short-range capability now includes scalar exchange, DMI, and
deterministic type-1 uniaxial anisotropy. One or two uniaxial axes per site
are accepted; cubic/type-2, combined type-7, random, and more general
anisotropies remain explicit rejections. Tensor exchange and the other
higher-order pair or multisite interactions remain outside the boundary.

`alat` is mandatory whenever `do_adaptive_cg=Y`. It must be explicitly
specified in `inpsd.dat` as a finite positive length in metres. The historical
default is not accepted because the continuum tensors have different length
dimensions and an implicit reduced-unit lattice parameter would silently
produce the wrong material constants.

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
optional supported atomistic initial phase has updated `emom`, `emomM`, and
`mmom`. A runner-independent handoff check rejects non-finite, zero, non-unit,
or mutually inconsistent moment state, then canonicalizes `emomM`, `emom2`,
`mmom2`, and `mmomi` from the accepted direction and magnitude arrays.
Hamiltonian input storage is retained until this handoff construction and
released immediately afterward. No coarse ownership or adaptive selector is
active during the initial phase.

Host `M/H/Q/Y/Z/G` preparation can hand off to a CUDA/HIP measurement because
the measurement device is initialized afterward. GPU `M/H` preparation is
also accepted when the initial phase does not request the measurement-only
adaptive FFT dipole. With `gpu_dipole_mode EWALD3D_FFT`, use `do_gpu_mc=N` so
the MC preparation remains on the host before the FFT owner is constructed.

Both CPU and GPU publish selector changes only after a complete corrector
stage. Coarse reconstruction is committed to the atomic direction array used
by ordinary measurements and output. GPU arrays remain alive through
preflight and upload; the C staging sentinel is cleared immediately after the
GPU owner has copied them.

## Field and energy ownership

Short-range scalar exchange and DMI bonds are owned once: a bond touching an
atomistic/interface block is evaluated by the atomistic operator with
projected coarse ghosts, while wholly coarse stencil terms are evaluated by
the coarse tensor operator. Deterministic uniaxial anisotropy is similarly
evaluated per atom in atomistic/interface blocks and as an energy density in
coarse blocks. The reported atomistic-bilinear, atomistic-onsite, and coarse
exchange/spiralization/anisotropy/external terms therefore form one hybrid
total rather than overlapping atomistic and continuum totals.

## SI material conversion

The production harness calls the typed CG-03 extractor directly on the
resolved UppASD neighbour lists. With directed pair energies
\(J_{ij}=\mu_B\mu_i\mu_j\,\mathtt{ncoup}_{ij}\), DMI vectors
\(\mathbf L_{ij}=\mu_B\mu_i\mu_j\,\mathtt{dm\_vect}_{ij}\), physical
displacements \(\mathbf r_{ij}=\mathtt{alat}\,\Delta\mathbf x_{ij}\), and
material-cell volume \(V\), the zero-regularization limits are

\[
A_{pq}={1\over4V}\sum_{ij}J_{ij}r_{ij,p}r_{ij,q},
\qquad
D_{kp}={1\over2V}\sum_{ij}L_{ij,k}r_{ij,p}.
\]

Thus \(A\) is in J m\(^{-1}\) and scales as `alat^-1`, while the
spiralization \(D\) is in J m\(^{-2}\) and scales as `alat^-2`. For each
accepted uniaxial axis, the resolved atomistic coefficients are converted
back to joules,

\[
K_{1,i}^{J}=\mu_B\,\mathtt{kaniso}_{1,i}\mu_i^2,\qquad
K_{2,i}^{J}=\mu_B\,\mathtt{kaniso}_{2,i}\mu_i^4,
\]

and summed over the material cell before division by \(V\), giving coarse
coefficients in J m\(^{-3}\), scaling as `alat^-3`. Both resolutions use

\[
E_i=K_{1,i}^{J}c_i^2+
K_{2,i}^{J}(2c_i^2-c_i^4),\qquad c_i=\mathbf m_i\cdot\mathbf e_i,
\]

with the field obtained from
\(\mathbf B_i=-(\mu_B\mu_i)^{-1}\partial E_i/\partial\mathbf m_i\).
The DMI continuum density is
\(\sum_{kp}D_{kp}[\mathbf m\times\partial_p\mathbf m]_k\), retaining the
CG-01/CG-03 handedness convention.

These long-wave moment formulas and their spin-spiral interpretation follow
the micromagnetic multiscale derivations in Zimmermann et al.,
*Phys. Rev. B* **99**, 214426 (2019),
https://doi.org/10.1103/PhysRevB.99.214426, and Schweflinghaus et al.,
*Phys. Rev. B* **94**, 024403 (2016),
https://doi.org/10.1103/PhysRevB.94.024403. The microscopic antisymmetric
pair interaction is the Moriya form,
https://doi.org/10.1103/PhysRev.120.91.

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
adaptive transition, mixed DMI/anisotropy, missing-`alat`, invalid-block,
CUDA/HIP static-mixed, and CUDA/HIP adaptive cases. Nonuniform seeded CPU/GPU fixtures compare state decisions,
per-term energies, atom/coarse field checksums, transition counts, and final
restart trajectories. A GPU `EWALD3D_FFT` mixed-state fixture proves that the
production convolution contributes a nonzero adaptive dipole energy and FFT
timing. CPU assertions compare active updates and prove that mixed and coarse
modes reduce short-range/integration work. GPU cases assert compact active
counts and are skipped only when the compiled backend reports that no device
is present. Additional executable cases cover `Initmag=2`, CPU/GPU
atomistic-SD handoff, GPU-MC-to-GPU-CG handoff, deterministic `Initmag=8`,
and validated `M/H/Q/Y/Z/G` atomistic preparation. `Q/Y/Z` cases retain
inhomogeneous direction checksums, while an `X` negative fixture proves that
replica exchange is rejected before its runner starts.

The ordinary inputs in `examples/AdaptiveCoarseGraining` cover a static mixed
mask and adaptive MAX_ANGLE selection. They contain no test-only hooks or
manually staged arrays and are executed by the same e2e test.
