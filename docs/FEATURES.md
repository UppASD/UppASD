# Features
An inventory of functionality and features currently contained in UppASD.

***Applications:***
- [Critical temperature][11] determination
- [Magnon spectra][7] simulations
- [Micromagnetic exchange and DMI stiffness][6] calculations
- [Current][9], [field][8], and [temperature-driven][5] dynamics

***Features:***
- [LLG][12], Metropolis, and Heat-bath algorithms
- Generalized Hamiltonian including
  - bilinear and [biquadratic Heisenberg][10] exchange
  - [Dzyaloshinskii-Moriya][15] exchange
  - magnetocrystalline anisotropies.
- Support for [arbitrary crystal structures][14] and [random alloys][13]
- [Spin-transfer][9] and spin-orbit torques
- [Temperature gradients][4]
- [Quantum statistics][3]
- [Multiscale (muASD) functionalities][16]
- [Spin-lattice dynamics (SLD)][17]

### Functionality available on GPU
- Spin dynamics measurement phase
- Tensorial exchange and single site anisotropy Hamiltonians
- The Depondt solver for the SLLG
- Copying of state variables back to CPU memory for measurements

### To GPU, high priority

- Measurements of magnetization m and energy E
- Initial phase
- Metropolis Monte Carlo
- Heat-bath Monte Carlo
- J-tensor support to Metropolis and Heat-bath MC
- Calculation of dipole-dipole interaction using the FFT method.
  - Implementations for both the FFT MKL and FFTW libraries
  - Implementations for HeFFTe libraries for heterogeneous architectures
  - Added capabilities to treat open boundaries conditions.

### To GPU, medium priority

- Spin-lattice dynamics (SLD)
- Spin-transfer torque (STT)
- Symplectic solvers
- multiscale (muASD) functionalities.
- Non-local damping
- Introduced **major new feature**, 3 Temperature model:
  - Added electron-only mode for 3TM.
  - Added variable specific heats for 3TM
  - Added 3TM example: `./SLDexamples/bccFe3TM`.
- Weights to cumulants to improve convergence.
- New function for real-space and time correlations.
- Triangulation routine for skyrmion number calculation.
- JSON output for cumulant data (for aiida-uppasd parsing).
- Minimal driver routines for ASD and MC for external use.
- Velocity correlations.
- Frequency convolution for S(q,w) (earlier only for DOS).
- Support for averaging of S(q,w)
- Calculation of the scalar S(q,w) scattering intensity.
- Cross-correlations: spin-displacement and displacement-velocity.
- Quantum Heath Bath Method
- [Dipole-dipole interaction via the macro-cell model](https://iopscience.iop.org/article/10.1088/0953-8984/28/6/066001/meta).
  - The capability of reading and writing the dipole-dipole tensor for both the macro dipole and brute force dipole interactions, which could be used to improve the speed of some calculations.
    - The capacity to write the dipolar fields in certain cells  produced by all the other cells. This should allow for a multi-scale approach, where the region of interest is studied with the full dipole-dipole interaction plus the far field contribution.
- Wang-Landau implementation introduced.

### To GPU, lower priority

- Temperature gradients
- Local field pulses
- Local temperature pulses
- Variant of induced moment treatment implemented.
- Metatype concept for projected correlations.
- ASD functionality to ip_mode SX.
- Topological center of mass as an observable
- Local angular momentum measurable for atoms (local frame).
- Temperature control to LSWT correlation functions.
- Polesya's original formulation for induced moments.
- Four-site ring exchange routines
- Rudimentary parallel tempering algorithm
- `do_bls_local_axis B` for correlations perpendicular to the local B-field.
- Capability for induced moment calculation is SD mode.
  - The calculation is based on the approach by Polseya and it treats the induced magnetic moments motion as being slaved by the motion of the nearest neighbour fixed moments.
  - Introduced a new array which contains whether a moment is induced or not while taking into account the chemical nature of the system.
  - Currently this type of dynamics is only taken into account in the Depondt solver.
  - Added the capability of having fixed moments and induced moments at the same time.
  - The selective Depondt solver will only update spins that 1) Are not induced 2) Are not fixed.
  - The magnitude and direction of the induced moments depends on the neighbouring fixed moments. This might cause problems at high temperatures.
- Added the generalized form for the SOT.
  - There are parameters controlling the damping-like torque and the field-like torque.
  - The sign of the parameters might still be in question.
  - There is also possibility for site dependent spin polarizations.

### Functionality that can stay on CPU, or global for CPU and GPU e.g. example input files

- Support for (C)ycloidal and (H)elical spirals in qminimizer.
- Support for angular frequency printing for S(q,w) outputs.
- Gradient of temperature rescaling parameter required for Cv in QHB.
- Spin spiral and skyrmion examples.
- Python interface for UppASD
- More general tensor J support to sparse routines.
- Temperature dependent specific heat treatment for SLD.
- CMake support.
- QHB example, reworked fccCo example.
- PGI/NVIDIA fortran to `make probe`.
- ELK/RSPt compatible input/output.
- `Maptype 3` for RPS Jij input files.
- Scripts that Parse SPRKKR data into ASD format :
  -The script also allows one to plot the Jij's as a function of distance, as well as some optional files.
  - If one has a `.sys` file and a `_SCF.out` the program will also generate momfiles and posfiles.
  - It is possible to write an exclusion list where one can write the atomic symbols of the atoms that one does not want to plot and which should be excluded from the ASD files
- Python scripts for BLS and S(q,w) plotting.
- Reading of site specific symmetry elements for inpsd.dat flag value 'sym 5'.
- Support for simultaneous C(q) and S(q,w) with `do_sc Y`.
- Custom type for Hamiltonian input data.
- Magnon DOS printout for NC-AMS.
- Introduced two-time autocorrelation into `mc_driver`.
- Edwards-Anderson model capability.
- 1q and 3q spin-spiral minimizer.
- Fixes for tensor exchange simulations and LSWT.
- `qpoints=='G'` to read old q-points on the `qpoints.out` format.
- `postQ.py` for plotting AMS according to `preQ.py` structure.
- SLD regressions tests.
- Support for `do_reduced Y` to NC-AMS.
- Makefile for win64 (cross-platform build).
- Support for sparse linear algebra for the evaluation of the effective field. Only for Jij and Dij and currently hardcoded to MKL.
- eta-functionality for DM-stiffness.
- Enabled the spin-ice loop algorithm for the initial phase.
- Implementation of AMS in random alloys.
- Chemically resolved AMS for random systems added.
- New S(q,w) sampling option.
- qpoints weights in AMS-MDOS calc if symmetry is used.
- Random stiffy, i.e. supercell averaged for random alloys.
- Accounting for PBC in `printhamiltonian`. Works for orthogonal systems, but needs verification for general lattices.
- Introduction of a timing routine for the dipole-dipole calculations. Added wrapper routines to handle the calculation of the dipole-dipole field and energy for the different methods implemented.
- `Third_party` folder that includes software from other sources. Included the OVF reader capability developed by Gideon Muller. Started to develop the I/O for the ovf format.
- Three-site (mi * (mj x mk)) interaction.
- Functionality of removing the parallel component of the effective field with keyword `perp Y`
- STT capability for the midpoint solver.
- Local Spin Fluctuations (LSF) introduced.
- Made AMS working for reduced Hamiltonian.
- Average energy in output file.
- New coordinate wrapping routines for stiffy and AMS. Works ok for ra_AMS.

---
(C) 2008-2024 [UppASD group][2]

[1]:https://global.oup.com/academic/product/atomistic-spin-dynamics-9780198788669
[2]:http://www.physics.uu.se/research/materials-theory/ongoing-research/uppasd/
[3]:https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013802
[4]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.014434
[5]:https://www.nature.com/articles/ncomms12430
[6]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.214424
[7]:http://iopscience.iop.org/article/10.1088/0953-8984/27/24/243202/meta
[8]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.224401
[9]:https://www.nature.com/articles/srep25685
[10]:https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.127204
[11]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.214439
[12]:http://iopscience.iop.org/article/10.1088/0953-8984/20/31/315203
[13]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.214410
[14]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.094411
[15]:https://www.nature.com/articles/ncomms5815
[16]:https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.013092
[17]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.104302
[logo]:https://github.com/UppASD/UppASD/blob/master/docs/uppasd_rot.png
