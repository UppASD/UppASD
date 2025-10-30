# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [6.1.0] - 2025-09-11

### Added

- Added **major new feature**, Chern number calculations for 3D systems and thermal magnon conductivity. [[Nastaran Salehi](https://github.com/Nastaran-Salehi)]
  - Added Chern number calculation subroutine for topological magnon properties
  - Added Kagome Chern number example: `./examples/SpinWaves/Kagome_ChernNumber`
- Added **major enhancement**, comprehensive GUI improvements for ASD_GUI. [[Anders Bergman](anders.bergman@physics.uu.se)]
  - Added PyPI installation support with `pip install asd_gui` command
  - Added interactive simulation capabilities and real-time visualization
  - Added comprehensive input file creation tools with structure templates
  - Added support for importing .cif, SPRKKR, and RS-LMTO file formats
  - Added J-file creation using neighbor vector search algorithms
  - Added HDR image-based lighting and advanced visualization controls
  - Added FXAA anti-aliasing and SSAO (Screen Space Ambient Occlusion) support
  - Added settings persistence and programmatic GUI updates
- Added **major enhancement**, CUDA/GPU acceleration improvements. [[Mariia Mohylna](mmogylnaya@gmail.com), [Arkadijs Slobodkins](aslobodkins@smu.edu)]
  - Added tensorial exchange (J-tensor) support to GPU calculations
  - Added uniaxial and cubic anisotropy implementation for GPU
  - Added reduced Hamiltonian option support for CUDA
  - Enhanced GPU memory management and error handling
  - Improved CUDA build system integration with CMake
- Added exponential temperature profile functionality for time-dependent simulations. [[Anders Bergman](anders.bergman@physics.uu.se)]
  - Added exponential temperature example: `./examples/SpecialFeatures/ExpTemperature`
- Added local autocorrelation measurements with macro cell support. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added modified temperature model for "3TM without lattice" (2-temperature model). [[Anders Bergman](anders.bergman@physics.uu.se)]
  - Added 2TM example: `./examples/SpecialFeatures/bccFe2TM`
- Added enhanced multiscale (muASD) functionalities with improved examples. [[Nastaran Salehi](https://github.com/Nastaran-Salehi), [Manuel Pereiro](manuel.pereiro@physics.uu.se)]
  - Added skyrmion-stress multiscale example
  - Added spin-states multiscale example  
  - Added wall-dislocation multiscale example
  - Improved STT calculations in multiscale simulations
- Added additional algorithms for randomizing exchange interactions (Jij) and DMI. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added scale factor support for local STT current densities. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added fallback routine for bosonic diagonalization in AMS calculations. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added relaxation functionality to spin-spiral minimizers. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for non-cubic macrocells in coarse-graining procedures. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added cell rotation functionality to q-space minimizer (qminimizer). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added cmake build support for coverage testing and improved build profiles:
  - Added AOCC (AMD Optimizing C/C++ Compiler) build profile
  - Added Cray compiler wrapper support (crayftn-ftn)
  - Added gfortran-ftn build profile for Cray environments
- Added comprehensive example scripts for running test suites:
  - Added scripts for SimpleSystems, PhaseDiagrams, Mappings examples
  - Added scripts for SpecialFeatures, SpinLattice, SpinWaves examples
  - Added script for running complete example suite
- Added GitHub Actions for automated building and testing. [[Anders Bergman](anders.bergman@physics.uu.se)]
  - Added wheel-building workflows for Python package distribution
  - Added Windows installer automation and cross-platform build support
  - Added status badges and CI/CD improvements

### Changed

- Updated CUDA implementation with modern C++ standards and improved memory management [[Arkadijs Slobodkins](aslobodkins@smu.edu)]
- Restructured GUI codebase with better modularization and UI organization [[Anders Bergman](anders.bergman@physics.uu.se)],[[Jonathan Chico](jonathanpchico@gmail.com)]
- Enhanced build system with better compiler detection and optimization flags [[Anders Bergman](anders.bergman@physics.uu.se)]
- Improved example organization and documentation structure
- Updated Python scripts for better Python 3 compatibility and modern dependencies

### Fixed

- Fixed CUDA compilation issues for modern GPU architectures and compute capabilities
- Fixed GUI stability issues and improved error handling in visualization components
- Fixed multiscale calculation accuracy in various solver implementations
- Fixed memory allocation issues in CUDA temperature field calculations
- Fixed cross-platform build compatibility issues for Windows and macOS

## [6.0.1] - 2022-10-08
 
### Added
- Added **major new feature**, multiscale (muASD) functionalities. [[Nikolaos Ntallis](nikos.ntallis@physics.uu.se), [Manuel Pereiro](manuel.pereiro@physics.uu.se)]
- Introduced **major new feature**, 3 Temperature model [[Anastasiia A. Pervishko](anastasiia.pervishko@physics.uu.se), [Anders Bergman](anders.bergman@physics.uu.se)]:
    - Added electron-only mode for 3TM. [[Anders Bergman](anders.bergman@physics.uu.se)]
    - Added variable specific heats for 3TM. [[Anders Bergman](anders.bergman@physics.uu.se)]
    - Added 3TM example: `./SLDexamples/bccFe3TM`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added weights to cumulants to improve convergence. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added new function for real-space and time correlations. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added triangulation routine for skyrmion number calculation. `skyno T`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added JSON output for cumulant data (for aiida-uppasd parsing). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added ASD functionality to ip_mode SX. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added minimal driver routines for ASD and MC for external use. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for (C)ycloidal and (H)elical spirals in qminimizer. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added topological center of mass as a measurement. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added symmetric anisotropic exchange as individual interaction. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added pypi installation support for ASD_GUI. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for angular frequency printing for S(q,w) outputs. [[Anders Bergman](anders.bergman@physics.uu.se)] 
- Added test functionality to Cmake compilation (ctest).  [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added Quantum Heath Bath Method [[Anders Bergman](anders.bergman@physics.uu.se), [Lars Bergqvist](lbergqv@kth.se)]
    - Added gradient of temperature rescaling parameter required for Cv in QHB. [[Lars Bergqvist](lbergqv@kth.se)]
- [Addition of a changelog](https://gitlab.com/UppASD/UppASD/-/issues/34). [[Jonathan Chico](jonathanpchico@gmail.com), [Anders Bergman](anders.bergman@physics.uu.se)]
- Added more general tensor J support to sparse routines. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added spin spiral and skyrmion examples. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added python interface for UppASD [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added temperature dependent specific heat treatment for SLD. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added cross-correlations: spin-displacement and displacement-velocity. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added metatype concept for projected correlations. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added velocity correlations. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added python scripts for BLS and S(q,w) plotting. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Reading of site specific symmetry elements for inpsd.dat flag value 'sym 5'. [[Johan Hellsvik](hellsvik@kth.se)]
- Added support for simultaneous C(q) and S(q,w) with `do_sc Y`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Cleanup and added new S(q,w) printout feature. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added custom type for Hamiltonian input data. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added local angular momentum measurable for atoms (local frame). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added temperature control to LSWT correlation functions. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added Polesya's original formulation for induced moments. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added Four-site ring exchange routines.  [[Natalya S. Fedorova](natalya.fedorova@list.lu), [Anders Bergman](anders.bergman@physics.uu.se), [Johan Hellsvik](hellsvik@kth.se)]
- Added Magnon DOS printout for NC-AMS. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Introduced 3 Temperature model [[Anastasiia A. Pervishko](anastasiia.pervishko@physics.uu.se), [Anders Bergman](anders.bergman@physics.uu.se)]:
    - Added electron-only mode for 3TM. [[Anders Bergman](anders.bergman@physics.uu.se)]
    - Added 3TM example: `./SLDexamples/bccFe3TM`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Introduced two-time autocorrelation into `mc_driver`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added Edwards-Anderson model capability. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Implemented 1q and 3q spin-spiral minimizer. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added J-tensor support to Metropolis MC. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixes for tensor exchange simulations and LSWT. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added frequency convolution for S(q,w) (earlier only for DOS). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for averaging of S(q,w) and updated regression checking routines/references. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added cmake support. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added QHB example, reworked fccCo example. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added PGI/NVIDIA fortran to `make probe`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Introduced ELK/RSPt compatible input/output. [[Johan Hellsvik](hellsvik@kth.se)]
    - `Maptype 3` for RPS Jij input files. [[Johan Hellsvik](hellsvik@kth.se)]
- Added a series of scripts that Parse SPRKKR data into ASD format [[Jonathan Chico](jonathanpchico@gmail.com)]:
    -The script also allows one to plot the Jij's as a function of distance, as well as some optional files.
    - If one has a `.sys` file and a `_SCF.out` the program will also generate momfiles and posfiles.
    - It is possible to write an exclusion list where one can write the atomic symbols of the atoms that one does not want to plot and which should be excluded from the ASD files
- Variant of induced moment treatment implemented. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Calculation of the scalar S(q,w) scattering intensity. [[Johan Hellsvik](hellsvik@kth.se)]
- Added `qpoints=='G'` to read old q-points on the `qpoints.out` format. [[Johan Hellsvik](hellsvik@kth.se)]
- Added `postQ.py` for plotting AMS according to `preQ.py` structure. [[Anders Bergman](anders.bergman@physics.uu.se)]
- SLD regressions tests added.  [[Johan Hellsvik](hellsvik@kth.se)]
- Rudimentary Parallel Tempering algorithm added. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added `do_bls_local_axis B` for correlations perpendicular to the local B-field. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for `do_reduced Y` to NC-AMS. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added makefile for win64 (cross-platform build). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added support for sparse linear algebra for the evaluation of the effective field. Only for Jij and Dij and currently hardcoded to MKL. [[Anders Bergman](anders.bergman@physics.uu.se)]
- [Added dipole-dipole interaction via the macro-cell model](https://iopscience.iop.org/article/10.1088/0953-8984/28/6/066001/meta). [[Jonathan Chico](jonathanpchico@gmail.com)]
    - Added the capability of reading and writing the dipole-dipole tensor for both the macro dipole and brute force dipole interactions, which could be used to improve the speed of some calculations.
    - Added the capacity to write the dipolar fields in certain cells  produced by all the other cells. This should allow for a multi-scale approach, where the region of interest is studied with the full dipole-dipole interaction plus the far field contribution.
- Introduced the calculation of the dipole-dipole interaction using the FFT method. [[Jonathan Chico](jonathanpchico@gmail.com)]
    - Implementations for both the FFT MKL and FFTW libraries.
    - Added capabilities to treat open boundaries conditions.
- Added capability for induced moment calculation is SD mode. [[Jonathan Chico](jonathanpchico@gmail.com)]
    - The calculation is based on the approach by Polseya and it treats the induced magnetic moments motion as being slaved by the motion of the nearest neighbour fixed moments.
    - Introduced a new array which contains whether a moment is induced or not while taking into account the chemical nature of the system.
    - Currently this type of dynamics is only taken into account in the Depondt solver.
    - Added the capability of having fixed moments and induced moments at the same time.
    - The selective Depondt solver will only update spins that 1) Are not induced 2) Are not fixed.
    - The magnitude and direction of the induced moments depends on the neighbouring fixed moments. This might cause problems at high temperatures.
- Added the generalized form for the SOT. [[Jonathan Chico](jonathanpchico@gmail.com)]
    - There are parameters controlling the damping-like torque and the field-like torque. 
    - The sign of the parameters might still be in question.
    - There is also possibility for site dependent spin polarizations.
- Added eta-functionality for DM-stiffness. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Enabled the spin-ice loop algorithm for the initial phase. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Implementation of AMS in random alloys. [[Lars Bergqvist](lbergqv@kth.se)]
    - Chemically resolved AMS for random systems added.
- Added new S(q,w) sampling option. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added qpoints weights in AMS-MDOS calc if symmetry is used. [[Lars Bergqvist](lbergqv@kth.se)]
- Added random stiffy, i.e. supercell averaged for random alloys. [[Lars Bergqvist](lbergqv@kth.se)]
- Accounting for PBC in `printhamiltonian`. Works for orthogonal systems, but needs verification for general lattices. [[Johan Hellsvik](hellsvik@kth.se)]
- Introduction of a timing routine for the dipole-dipole calculations. Added wrapper routines to handle the calculation of the dipole-dipole field and energy for the different methods implemented. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Addition of the `Third_party` folder to include software from other sources. Included the OVF reader capability developed by Gideon Muller. Started to develop the I/O for the ovf format. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Implemented three-site (mi * (mj x mk)) interaction. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added possibility of removing parallel components of B_eff. [[Anders Bergman](anders.bergman@physics.uu.se)]
    - Added the possibility of removing the parallel component of the effective field with keyword `perp Y`
- Added STT capability for the midpoint solver. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Wang-Landau implementation introduced. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Local Spin Fluctuations (LSF) introduced. [[Lars Bergqvist](lbergqv@kth.se)]
- Made AMS working for reduced Hamiltonian. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Added average energy in output file. [[Lars Bergqvist](lbergqv@kth.se)]
- Added new coordinate wrapping routines for stiffy and AMS. Works ok for ra_AMS. [[Anders Bergman](anders.bergman@physics.uu.se)]

### Changed
- Improved m=z handling for MC HB. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Improved support for non-tensor interactions (anisotropy+SA) in nc-AMS. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Changed logfile format from YAML to json. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Cmake now builds library and executable. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Unified I/O for magnetic configurations. Now **all** magnetic configurations (moments, restartfile) have the same structure. [[Jonathan Chico](jonathanpchico@gmail.com)]
- [Renaming the examples and tests folders](https://gitlab.com/UppASD/UppASD/-/issues/32) so that they follow the current standards. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Removed symmetrization of K in NC-LSWT routines. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Reimplemented G(k,t). Accessible with `do_sc T`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Revision of SLDexamples. [[Johan Hellsvik](hellsvik@kth.se)]
- Made the memory profiling file `meminfo` optional (`do_meminfo 0/1`). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Improved (uniaxial) MAE support for LSWT. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Changed Verlet method for zero lattice damping. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Changes to S(q,w) tensor calculations. Changed real transform in S(q) to complex transform. Updated regression references accordingly. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Reactivation of point group symmetry in `setup_nm_nelem`. [[Johan Hellsvik](hellsvik@kth.se)]
- Added conditional compilation of OVF support. Off by default. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Changed back so that `posfiletype` only affects `posfile`, not `jfile`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Updated default CUDA compute capability to 3.0. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Increased sanity checks for `inputhandler`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- `read_elk_geometry` sets `do_ams='N'` if system is noncollinear. [[Johan Hellsvik](hellsvik@kth.se)]
- Changes to `read_llphonopydata`. [[Johan Hellsvik](hellsvik@kth.se)]
- Updated the convention for direct q-vectors. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Change of `qpoints==H` functionality.  [[Anders Bergman](anders.bergman@physics.uu.se)]
- Auxiliary STT fields now allocated according to different values of the `ipSDEalgh` flag. [[Johan Hellsvik](hellsvik@kth.se)]
- For `qpoints=='C'`, writes q points in direct coordinates. [[Johan Hellsvik](hellsvik@kth.se)]
- Added wrapper routine for the calculation of the geometry of the cluster. [[Jonathan Chico](jonathanpchico@gmail.com)]
    -This routine should allow for the capacity of having atoms in the cluster that are not in the magnetic layers, without the need of inserting  vacuum  layers. Thus this requires the modification of the number of atoms in the system.
    -Fixed an issue which would cause deposited atoms in the cluster method to not properly work with the midpoint solver.
- Change to have mRyd as the fundamental unit for energy in LD and SLD subroutines. [[Johan Hellsvik](hellsvik@kth.se)]
- Changed the `abs(w)` usage in `ams()`. Removed `taniso==3` (CMA) from the obsolete file `applyhamiltonian.f90`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Re-instated the first generation temperature gradient routines. So far confirmed stable only for 1d gradients. Also updated the default `sc_local_axis_mix` to zero. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Updated `examples/SLDexamples/bccFeSLD66/` to use Jij with 6 lattice constants cutoff for r_ij. [[Johan Hellsvik](hellsvik@kth.se)]
- Format change from text to aml for `inp.simid.out`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Improvements and generalization to AF-DLM in the renormalization exchange method. [[Lars Bergqvist](lbergqv@kth.se)]
- Started using `effective_field()` for HeatBath MC. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Updated python-scripts to python3 compatibility. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Restored type resolved AMS using super-VCA of random alloys. [[Lars Bergqvist](lbergqv@kth.se)]
- Export of geometry to POSCAR format. [[Johan Hellsvik](hellsvik@kth.se)]
    - Added flag `do_prn_poscar`.
- Changed the format of the neighbour list to be better suited for large systems. [[Jonathan Chico](jonathanpchico@gmail.com)]
    - Changes to the dmdata file so that it now has the same structure as the structfile.
    - This file contains the bonding vectors, atom types, atom number and interaction strength in mRy
- Changed format for site-fields when using `locfield`. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed missing 1/(1+a^2) in the Heun_proper solver. Added STT terms to the Heun solver. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Refactoring of the GNEB routines [[Jonathan Chico](jonathanpchico@gmail.com)]:
    - Added a wrapper for the calculation of the Hessian in GNEB 
    - Added test suite for the GNEB functionality via the `make tests-gneb` command. This functionality will test the VPO minimizer, the calculation of eignevalues, the prefactors obtained from the Hessians as well as parameters comming from the energy barrier calculation.
    -Have also modified the output from the GNEB routines to be more UppASD and testing compliant, as well as modified the names of several of the output files so that they follow the UppASD convention. 
- Moved out step size in Gaussian rotation trial move in order to prepare for adaptive sampling. [[Lars Bergqvist](lbergqv@kth.se)]
- Adjusted cutoff energy in MDOS calc. [[Lars Bergqvist](lbergqv@kth.se)]
- Changed the way in which the energy is being calculated so that it is consistent with the effective field calculation. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Changed datatype for `do_reduced` to character (Y/N) Preprocessed sparse mkl-calls. Added `do_reduced Y` to the kagome test. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Introduced the Hamiltonian type. [[Anders Bergman](anders.bergman@physics.uu.se)]
### Fixed
- Fixed [issue with short simid strings](https://gitlab.com/UppASD/UppASD/-/issues/31). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed CUDA compilation for now. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed topological charge calculation. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed a bug in initializing initial and final states for GNEB calculations; made calculations of hessians more efficient. [[Pavel Bessarab](bessarab@hi.is)]
- Revision to `read_llphonopydata`. [[Johan Hellsvik](hellsvik@kth.se)]
- Fixed bounds in magnon dos routine. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed sqw for shorter simulation than sampling times. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Introduced lvec for lattice L. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Changed G(k) routine. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed issue with PGF90 compilers. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed q-point array limits in nc-s(q,w) subroutine. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Optimized SLD for diagonal Aijk cases (mml_diag T). [[Anders Bergman](anders.bergman@physics.uu.se)]
- Corrected the use of `dm_scale` for rescaling DM interactions. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed lattice temperature for SLD initial phase. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed for NC-ams for large supercells. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Fixed bug for reading of harmonic lattice potential on phonopy format. [[Johan Hellsvik](hellsvik@kth.se)]
- Inactivation of `timesym` and `invsym` in `get_fullnnlist_nelem`. [[Johan Hellsvik](hellsvik@kth.se)]
- Lattice measurements; correction of garlic NaN for T=0 K. [[Johan Hellsvik](hellsvik@kth.se)]
- Correction of bug in calculate_ams after bug being reported by Roberto Díaz Pérez. For systems without inversion symmetry, the elements of the A matrix can take complex values. The matrix diagonalization is now performed for the complex matrix, and not for its real part. [[Johan Hellsvik](hellsvik@kth.se)]
- Fixed special q-point cases for qminimizer. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Corrected **MAJOR** bug in `magnetizationinit` for random alloys and `initmag 1`. [[Lars Bergqvist](lbergqv@kth.se)]
- OMP fix in random_ams causing random segfaults. [[Lars Bergqvist](lbergqv@kth.se)]
- Corrected Tc MFA and RPA for random alloys. [[Lars Bergqvist](lbergqv@kth.se)]
- Fixed a bug with `the_path` subroutine. [[Pavel Bessarab](bessarab@hi.is)]
- Corrected Tc-MFA estimate for random alloys. [[Lars Bergqvist](lbergqv@kth.se)]
- Corrected printout of J(q) in AMS. [[Lars Bergqvist](lbergqv@kth.se)]
- Correction to the parameters to ensure that the dipole-dipole interaction has the proper magnitude. [[Jonathan Chico](jonathanpchico@gmail.com)]
- Corrected scaling of the lattice vectors for POSCAR export. [[Johan Hellsvik](hellsvik@kth.se)]
- Corrected MC implementation for chiral term. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Corrected possible bug in Heat bath. [[Lars Bergqvist](lbergqv@kth.se)]
- Fixed call to temperature gradient routines. [[Anders Bergman](anders.bergman@physics.uu.se)]
### Ongoing
- Temperature models for SLD. [[Anders Bergman](anders.bergman@physics.uu.se)]
- General four-spin interaction term. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Re-formulation of induced moment functionalities. [[Anders Bergman](anders.bergman@physics.uu.se)]
### Deprecated

### Removed
- Removed default resetting of u and p in sld_driver. [[Anders Bergman](anders.bergman@physics.uu.se)]
- Removed old (`UppASD < 5.0`) neighbour list output format.  [[Jonathan Chico](jonathanpchico@gmail.com)]
- Removed quasibls functionality and test. [[Anders Bergman](anders.bergman@physics.uu.se)]

## [5.0.0] - 2017-11-23
To check all the features available in `UppASD` please visit the `README.md` where they are presented. All features not listed here are also described in the documentation.
