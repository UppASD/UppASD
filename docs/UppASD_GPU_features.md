---
marp: true
style: |
  section h1 {
    text-align: center;
    }
  .column {
    float: left;
    width: 50%;
    outline: 20px solid #FFF;
    border: 20px solid #AAF;
    background-color: #AAF;
    }
  .row:after {
    display: table;
    clear: both;
    }
---

<!-- paginate: true -->

# Inventory of UppASD functionality in Fortran source code

- Listing of the functionality implemented in the Fortran code. 
- Items sorted to have high, medium, or low priority for CUDA implementation.
- For a first iteration of the inventory, see

https://github.com/johanhellsvik/UppASD/blob/inventory/docs/FEATURES.md

---

# UppASD functionality
- [Critical temperature][11] determination
- [Magnon spectra][7] simulations
- [Micromagnetic exchange and DMI stiffness][6] calculations
- [Current][9], [field][8], and [temperature-driven][5] dynamics
- [LLG][12], Metropolis, and Heat-bath algorithms
- Generalized Hamiltonian
- Support for [arbitrary crystal structures][14] and [random alloys][13]
- [Spin-transfer][9] and spin-orbit torques
- [Temperature gradients][4]
- [Quantum statistics][3]
- [Multiscale (muASD) functionalities][16]
- [Spin-lattice dynamics (SLD)][17]

---

# Functionality currently available on GPU
- Spin dynamics measurement phase
- Tensorial exchange and single site anisotropy Hamiltonians
- The Depondt solver for the SLLG
- Copying of state variables back to CPU memory for measurements

---

# To GPU, high priority - short term

- Measurements of magnetization m and energy E
- Initial phase
- Metropolis Monte Carlo
- Heat-bath Monte Carlo
- J-tensor support to Metropolis and Heat-bath MC

---

# To GPU, high priority - longer term

- Spin-lattice dynamics (SLD)
- Spin-transfer torque (STT)
- Symplectic solvers
- multiscale (muASD) functionalities.
- Quantum Heath Bath Method
- Non-local damping
- Calculation of dipole-dipole interaction using the FFT method. 
    - Implementations for both the FFT MKL and FFTW libraries
    - Implementations for HeFFTe libraries for heterogeneous architectures 
    - Capabilities to treat open boundaries conditions.

---

# To GPU, medium priority I
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

---

# To GPU, medium priority II
- Frequency convolution for S(q,w) (earlier only for DOS). 
- Support for averaging of S(q,w)
- Calculation of the scalar S(q,w) scattering intensity.
- Cross-correlations: spin-displacement and displacement-velocity.
- [Dipole-dipole interaction via the macro-cell model](https://iopscience.iop.org/article/10.1088/0953-8984/28/6/066001/meta). 
    - The capability of reading and writing the dipole-dipole tensor for both the macro dipole and brute force dipole interactions, which could be used to improve the speed of some calculations.
    - The capacity to write the dipolar fields in certain cells  produced by all the other cells. This should allow for a multi-scale approach, where the region of interest is studied with the full dipole-dipole interaction plus the far field contribution.
---

# To GPU, medium priority III


- Wang-Landau implementation introduced. 
- Longitudinal spin fluctuations (LSF)
- Random alloy
