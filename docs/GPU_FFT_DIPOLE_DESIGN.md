# GPU dipole solver: particle--mesh Ewald design

## Status and authority

This is the CV6 production design.  It is a new GPU macrocell
particle--particle/particle--mesh Ewald (P3M/PME) solver for magnetic dipoles;
it is **not** a port of `dipolemanager.f90`, `fftdipole_*`, or `ewaldmom.f90`.
Those CPU paths are useful diagnostic comparisons, but neither their tensor
convention nor their numerical result is an acceptance oracle.

The correctness authority is (in order): an independently implemented
high-accuracy Ewald evaluator, analytic limiting cases, the stated boundary
condition, and published Ewald/PME theory.  CPU captures are regression data
only and must be labelled as such.

### Working notes

- **2026-07-22 — CV6.0 started.** Added
  `tests/dipole_validation/periodic_ewald_reference.py`: a compact fp64,
  direct real-plus-reciprocal Ewald evaluator for the block-one macrocell
  baseline.  It omits `k=0` (tin-foil convention), applies the point-dipole
  self correction, and is intentionally slow and independent of UppASD/GPU
  code.  Its alpha-invariance and lattice-translation self-checks live in
  `test_periodic_ewald_reference.py`.  The checks pass (alpha invariance and
  lattice translation invariance) and now run first from
  `tests/dipole_validation/run.py`.
- **2026-07-22 — independent-oracle regression values.** Added fixed fp64
  values for a cubic one-macrocell case and a skew two-macrocell case.  The
  latter exercises the reciprocal cell matrix rather than an orthogonal-box
  shortcut.  The same test also evaluates four independent ensembles, defining
  the required batched-PME comparison pattern.  These are dimensionless
  dipole-tensor values; the UppASD `mu0*muB/(4*pi*alat^3)` prefactor is applied
  only at the production backend boundary.  All periodic-oracle checks and the
  five existing finite analytic cases pass in the current workspace.

There is exactly one source representation: **macrocells**.  With macro block
size one, each atom is one macrocell and this same solver is atomically
resolved.  Larger blocks are a controlled coarse-grained approximation; there
is no separate atomistic dipole backend.  The macrocell map, membership and
centres are exported from Fortran, never reconstructed from floating-point
coordinates in C++.  Open finite samples and periodic slabs are separate
explicitly selected modes.  A mode is never inferred from the exchange FFT
backend or silently substituted for another mode.

## Physical contracts

### 3D periodic PME -- first production mode

The default contract is a 3D-periodic lattice of macrocell magnetic moments
with the **conducting (tin-foil) exterior convention**: the reciprocal `k=0`
mode is omitted.  When the macro block is one this is the point-dipole Ewald
problem.  The selected exterior dielectric/surface term is part of the
Hamiltonian, not a numerical detail.  A vacuum/spherical-sample surface term
may be added later as a separately named option after an analytic derivation
and tests; it must never be obtained by copying the current CPU implementation.

The energy is defined once as the Ewald-split dipole Hamiltonian, including
real-space, reciprocal-space, assignment/deconvolution, self, and any explicit
surface terms.  The field is the negative derivative of that same discrete
energy with respect to each magnetic moment.  This single definition prevents
the common failure in which the reported energy and the integrated field have
different interpolation or self conventions.

### Open finite sample -- second mode

This mode uses the macrocell source grid and a zero-padded linear convolution
on `(2*G1-1, 2*G2-1, 2*G3-1)`.  With block size one, `G` has atomic resolution
and its reference is the analytic finite point-dipole sum.  It has no Ewald
surface term and is not an approximation to the periodic PME mode.

### 2D-periodic slab -- later mode

Slab Ewald/PME needs a separately derived two-dimensional periodic Green
function or a validated electrostatic-layer correction.  Adding vacuum padding
to 3D PME is not sufficient.  This mode remains unavailable until its `k=0`,
surface, and layer-correction convention have analytic tests.

## GPU algorithm: 3D P3M/PME

1. **Macrocell geometry.** Form each macro moment
   `M_c=sum(i in c) emomM_i`, its centre and finite-cell metadata from the
   exported map.  Form the full cell matrix `H=[A B C]`, its reciprocal matrix,
   volume and fractional macrocell coordinates on the host in fp64.  Grid
   dimensions are chosen independently along the three fractional axes; skew
   cells are supported through `H`, not by treating Cartesian coordinates as
   orthogonal.  Block size one supplies one source per atom.
2. **Near field.** Build a GPU cell/Verlet list for the screened real-space
   dipole tensor between macrocells.  Evaluate the analytic Ewald real-space
   tensor within a cutoff, excluding only the exact self pair.  This is the
   particle--particle part and is independent of the exchange neighbour list.
3. **Mesh deposition.** Spread each vector moment to a periodic mesh with a
   cardinal B-spline of configurable order (initially 4 or 6).  The exact same
   weights, in adjoint order, interpolate the mesh field back to macrocells,
   then scatter the cell-constant field to member atoms.  This adjointness is
   asserted in a standalone test.
4. **Reciprocal field.** Use batched R2C/C2R transforms of the three mesh
   components.  Apply the dipolar Ewald reciprocal Green tensor and the PME
   influence/deconvolution function at every nonzero reciprocal mode.  Its
   formula is derived in the implementation note before coding, then tested
   against the independent Ewald evaluator; do not infer a sign or transpose
   from the old Fortran arrays.
5. **Corrections.** Apply the analytic self field/energy and the declared
   `k=0` surface convention.  Coarse macrocells additionally need an explicit
   finite-cell form factor and self-demagnetising correction; for block size
   one the point self correction is recovered.  The reciprocal mesh and
   real-space contributions are accumulated into `beff`; no term overwrites
   exchange, anisotropy or external fields.  Measurement reduces
   `-1/2 sum_c M_c.B_c` once per macrocell and ensemble.

All persistent grids, half spectra, B-spline metadata, neighbour-list data and
FFT workspace are accounted for before allocation.  Plans use the work stream,
are batched over ensembles, and are destroyed before their tensors.

## Macrocell-resolution rule

Macrocell size is a numerical resolution control of this one solver.  The
block-one result is the acceptance baseline for the macro formulation; coarse
results must converge towards it as block size is reduced.  A coarse cell must
not be treated as a point dipole without declaring that approximation: its
finite-cell form factor and self term define the approximation.  CPU
`do_dip=2` agreement is informative but is never a proof of correctness.

## Inputs, ownership, and staging

Introduce an explicit GPU dipole selection rather than overloading
`do_dip`:

```
gpu_dipole_mode       OFF | PME3D | OPEN_FFT | SLAB_PME
gpu_dipole_surface    TINFOIL | VACUUM_SPHERE   # PME3D only
gpu_dipole_alpha      auto | positive value
gpu_dipole_rcut       auto | positive value
gpu_dipole_mesh       auto | n1 n2 n3
gpu_dipole_spline     4 | 6 | ...
```

The Fortran bridge exports the selected mode and numerical parameters along
with `alat`, `H`, macrocell map/membership/counts/centres and magnetic moments.
It rejects invalid combinations before allocation.  Once selected, the GPU
backend owns both the dipole field and dipole energy for the GPU step; no CPU
dipole calculation may also be applied.

## Validation gates

1. **Independent finite oracle.** Retain `check_fields.py` for isolated
   block-one macrocells, axial/transverse pairs, and thin films.  It establishes
   units, tensor sign, factor and `-1/2 sum mu.B` energy independently of
   UppASD.
2. **Independent periodic oracle.** Add a small fp64 Python/C++ reference
   evaluator that directly sums real and reciprocal Ewald terms with increasing
   cutoffs until converged.  It must not call UppASD routines or reuse GPU
   tensor code.  Cover block-one macrocells in cubic and triclinic cells,
   M=1/4, all moment orientations, and the declared tin-foil `k=0` convention.
3. **PME convergence.** At fixed physical system, sweep alpha, real cutoff,
   mesh spacing and spline order.  Field and energy must converge to the
   independent periodic oracle at a documented target (fp64 target `1e-10` on
   small systems, with a separately measured mesh-error budget for production
   grids).
4. **Analytic/literature checks.** Check the isolated-pair result, uniform
   thin-film easy-plane response for the open mode, and the documented surface
   response for the selected 3D Ewald exterior.  Cite and retain the numerical
   parameters/results in `docs/`.
5. **Macrocell convergence and engineering.** Sweep block size from one to
   coarse blocks, including incomplete edges, and document field/energy error
   against the block-one independent-Ewald reference.  Also exercise
   non-32-aligned N, multi-basis/skew cells, M=4, coexistence with sparse and
   exchange FFT paths, CUDA sanitizer and HIP execution.  CPU `do_dip=2/3`
   captures may be reported beside these results but cannot relax any
   analytical or independent-Ewald acceptance criterion.

## Literature basis

- U. Essmann *et al.*, “A smooth particle mesh Ewald method,” *J. Chem.
  Phys.* **103**, 8577 (1995), DOI
  [10.1063/1.470117](https://doi.org/10.1063/1.470117): B-spline PME,
  reciprocal mesh treatment and `N log N` scaling.
- S. W. de Leeuw, J. W. Perram and E. R. Smith, “Simulation of electrostatic
  systems in periodic boundary conditions. I,” *Proc. R. Soc. A* **373**, 27
  (1980), DOI [10.1098/rspa.1980.0135](https://doi.org/10.1098/rspa.1980.0135):
  dipolar periodic sums and their boundary/surface dependence.
- M. C\u00e9rda, V. Ballenegger and C. Holm, “Particle-particle particle-mesh
  method for dipolar interactions,” *J. Chem. Phys.* **135**, 184110 (2011),
  DOI [10.1063/1.3652921](https://doi.org/10.1063/1.3652921): dipolar P3M
  mesh-error analysis and implementation guidance.
