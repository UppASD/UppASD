# Historical CV6 GPU FFT/PME dipole handoff notes

> Superseded for implementation decisions by
> `FFT-dipole_implementation_plan.md` and `GPU_FFT_DIPOLE_DESIGN.md`.
> In particular, this handoff's staged `PME3D` spelling is obsolete: the
> implemented selector is `EWALD3D_FFT`.  Retain this file only for its audit
> history; do not use its mode contract as current documentation.

## Goal and physical contract

CV6 is intended to add a production GPU dipole solver to UppASD. It is not a CUDA/HIP port of the legacy CPU `dipolemanager.f90`, FFT dipole, or Ewald implementations.

The selected first production mode is:

```text
gpu_dipole_mode = PME3D
gpu_dipole_surface = TINFOIL
```

Its physics is a 3D-periodic lattice of macrocell moments with conducting/tin-foil boundary conditions:

- The reciprocal `k=0` contribution is omitted.
- The field convention is:

```text
B = Hessian(1/r) M
```

- Energy is:

```text
E = -½ Σcell Mcell · Bcell
```

- UppASD’s physical prefactor belongs at the production boundary:

```text
mu0 * muB / (4*pi*alat^3)
```

The Python oracle works in dimensionless dipole-tensor units and deliberately omits that UppASD prefactor.

The intended Ewald decomposition is:

```text
Btotal = Breal + Breciprocal + Bself [+ Bsurface]
```

For tin-foil PME3D:

```text
Bsurface = 0
```

The point self term under the adopted convention is:

```text
Bself = + 4 alpha^3 / (3 sqrt(pi)) M
```

The screened real-space tensor is the Hessian of:

```text
erfc(alpha r) / r
```

Written as:

```text
Treal(r) = A r r^T + B I
```

with:

```text
A =
  3 erfc(alpha r) / r^5
+ 6 alpha exp(-alpha^2 r^2) / (sqrt(pi) r^4)
+ 4 alpha^3 exp(-alpha^2 r^2) / (sqrt(pi) r^2)

B =
- erfc(alpha r) / r^3
- 2 alpha exp(-alpha^2 r^2) / (sqrt(pi) r^2)
```

The reciprocal-space tensor in the independent oracle is:

```text
Treciprocal(k) =
- 4 pi / V
  * exp(-k^2 / (4 alpha^2))
  / k^2
  * k k^T
```

for every nonzero reciprocal lattice vector `k`.

The reciprocal cell matrix is:

```text
Brec = 2*pi*(H^-1)^T
```

where the full periodic supercell matrix is:

```text
H = [N1*C1, N2*C2, N3*C3]
```

not simply the primitive UppASD `C1/C2/C3` vectors.

## Macrocell representation

There is one intended source representation: macrocells.

For a regular periodic system:

```text
M(g,a) = sum of atom moments in coarse cell g, basis channel a
```

where:

- `g` is a coarse Bravais-grid cell;
- `a = 1..NA` is the basis channel;
- `NA=1` is the first validation slice;
- block size one gives one source per atom for `NA=1`;
- regular `NA>1` must preserve distinct basis channels.

The regular-grid spectral operator should ultimately be a block convolution:

```text
B(g,a) = Σg',b K(g-g', a, b) M(g',b)
```

with a `3×3×NA×NA` tensor block at each grid displacement.

Its FFT representation is:

```text
Bhat(k,a) = Σb Khat(k,a,b) Mhat(k,b)
```

The current GPU spectral batch convention is:

```text
field/moment batch =
component + 3*(basis_channel + NA*ensemble)
```

The current kernel batch convention is:

```text
kernel batch =
row + 3*(column + 3*(target_basis + NA*source_basis))
```

This convention must be preserved by future kernel construction and multi-basis validation.

## Independent correctness authority

The primary correctness authority is:

```text
tests/dipole_validation/periodic_ewald_reference.py
```

It is deliberately independent of UppASD and GPU code.

It computes a direct fp64 real-plus-reciprocal Ewald sum with:

- explicit real image ranges;
- explicit reciprocal image ranges;
- tin-foil `k=0` omission;
- analytic point self correction;
- field convention `B = Hessian(1/r) M`;
- energy `-½ Σ M·B`.

Important oracle checks already present:

- alpha invariance after adequate real/reciprocal convergence;
- lattice-translation invariance;
- fixed cubic one-macrocell value;
- fixed skew-cell two-macrocell values;
- independent `M=4` ensemble checks;
- finite analytic CPU dipole checks.

The oracle has already caught a crucial flaw in the current GPU reciprocal builder.

## Implemented code and wiring

### Fortran macrocell and bridge work

Important Fortran-side files:

- `source/CoarseGraining/macrocells.f90`
- `source/uppasd.f90`
- `source/chelper.f90`
- `source/Input/inputdata.f90`
- `source/Input/inputhandler.f90`

Implemented/staged:

- `create_pme_macrocell_layout` creates a separate basis-resolved PME layout.
- Legacy CPU `do_dip=2` macrocell data remains separate; CV6 should not change legacy CPU physics.
- The PME layout exports:
  - `pme_Num_macro`;
  - `pme_macro_grid`;
  - `pme_cell_index`;
  - `pme_macro_nlistsize`;
  - `pme_macro_center`;
  - `pme_macro_min_coord`;
  - `pme_macro_max_coord`.

Input fields added:

```text
gpu_dipole_mode     OFF | PME3D | OPEN_FFT | SLAB_PME
gpu_dipole_surface  TINFOIL | VACUUM_SPHERE
gpu_dipole_alpha    0 | positive
gpu_dipole_rcut     0 | positive
gpu_dipole_mesh     0 0 0 | n1 n2 n3
```

Current semantics:

```text
0 for alpha/rcut/mesh = future auto-selection
```

No auto-parameter policy is implemented yet.

Important bug fixed:

- Initial mode/surface enum values were local variables in `FortranData_Initiate`, so C++ held dangling pointers after that routine returned.
- They now live persistently in `InputData`.
- Default mode is reliably `OFF`.

This was required because regular CUDA simulations were incorrectly hitting:

```text
GPU dipole mode was requested, but its field/energy operator is not yet available
```

### C++/GPU bridge and storage

Important files:

- `source/gpu_files/fortranData.hpp`
- `source/gpu_files/fortranData.cpp`
- `source/gpu_files/gpuSimulation.cpp`
- `source/gpu_files/gpuHamiltonianCalculations.hpp`
- `source/gpu_files/gpuHamiltonianCalculations.cpp`
- `source/gpu_files/gpuDipoleConvolution.hpp`
- `source/gpu_files/gpuDipoleConvolution.cpp`
- `source/gpu_files/gpuFftWrapper.hpp`

`FortranData` now exposes pointers for:

```text
gpu_dipole_mode
gpu_dipole_surface
gpu_dipole_alpha
gpu_dipole_rcut
gpu_dipole_mesh
```

The bridge keeps GPU dipole selection independent from legacy `do_dip`.

Current production guard in `GpuHamiltonianCalculations::initiate`:

- Non-`OFF` GPU dipole mode fails fast.
- This is intentional until the field, energy, and validation path are complete.
- Legacy `do_dip=2` is still used only to stage macrocell data and current dormant buffers.

### `GpuDipoleConvolution`

This class is the central owner for the future GPU solver.

It currently owns:

- descriptor and boundary/discretization settings;
- padded-grid and half-spectrum layout;
- explicit forward, backward, and kernel FFT plans;
- explicit shared FFT workspace;
- real moment field grid;
- real field grid;
- moment spectrum;
- field spectrum;
- spectral kernel;
- full supercell matrix `H`;
- reciprocal matrix `2*pi*(H^-1)^T`;
- cell volume;
- point-Ewald energy storage.

It has memory estimators for:

- persistent real/spectral/kernel storage;
- geometry;
- point-energy storage;
- FFT work area using FFT-library `EstimateMany`.

`gpuSimulation.cpp` includes the CV6 FFT storage and workspace in the GPU memory preflight budget.

### Geometry conventions implemented

The descriptor is passed primitive `C1/C2/C3`.

It constructs:

```text
H = [N1*C1, N2*C2, N3*C3]
```

and uploads it.

It computes and uploads:

```text
Brec = 2*pi*(H^-1)^T
V = det(H)
```

The current builder requires a right-handed positive-volume cell.

This geometry correction matters especially for skew/triclinic cells. Using Cartesian grid indices as if the cell were orthogonal is incorrect.

### Device primitives currently implemented

All of these are dormant while `gpu_dipole_mode=OFF`.

1. Macro moment reduction

In `GpuHamiltonianCalculations::UpdateMacroMoments`:

```text
emomM(atom,ensemble) -> macroMoments(component,macrocell,ensemble)
```

using the CPU-exported one-based macrocell map.

2. Macro moment FFT packing

`GpuDipoleConvolution::packMacroMoments`

Maps:

```text
[component, macrocell, ensemble]
```

to FFT channels:

```text
[cell, component + 3*(basis + NA*ensemble)]
```

For `NA=1`, this is a straightforward regular-grid packing.

3. Forward transform

```text
forwardTransformMoments()
```

executes batched R2C.

4. Spectral contraction

```text
applySpectralKernel()
```

computes the complex `3×3×NA×NA` block multiplication.

5. Inverse transform

```text
inverseTransformFields()
```

executes C2R.

Important convention:

- No additional `1/N` normalization is applied after C2R.
- The intended lattice Ewald kernel must be defined for the raw inverse FFT sum.

6. Screened real-space primitive

```text
addRealSpaceField(alpha, cutoff, image_extent)
```

Current limitations:

- only `NA=1`;
- direct all-pairs macrocell loop;
- direct explicit periodic image loop;
- uses macrocell centres;
- excludes only exact self `(same cell, zero image)`;
- no Verlet/cell list yet;
- intended for validation, not performance.

7. Point self correction

```text
addPointSelfField(alpha)
```

adds:

```text
+4 alpha^3/(3 sqrt(pi)) M
```

8. Point-Ewald composition

```text
evaluatePointEwald(macro_moments, alpha, cutoff, image_extent)
```

currently composes:

```text
buildReciprocalEwaldKernel(alpha)
packMacroMoments(...)
forwardTransformMoments()
applySpectralKernel()
inverseTransformFields()
addRealSpaceField(...)
addPointSelfField(alpha)
reducePointEwaldEnergy()
```

9. Field scatter

```text
scatterPointFields(beff, one_based_cell_index, atom_count)
```

adds, rather than overwrites, the grid field into atom `beff`.

It is not currently invoked from production simulation.

10. Energy reduction and readback

```text
reducePointEwaldEnergy()
pointEwaldEnergies()
pointEwaldFields()
```

Energy is reduced as:

```text
-½ Σcell Mcell · Bcell
```

once per macrocell and ensemble.

Readback synchronizes the work stream first.

## Current external temporary harness

Locally added, external-to-CMake test-phase files:

- `tests/dipole_validation/gpu_ewald_driver.cpp`
- `tests/dipole_validation/run_gpu_ewald.py`

The driver:

- creates a `1×1×1`, `NA=1`, cubic cell;
- uploads one macrocell centre and one moment;
- calls `evaluatePointEwald`;
- reads back fields and energy;
- prints JSON.

The Python script:

- compiles the driver with `nvcc`;
- links it to `build/libasdlib.a`;
- runs it;
- compares/prints against the independent Python oracle.

This harness is intentionally outside CMake/CTest. It is not a production regression test yet.

## Important numerical failure found

The external GPU harness was run successfully and returned:

```text
GPU fields:
[ 0.5506404409927862,
 -0.11012808819855724,
  0.22025617639711448 ]

GPU energy:
-0.3303842645956717
```

The independent oracle gives:

```text
Oracle fields:
[ 0.1551403779550518,
 -0.031028075591010243,
  0.06205615118202068 ]

Oracle energy:
-0.09308422677303108
```

This is not a small convergence error. The GPU result is incorrect.

## Root cause of the failure

The current `buildReciprocalEwaldKernel(alpha)` is physically incomplete.

It directly samples the reciprocal Ewald tensor only on the FFT grid modes:

```text
-4*pi*exp(-k^2/(4 alpha^2))*k*k^T/(V*k^2)
```

and omits `k=0`.

For a `1×1×1` source grid, the FFT has only the zero mesh mode. Tin-foil correctly removes that mode.

Therefore the current GPU operator includes:

```text
Breal + Bself
```

but has essentially no reciprocal contribution.

The point self term alone is large and positive:

```text
4 alpha^3/(3 sqrt(pi)) M
```

and dominates the observed incorrect result.

The independent Ewald field includes nonzero reciprocal vectors outside the one-cell FFT mesh. Those missing terms provide the necessary compensating contribution.

This demonstrates that simply treating FFT-grid reciprocal modes as the complete Ewald reciprocal sum is wrong for this formulation.

## Required kernel redesign

The current reciprocal-only builder must not be enabled.

For the regular macrocell convolution formulation, build a complete periodic lattice tensor in real displacement space:

```text
K(Δg, a, b) =
    converged real-space Ewald contribution
  + converged reciprocal-space Ewald contribution
  + point self term only when Δg=0 and a=b
```

Then FFT that complete displacement-space kernel:

```text
Khat(k, a, b) = FFT[K(Δg, a, b)]
```

This has several benefits:

- For a `1×1×1` grid, the zero displacement kernel contains the full periodic self-image response, so the one-cell case is meaningful.
- The FFT performs the periodic lattice convolution of the already complete Ewald tensor.
- The macrocell grid size is not incorrectly treated as a reciprocal Ewald cutoff.
- The eventual `NA>1` basis-block kernel has a natural displacement-space representation.

For the initial `NA=1` slice:

1. Build host-side tensor `K(Δg)` for all periodic grid displacements.
2. For each displacement:
   - use macrocell centres or regular-grid geometry;
   - sum screened real-space images to convergence;
   - sum reciprocal vectors to convergence;
   - add self correction only at exactly zero source-target displacement.
3. Upload the resulting 9 real tensor components per displacement.
4. Use the existing kernel R2C plan to transform it into `kernel_fft`.
5. Run the existing pack/forward/contract/inverse path.
6. Compare fields and energies with the independent oracle.
7. Sweep:
   - alpha;
   - real image/cutoff;
   - reciprocal image cutoff;
   - grid size;
   - cubic and skew cells;
   - `M=1` and `M=4`.

Only after that should production runtime dispatch be enabled.

## Future performance work

The first full-kernel builder may be host-built and slow. That is acceptable for initialization and validation.

After correctness is established:

- move tensor construction to GPU if beneficial;
- replace direct real-space all-pairs/image loops with a GPU cell/Verlet list;
- cache/reuse kernel by geometry and Ewald parameters;
- extend from `NA=1` to `NA>1`;
- add finite-cell form factors and macrocell self-demagnetization for coarse blocks;
- wire energy into UppASD energy columns;
- enable `gpu_dipole_mode=PME3D` only after regression acceptance.

## Safety/status summary

Safe to use now:

```text
gpu_dipole_mode OFF
```

This is the default and normal CUDA calculations should run without invoking unfinished dipole physics.

Not safe to enable:

```text
gpu_dipole_mode PME3D
```

The runtime intentionally rejects it, and the currently built reciprocal-only kernel is numerically incorrect as demonstrated by the external harness.

## Relevant commits

Recent CV6 commits include:

```text
92d83c9  gpu: add CV6 dipole FFT layout contract
bf792e1  gpu: own CV6 dipole FFT resources
add9237  gpu: stage CV6 dipole cell matrix
e5337aa  gpu: add explicit CV6 dipole input contract
c7398f4  gpu: budget CV6 dipole FFT workspace
9c2bc0b  gpu: stage CV6 dipole supercell geometry
1b7d834  gpu: default CV6 dipole mode to off
e612782  gpu: stage CV6 reciprocal cell geometry
92a0fae  gpu: pack CV6 macro moments for FFT
7f31360  gpu: add CV6 dipole FFT execution primitives
592a84e  gpu: contract CV6 dipole spectral kernel
67637bc  gpu: build CV6 reciprocal Ewald kernel
50155c5  gpu: add CV6 dipole point self correction
0ba411e  gpu: add CV6 screened real-space primitive
657dd18  gpu: compose CV6 point Ewald field path
088fa3e  gpu: scatter CV6 point Ewald fields
7b4c9a4  gpu: reduce CV6 point Ewald energy
c4c1932  gpu: compose CV6 point Ewald energy
66845f1  gpu: synchronize CV6 Ewald energy readback
f449246  gpu: expose CV6 point Ewald fields
```

The temporary external harness files should be reviewed and committed separately if retained.
