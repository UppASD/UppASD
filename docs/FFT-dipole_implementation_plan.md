# FFT dipole implementation and validation plan

## 1. Purpose and current conclusion

This document turns the handoff in `docs/FFT-dipole_status.md` into an
implementation plan. It is based on a code audit of the current Fortran
macrocell path, the C/Fortran bridge, `GpuDipoleConvolution`, the GPU
Hamiltonian and measurement paths, the independent Python Ewald evaluator,
and both legacy CPU FFT implementations.

The central conclusion is:

> Do not repair `buildReciprocalEwaldKernel()` by merely adding a few more
> reciprocal modes. Replace the current runtime split with a complete,
> precomputed periodic lattice Green-function kernel and use the FFT only to
> apply that complete discrete convolution.

For UppASD's fixed regular lattice this is simpler and more accurate than a
general particle-mesh method:

```text
initialization:
    converged periodic Ewald tensor K(displacement, target basis, source basis)
    -> FFT and normalize K once

each field evaluation:
    atom moments -> macro/basis moments
    -> forward FFT
    -> block spectral contraction
    -> inverse FFT
    -> physical-unit scaling
    -> add field to atoms
```

There is no real-space all-pairs Ewald calculation on every spin-dynamics
step in this design. Real, reciprocal, point-self, and the selected surface
term all belong to one initialization-time discrete operator.

The first production slice should be deliberately narrow:

```text
execution:          GPU spin dynamics only
boundary:           3D periodic
surface:            tin foil / conducting
macro block:        1 x 1 x 1
basis:              NA = 1
precision gate:     fp64
source model:       point dipoles at the atom/macrocell positions
```

Only after that slice passes independent field, energy, indexing, unit, and
runtime tests should the implementation expand to `NA>1`, coarse blocks,
fp32, Monte Carlo, open boundaries, or slabs.

The maintainer decisions for WP0 are now frozen:

```text
public selector:      gpu_dipole_mode EWALD3D_FFT
default selector:     gpu_dipole_mode OFF
legacy CPU selector:  do_dip 0 is mandatory when EWALD3D_FFT is selected
accuracy selector:    gpu_dipole_tol, default 1.0e-10
construction inputs:  selected automatically from geometry and tolerance
```

`EWALD3D_FFT` is intentionally short enough for the currently declared
16-character Fortran input field while identifying both the 3D-periodic Ewald
physics and the FFT application method.

## 2. Frozen algorithm name and input decision

The staged name `PME3D` is misleading for the preferred regular-lattice
algorithm. Classical PME/P3M has:

- a mesh chosen independently from the particle positions;
- particle-to-mesh assignment;
- a lattice influence function;
- mesh-to-particle interpolation adjoint to assignment;
- mesh aliasing and mesh self-interaction corrections; and
- a direct, short-range real-space calculation at runtime.

The first UppASD solver needs none of those operations because the spins are
fixed on a regular Bravais lattice and can be represented as FFT channels
directly. It is an exact discrete periodic Ewald/lattice-Green-function FFT
solver, up to the explicitly converged initialization cutoffs.

The accepted public input is:

```text
gpu_dipole_mode EWALD3D_FFT
```

Replace the staged `PME3D` spelling rather than retain it as a physics alias.
If an existing development input uses `PME3D`, reject it with a migration
message directing the user to `EWALD3D_FFT`. A separate PME mesh is not
implemented, so `gpu_dipole_mesh != 0` remains unsupported by this mode. Do
not silently pretend that the macro grid is a general PME mesh.

A true dipolar P3M/PME implementation remains a useful later fallback for
irregular macrocell centres. It must be a separate work package with
assignment order, optimized influence function, adjoint interpolation, and
mesh self-correction. The relevant dipolar P3M reference is Cerdà,
Ballenegger, and Holm, J. Chem. Phys. 135, 184110 (2011),
<https://doi.org/10.1063/1.3657407>. The DOI currently written in
`GPU_FFT_DIPOLE_DESIGN.md` is incorrect and should be fixed.

## 3. Audit of the current repository

### 3.1 Components that are good foundations

The following work should be retained:

- The default GPU dipole mode is `OFF`, and non-`OFF` currently fails fast.
- GPU dipole selection has its own input contract instead of silently
  overloading the exchange convolution switch.
- The input mode/surface enum storage is persistent across the C ABI.
- The legacy CPU `do_dip=2` macrocell layout is separate from the new
  basis-resolved PME layout.
- `create_pme_macrocell_layout()` orders macrocells with basis as the fastest
  channel:

  ```text
  macro = basis + NA * cell
  ```

- The full periodic supercell matrix is formed as:

  ```text
  H = [N1*C1, N2*C2, N3*C3]
  ```

  and not confused with the primitive cell.

- The reciprocal columns are formed as `2*pi*(H^-1)^T`.
- FFT plans, persistent buffers, and a shared explicit FFT workspace have a
  clear owner and stream.
- The moment batch ordering is suitable:

  ```text
  component + 3*(basis + NA*ensemble)
  ```

- The kernel batch ordering is suitable:

  ```text
  row + 3*(column + 3*(target_basis + NA*source_basis))
  ```

- Macro moments are recomputed from current device moments rather than copied
  once from stale host data.
- The field scatter is additive, which is required for coexistence with
  exchange, anisotropy, DMI, and external fields.
- The independent Python Ewald evaluator has the correct conceptual
  decomposition and already caught a genuine architecture error.

### 3.2 Components that must be replaced or restricted

#### Incomplete reciprocal kernel

`buildReciprocalEwaldKernel()` samples only one reciprocal vector represented
by each FFT mode. The exact reciprocal Ewald sum contains every reciprocal
vector. Vectors separated by an FFT grid period alias into the same discrete
mode and must all be included.

For a `1 x 1 x 1` grid, every reciprocal vector aliases to the sole FFT bin.
Tin-foil removes only the physical vector `k=0`; it does not make that FFT bin
zero. This is why the current one-cell result lacks the full reciprocal
response.

#### FFT normalization

cuFFT and hipFFT use unnormalised forward and inverse transforms. With

```text
Mhat(q) = sum_g M(g) exp(-i q.g)
raw_inverse(F)(g) = sum_q F(q) exp(+i q.g),
```

the raw inverse of `FFT(K) * FFT(M)` is:

```text
Ngrid * convolution(K, M).
```

Therefore a displacement-space kernel cannot simply be FFT'd and stored.
Exactly one of the following is required:

1. scale the inverse field by `1/Ngrid`; or
2. store `FFT(K)/Ngrid`.

The recommended contract is option 2. Scale the kernel spectrum once at
construction, retain the raw C2R result at runtime, and test this with a
delta-source convolution before testing any dipole physics.

#### Unsupported partial edge blocks

The PME macrocell builder rounds the grid dimensions upward when an atomistic
extent is not divisible by its block size. The final block then has a
different width and centre. Such centres do not form a cyclic regular lattice,
so their interaction is not a translation-invariant block convolution.

The first periodic FFT implementation must require:

```text
N1 % block_size_x == 0
N2 % block_size_y == 0
N3 % block_size_z == 0
```

Partial-edge cases remain valuable layout tests, but they must be rejected by
the regular periodic FFT solver until an irregular P3M or another explicit
fallback exists.

#### GPU mode is still coupled to legacy `do_dip=2`

Although the input is described as independent, the current code:

- creates the PME map only if `do_macro_cells='Y'` or `do_dip=2`;
- stages the map only if `do_dip=2`;
- budgets it only if `do_dip=2`; and
- initializes `GpuDipoleConvolution` only if `do_dip=2`.

This would force users to enable the legacy CPU macrocell Hamiltonian merely
to allocate GPU data. It also risks duplicate physics and allocates the dense
CPU `Qdip_macro` matrix.

The correct ownership rule is:

```text
gpu_dipole_mode != OFF:
    create and export the GPU macrocell layout
    require legacy do_dip == 0

legacy do_dip != 0:
    use only the selected CPU dipole implementation
```

The two paths should be mutually exclusive during a GPU run unless a future
diagnostic mode explicitly computes one without applying it.

#### Physical prefactor is not available in the GPU descriptor

The production field must be in Tesla:

```text
field_prefactor = (1.0e-7 * mu_B) / alat^3
                = mu0 * mu_B / (4*pi*alat^3).
```

`alat` is not currently exported through `FortranData` or stored in
`SimulationParameters`. The dormant solver consequently produces only
dimensionless tensor units. Add `alat` to the bridge, compute the prefactor in
fp64, and apply it exactly once at the production boundary.

The kernel oracle should remain dimensionless. Unit conversion belongs in an
independent test layer so an erroneous physical constant cannot be hidden in
the Ewald comparison.

#### Production energy integration is absent

`point_energy` is a diagnostic total per ensemble. UppASD GPU measurement
currently stores six energy columns:

```text
exchange, anisotropy, DM, tensor, external, total
```

and `EnergyData::Dip` is unimplemented. The production path must:

- reduce `-0.5 * sum_macro M_macro . B_macro`;
- apply the same Tesla prefactor as the field;
- divide by the atom count for UppASD's per-atom energy convention;
- add the result to the total energy;
- expose it as the `Dip` column; and
- multiply by `mu_B/mRy` only in the existing measurement conversion layer.

Do not make the diagnostic oracle readback and production per-atom energy use
the same ambiguously named buffer.

#### `beff` and `eneff` must both be handled

The current scatter only updates `beff`. A bilinear dipole Hamiltonian has the
same field contribution in the effective and energy-field paths. Production
wiring must update both arrays where the GPU integrator/Monte Carlo code
expects them, without overwriting the existing Hamiltonian terms.

#### Monte Carlo is not in the first safe slice

Spin dynamics recomputes a global field at well-defined Hamiltonian
evaluations. The GPU Monte Carlo implementation performs trial updates in a
different schedule. A long-range field changes after every accepted moment
change, so using a field computed once for a parallel batch can violate the
intended Markov chain.

Until a correct long-range Monte Carlo update scheme is designed and tested,
non-`OFF` GPU dipoles must be rejected for GPU MC. Do not assume that sharing
`GpuHamiltonianCalculations` makes the runtime semantics correct.

#### Memory accounting is low

`GpuDipoleFftLayout::persistentBytes()` counts one real field grid although
both `moments_real` and `fields_real` are allocated. The construction buffer
is described but is not included in `estimateBytes()`. The peak budget must
include:

```text
all persistent arrays
+ FFT workspace
+ every construction array alive at the same time
+ energy and geometry buffers
```

The estimate and actual allocation list should be tested from the same layout
object to prevent drift.

#### Open/slab layouts are not executable yet

For open and 2D-periodic modes, `layout.real_cells` is the padded FFT size,
but `packMacroMoments()` currently expects the source tensor to contain that
many physical macrocells. It does not zero and embed the smaller active grid.
This is another reason to keep those modes rejected while the 3D periodic
slice is completed.

### 3.3 Legacy CPU issues discovered during this audit

These should not distract from the GPU vertical slice, but they define the
later CPU work:

- `do_dip=1` is a finite open-boundary direct sum and is the only CPU path
  currently validated analytically.
- `do_dip=2` is a finite dense macrocell projection; it does not implement a
  periodic Ewald sum.
- `do_dip=3` builds a zero-padded finite/open convolution regardless of the
  boundary flags. It should not be used as a periodic GPU oracle.
- In `fftdipole_fftw.f90`, the temporary kernel array is allocated as
  `(N1_pad,N2_pad,N2_pad)` instead of `(N1_pad,N2_pad,N3_pad)`.
- The FFTW output unpack order uses an `I3`-contiguous formula while input
  packing uses `I1` contiguous. This differs from the MKL implementation and
  is wrong for non-cubic grids.
- CPU FFTW and MKL should eventually share one indexing helper and one
  backend-neutral test suite rather than be repaired independently by visual
  inspection.

## 4. Revalidated physical contract

### 4.1 Geometry and units

All coordinates used by the dimensionless kernel are in UppASD lattice units.
Let:

```text
H = [h1, h2, h3] = [N1*C1, N2*C2, N3*C3]
V = det(H) > 0
Brec = 2*pi*(H^-1)^T.
```

The physical field is:

```text
B_Tesla = (1.0e-7 * mu_B / alat^3) * B_dimensionless.
```

Moments in `emomM` are numerical magnetic moments in units of `mu_B`.

### 4.2 Sign convention

For source moment `M` and displacement `r = r_target-r_source`, use:

```text
B(r) = Hessian(1/r) M
     = (3 r r^T/r^5 - I/r^3) M.
```

The energy of a symmetric linear field operator is:

```text
E = -0.5 * sum_i M_i . B_i.
```

Consequently:

```text
-dE/dM_i = B_i.
```

The finite-difference energy derivative is a required validation, not merely
a documentation identity.

### 4.3 Ewald decomposition

For `alpha > 0`, the screened real tensor is:

```text
Treal(r) = A(r) r r^T + B(r) I

A(r) =
    3 erfc(alpha*r)/r^5
  + 6 alpha exp(-alpha^2*r^2)/(sqrt(pi)*r^4)
  + 4 alpha^3 exp(-alpha^2*r^2)/(sqrt(pi)*r^2)

B(r) =
   -erfc(alpha*r)/r^3
   -2 alpha exp(-alpha^2*r^2)/(sqrt(pi)*r^2).
```

The reciprocal tensor is:

```text
Trec(k) =
    -4*pi/V
    * exp(-k^2/(4*alpha^2))/k^2
    * k k^T,                 k != 0.
```

Under the adopted field sign, the point-self field is:

```text
Tself = +4*alpha^3/(3*sqrt(pi)) I.
```

For tin-foil boundaries:

```text
Tsurface = 0.
```

The corresponding conventional Ewald self energy is negative,
`-2*alpha^3/(3*sqrt(pi))*|M|^2`, which is consistent with
`E=-0.5*M.Bself`.

For a vacuum spherical exterior, the surface energy would be:

```text
Esurface = 2*pi/(3*V) * |sum_i M_i|^2
```

and its field would be:

```text
Bsurface = -4*pi/(3*V) * sum_i M_i.
```

That option is not part of the first implementation and must continue to fail
fast.

### 4.4 Complete periodic displacement kernel

For regular coarse-grid cell index `g` and basis channel `a`, define:

```text
r(g,a) = R(g) + tau(a).
```

For block size one:

```text
R(g) = g1*C1 + g2*C2 + g3*C3
tau(a) = Bas(:,a).
```

For a uniform divisible coarse block, `R(g)` uses the coarse translation
vectors and `tau(a)` includes the constant centre offset of the block/basis
channel.

Let `d = g_target-g_source` modulo the coarse grid. Define:

```text
K(d,a,b) =
    sum_real_images Treal(R(d)+tau(a)-tau(b)-H*n)
  + sum_k_nonzero Trec(k) cos/sin phase for R(d)+tau(a)-tau(b)
  + delta(d=0,a=b) Tself
  + Tsurface.
```

This is a real `3 x 3 x NA x NA` block for every grid displacement. It must
satisfy the reciprocity identity:

```text
K(d,a,b) = transpose(K(-d,b,a)).
```

For a Hessian point tensor the Cartesian block is itself symmetric, but the
full basis/displacement reciprocity check is still required because it catches
target/source and phase indexing errors.

### 4.5 Discrete convolution and stored spectrum

The desired field is:

```text
B(g,a) = sum_h,b K(g-h,a,b) M(h,b).
```

For `Ngrid = G1*G2*G3`, construct:

```text
Khat_stored(q,a,b) = FFT[K(:,a,b)] / Ngrid.
```

Runtime then uses:

```text
Mhat = raw_R2C(M)
Bhat(q,a) = sum_b Khat_stored(q,a,b) Mhat(q,b)
B = raw_C2R(Bhat).
```

No further normalization is applied.

The construction should contain an explicit function named along the lines
of `normalizeKernelSpectrum(1/Ngrid)`. A comment saying that the kernel is
"defined for the raw inverse" is not precise enough to prevent a regression.

### 4.6 Why reciprocal aliases matter

For a regular grid, reciprocal vector indices can be written:

```text
n = q + (G1*l1, G2*l2, G3*l3),
```

where `q` is an FFT bin and `l` ranges over all integer alias images. Every
such vector has the same phase on Bravais grid translations, while basis
offsets retain their individual phase.

The normalized reciprocal spectrum can therefore be built directly as:

```text
Khat_rec_stored(q,a,b) =
    sum_l, n!=0 Trec(Brec*n)
        * exp(i * (Brec*n) . (tau(a)-tau(b))).
```

For a `1 x 1 x 1` grid, `q=0` and the sum still contains every `n!=0`.
This is the precise correction to the current failed builder.

### 4.7 Two kernel builders

Implement two builders with different purposes.

#### Builder A: independent/reference displacement builder

This is the first correctness implementation:

1. Loop over every displacement and basis pair.
2. Sum real images explicitly.
3. Sum reciprocal images explicitly.
4. Add self only for zero displacement and identical basis.
5. Check convergence by increasing both image ranges.
6. Store the complete real tensor.
7. FFT it and divide its spectrum by `Ngrid`.

It may be slow and host-only. It is acceptable for small validation systems
and should be written for clarity in fp64.

#### Builder B: production alias-summed builder

After Builder A and the GPU convolution pass:

1. Build the real-space displacement tensor in fp64.
2. Add the point-self block.
3. FFT this tensor and divide by `Ngrid`.
4. For each FFT bin, sum all reciprocal alias vectors to the requested
   convergence tolerance, including basis phases.
5. Add the reciprocal alias sum directly to the stored spectrum.
6. Enforce/check Hermitian symmetry and basis reciprocity.

Builder B avoids an `O(Ngrid * Nk)` explicit reciprocal sum during
initialization while remaining algebraically equivalent to Builder A.

Do not optimize Builder A into Builder B before the delta-convolution and
oracle gates pass. They are intended to be two independently comparable
formulations.

## 5. Macrocell model

### 5.1 Block-one model

For block size one, each basis channel contains one atom. This is the exact
point-dipole model on the regular crystal lattice and is the acceptance
baseline.

`NA=1` is only the first implementation slice, not a physics limitation.

### 5.2 Multi-basis block convolution

For `NA>1`, retain distinct channels:

```text
M(g,a) = atom moment at cell g, basis a          # block one
Bhat(q,a) = sum_b Khat(q,a,b) Mhat(q,b).
```

The GPU contraction indexing already supports this shape. The missing work
is kernel construction, macro-to-grid mapping validation, multi-basis field
scatter, and energy reduction over all `(g,a)`.

Required multi-basis tests include:

- a two-basis orthogonal cell with unequal basis moments;
- a two-basis skew cell;
- swapping basis labels and permuting the corresponding input moments;
- a source impulse in each basis channel and Cartesian component; and
- `M=4` ensembles with different moments.

### 5.3 Coarse blocks: use an explicit projected Hamiltonian

Treating a coarse cell as a point moment at its centre is an approximation
that omits the atomistic interactions inside the cell. A more defensible
coarse model is the block projection of the block-one periodic tensor.

Assume every atom in a macro channel shares a uniform moment:

```text
m_i = M_A / n_A,     i in macro channel A.
```

Define the projected block:

```text
Kcoarse(A,B) =
    1/(n_A*n_B)
    * sum_{i in A, j in B} Katom(i,j),
```

with the physical point self excluded for `i=j` but all distinct intra-cell
pairs and periodic images retained.

Then:

```text
B_A = sum_B Kcoarse(A,B) M_B
Ecoarse = -0.5 * sum_A M_A . B_A
```

is exactly the atomistic Hamiltonian restricted to cell-uniform moments and
average target fields. This naturally supplies the finite-cell form factor
and self-demagnetizing block instead of adding an ad hoc centre-point term.

For uniform divisible blocks this projection remains translation invariant
and can be FFT'd. Implement it only after block-one `NA>1` is correct.

The existing CPU `do_dip=2` averaging formula is conceptually related for a
finite sample, but it is not a periodic correctness oracle.

## 6. Runtime blueprint

### 6.1 Initialization

The production initialization sequence should be:

```text
validate input combination
validate regular grid and block divisibility
create/export GPU-only basis-resolved layout
validate macro_count == G1*G2*G3*NA
validate every atom maps to exactly one macro channel
validate channel centres against expected fractional regular positions
form H, Brec, V, and field prefactor in fp64
read gpu_dipole_tol and select internal Ewald construction parameters
build complete normalized kernel spectrum
run optional initialization invariants
create/reuse FFT plans and buffers
mark dipole backend ready
```

The backend-ready flag must be separate from descriptor/plans-ready. A valid
FFT allocation is not evidence that a physical kernel exists.

### 6.2 Each spin-dynamics Hamiltonian evaluation

The sequence should be:

```text
zero macro moments
reduce current emomM -> macro moments
pack macro moments into FFT batches
forward R2C
block spectral contraction
inverse C2R
scale dimensionless field by physical prefactor
scatter/add to beff
scatter/add to eneff where required

if measuring energy:
    reduce -0.5 * sum_macro M_macro . B_macro
    divide by Natom
    add to dipole energy column
    add to total energy column
```

Use the same post-C2R macro field for scatter and energy. Never reconstruct
energy from a different reciprocal expression.

### 6.3 Stream and synchronization rules

- All kernels and FFT executions use the Hamiltonian work stream.
- No host synchronization occurs during normal field evaluation.
- Diagnostic readback functions synchronize explicitly and remain test-only.
- Kernel construction must finish before `kernel_ready` becomes visible.
- Measurement kernels on another stream need an explicit event dependency if
  the existing simulation does not already serialize those streams.

### 6.4 Failure behavior

Reject, with a specific message:

- non-periodic boundaries in periodic mode;
- non-tin-foil surface in the first slice;
- `do_dip != 0` together with a GPU dipole mode;
- non-finite or non-positive `gpu_dipole_tol`;
- nonzero legacy `gpu_dipole_alpha` or `gpu_dipole_rcut` overrides for
  `EWALD3D_FFT`;
- nonzero `gpu_dipole_mesh` for the regular-lattice solver;
- non-divisible coarse blocks;
- irregular/missing centres;
- wrong macro count or population;
- MC execution;
- single precision before its acceptance gate;
- `NA>1` before the multi-basis gate;
- kernel construction that does not converge; and
- a kernel spectrum that violates finite/Hermitian/reciprocity checks.

Silent fallback is not acceptable for a requested long-range Hamiltonian.

## 7. Input contract and parameter selection

### 7.1 First enabled contract

The accepted user-facing contract is:

```text
gpu_dipole_mode     EWALD3D_FFT
gpu_dipole_surface  TINFOIL
gpu_dipole_tol      1.0e-10             # optional; this is the default
gpu_dipole_mesh     0 0 0
do_dip              0
BC                  P P P
```

`gpu_dipole_mode` is the GPU-specific selector and defaults to `OFF`.
Selecting `EWALD3D_FFT` requires legacy `do_dip=0`; this prevents the CPU and
GPU dipole Hamiltonians from being constructed or applied together.

`gpu_dipole_tol` is a dimensionless convergence target for the complete
periodic tensor constructed in fp64. Define convergence with a mixed
component-wise criterion:

```text
max_component |K_new-K_old|
    <= gpu_dipole_tol * max(1, max_component |K_new|)
```

and require the criterion for two consecutive real/reciprocal shell
expansions. The default `1.0e-10` targets ten-digit kernel construction while
leaving headroom above fp64 roundoff. Reject a non-finite or non-positive
tolerance. Requests too close to fp64 machine precision must fail with a
useful lower-bound message instead of looping indefinitely.

Do not judge convergence only after enlarging real and reciprocal ranges
together, because their errors could cancel accidentally. Hold the other
range fixed and bound the real and reciprocal tail changes independently,
with an explicit error-budget split (initially `tol/4` for each tail,
`tol/4` for the second-alpha cross-check, and `tol/4` numerical headroom).

The existing staged `gpu_dipole_alpha` and `gpu_dipole_rcut` inputs are not
part of the public `EWALD3D_FFT` contract. Keep their zero defaults during the
transition and reject nonzero values for this mode rather than silently
ignoring them. Explicit alpha, real cutoff, and reciprocal extents may remain
available through a test-only construction API for oracle comparisons and
diagnostics.

`gpu_dipole_mesh` remains reserved for a future true P3M implementation. It
does not control the exact regular-lattice convolution grid and must be all
zero for `EWALD3D_FFT`.

### 7.2 Automatic construction policy

Automatic construction is required before the first production mode can be
enabled; users do not tune Ewald split parameters. The builder should:

1. derive direct- and reciprocal-space length bounds from the full skew
   matrices `H` and `Brec`, not only from Cartesian box lengths;
2. choose a small deterministic set of internal alpha candidates expected to
   balance real and reciprocal setup work;
3. use analytic exponential-tail estimates to seed real and reciprocal shell
   extents;
4. expand and test the real and reciprocal tails independently until each
   passes the two-consecutive-shell budget;
5. verify the converged result with a second alpha candidate during the
   reference-builder stage;
6. choose the converged candidate with the lowest measured/estimated setup
   work; and
7. report the chosen alpha, shell extents, residual, and setup time as
   diagnostics without exposing them as required user inputs.

Alpha is a numerical split parameter, not part of the Hamiltonian. The
complete kernel must be invariant under changing it within
`gpu_dipole_tol`. Builder B may later optimize this search, but it must
preserve the same external tolerance contract and agree with Builder A.

## 8. Code-level blueprint

### 8.1 Fortran layout and lifecycle

Files:

```text
source/CoarseGraining/macrocells.f90
source/uppasd.f90
source/chelper.f90
source/Input/inputdata.f90
source/Input/inputhandler.f90
```

Required changes:

- Create the GPU layout when `gpu_dipole_mode != OFF`, independently of
  `do_macro_cells` and `do_dip`.
- Do not create the legacy dense macrocell layout solely for GPU use.
- Add a GPU-layout-only cleanup path; current cleanup reaches
  `release_pme_macrocell_layout()` only through legacy macrocell cleanup.
- Replace the staged `PME3D` enum spelling with `EWALD3D_FFT`.
- Add persistent `gpu_dipole_tol` input storage with default `1.0e-10` and
  export it through the bridge.
- Export `alat`.
- Validate mode, tolerance, legacy-override zeros, surface, BC, block
  divisibility, and `do_dip` conflict before the GPU allocates memory.
- Export or derive the uniform block translations/basis offsets explicitly.
  Do not rely only on averaged floating-point centres without checking them.
- Keep all bridge scalars in persistent module storage.

### 8.2 Bridge and simulation parameters

Files:

```text
source/gpu_files/fortranData.hpp
source/gpu_files/fortranData.cpp
source/gpu_files/gpuStructures.hpp
source/gpu_files/gpuSimulation.cpp
```

Required changes:

- Add `alat`.
- Add `gpu_dipole_tol` to the bridge, simulation parameters, descriptor, and
  construction diagnostics.
- Stage GPU macro data based on GPU mode, not `do_dip=2`.
- Remove the incorrect `macroCount > N` assumption if a later representation
  can have empty/deposition channels; retain equality/population checks for
  the current basis-resolved map.
- Make one function construct the descriptor used by both allocation and
  memory preflight so those paths cannot diverge.
- Correct the persistent real-grid byte count.
- Include construction peak memory.
- Add overflow-checked addition/multiplication to the top-level budget; the
  current top-level sum can itself overflow.

### 8.3 Kernel construction API

Suggested new host-side types:

```cpp
struct DipoleKernelSettings {
    double tolerance = 1.0e-10;
};

// Internal/test-visible values selected automatically from settings and
// geometry; these are not required production input keywords.
struct DipoleEwaldParameters {
    double alpha;
    std::array<int,3> real_images;
    std::array<int,3> reciprocal_images;
    double real_cutoff;
    double tolerance;
};

struct DipolePeriodicGeometry {
    std::array<double,9> H;
    std::array<double,9> Brec;
    double volume;
    GpuDipoleGridShape grid;
    unsigned basis;
    std::vector<std::array<double,3>> basis_offsets;
};

struct DipoleKernelDiagnostics {
    double max_reciprocity_error;
    double max_hermitian_error;
    double max_alpha_difference;
    int real_shells;
    int reciprocal_shells;
};
```

Suggested functions:

```cpp
std::vector<double> buildPeriodicEwaldDisplacementKernel(...);
std::vector<std::complex<double>> buildAliasSummedSpectrum(...);
DipoleKernelDiagnostics validatePeriodicKernel(...);
```

Keep these functions free of CUDA/HIP calls so they can be unit tested on
CPU-only CI. Convert to `real` only at the upload boundary.

### 8.4 `GpuDipoleConvolution`

Replace the public runtime API:

```text
buildReciprocalEwaldKernel
addRealSpaceField
addPointSelfField
evaluatePointEwald
```

with a clearer state machine:

```cpp
bool initiatePlansAndStorage(...);
void buildReferencePeriodicKernel(...);  // validation/small systems
void buildProductionPeriodicKernel(...); // later alias builder
bool kernelReady() const;
void evaluate(const GpuTensor<real,3>& macro_moments);
void addFields(...);
void accumulateEnergy(...);
```

The old split primitives may remain temporarily as private diagnostics, but
they must not be callable by production dispatch.

Add a transient real kernel tensor:

```text
[real_cell, kernel_batch]
```

execute the existing batched kernel R2C plan, then scale all spectral elements
by `1/Ngrid`. Release the transient tensor after construction.

Add explicit state flags:

```text
plans_ready
kernel_ready
physical_scale_ready
```

### 8.5 GPU Hamiltonian integration

Files:

```text
source/gpu_files/gpuHamiltonianCalculations.hpp
source/gpu_files/gpuHamiltonianCalculations.cpp
```

Required changes:

- Replace the unconditional non-`OFF` rejection only after all vertical-slice
  gates pass.
- Select the dipole backend from `gpu_dipole_mode`, never from inferred BC
  alone.
- Keep `refreshMacroMoments` tied to the GPU dipole backend.
- Evaluate the dipole after the base Hamiltonian has written `beff/eneff` so
  its scatter is additive.
- On measurement steps, put dipole energy in its own column and total.
- Reject MC at the call-site or simulation initialization, not halfway through
  a run.

### 8.6 Measurement integration

Files:

```text
source/gpu_files/gpuStructures.hpp
source/gpu_files/measurement/kernels.hpp
source/gpu_files/measurement/kernels.cpp
source/gpu_files/measurement/gpuMeasurement.cpp
source/gpu_files/measurement/measurementData.h
```

Recommended minimal change:

- Expand `energyM` with a dipole column while preserving the existing total
  index or use named constants instead of raw indices.
- Increase `N_ENERGY_TYPES`.
- Map the new column to `EnergyData::Dip` and `std_Dip`.
- Ensure total already contains the dipole contribution exactly once.
- Update Binder energy use only if it reads a total that previously omitted
  dipoles.
- Add a test with exchange disabled and dipoles enabled, and another with
  both enabled, to catch missing and double-counted totals.

## 9. Validation architecture

### 9.1 Test layers

Use four independent layers:

1. **Formula/oracle tests:** Python, no UppASD or GPU code.
2. **Host kernel tests:** C++ or Python/C++ driver, no GPU required.
3. **GPU operator tests:** standalone CMake/CTest executable using the real
   GPU class.
4. **End-to-end UppASD tests:** input files, Tesla fields, mRy/atom energy,
   interaction with the rest of the GPU Hamiltonian.

CPU output is never the periodic correctness oracle.

### 9.2 Strengthen the Python oracle

Add:

- automatic convergence checks using consecutive real/reciprocal shell
  expansions;
- numerical finite-difference checks of `-dE/dM = B`;
- reciprocal identity `H^T*Brec = 2*pi*I`;
- explicit pair reciprocity;
- alpha sweeps for cubic, rectangular, and skew cells;
- non-cubic image extents;
- a second implementation of the complete displacement kernel;
- basis-offset cases rather than only arbitrary two-particle positions; and
- clear total-versus-per-particle energy metadata.

Golden values protect against accidental edits but are not themselves an
independent derivation. Retain them while adding the derivative and
cross-formulation checks.

### 9.3 Mandatory delta-convolution tests

Before comparing Ewald values, fill a small real kernel with recognizable
non-symmetric values that still obey real FFT requirements. Put a unit impulse
in one source basis/component and verify every output entry.

Cover:

```text
grid:       1x1x1, 2x1x1, 2x3x1, 3x2x2
basis:      1 and 2
ensemble:   1 and 4
source:     every basis and Cartesian component
```

This isolates:

- `n1`-fast indexing;
- R2C half-spectrum layout;
- kernel batch ordering;
- target/source basis ordering;
- C2R normalization;
- ensemble isolation; and
- wraparound displacement sign.

No physics test can substitute for this gate.

### 9.4 Periodic physics matrix

Minimum fp64 periodic cases:

| Case | Grid/basis | Cell | Purpose |
|---|---:|---|---|
| one cubic source | `1x1x1`, `NA=1` | cubic | reciprocal aliases at FFT bin zero |
| two axial cells | `2x1x1`, `NA=1` | orthogonal | displacement sign |
| non-cubic grid | `2x3x1`, `NA=1` | orthogonal | indexing/normalization |
| skew two-source | small, `NA=1` | triclinic | reciprocal geometry |
| two-basis cell | `2x2x2`, `NA=2` | orthogonal | block contraction |
| skew two-basis | small, `NA=2` | triclinic | complex basis phases |
| four ensembles | any above, `M=4` | both | batch isolation |
| zero moments | any | any | exact zero/no NaN |
| global sign flip | any | any | field odd, energy even |

For every case compare:

- each macro field component;
- total dimensionless energy;
- `-dE/dM`;
- alpha invariance;
- translation invariance;
- reciprocity;
- host reference kernel versus alias-summed kernel; and
- GPU versus host kernel application.

### 9.5 Initial numerical tolerances

Use explicit, separately reviewed tolerances:

```text
oracle shell convergence, fp64:      <= 1e-13 component change
host builder vs oracle, default tol: <= 2e-10 relative/absolute scale
GPU fp64 vs host convolution:        <= 5e-11 on small cases
end-to-end fp64 vs oracle:            <= 3e-10 at default tolerance
energy/field identity:               <= 1e-9 with tuned finite difference
reciprocity/Hermitian residual:       <= 1e-12 host fp64
GPU fp32 later:                      measured target, initially <= 5e-5
```

Use a mixed absolute/relative comparison near zero. Tighten or relax only
from measured convergence evidence, never to make one unexplained case pass.

### 9.6 End-to-end gates

An enabled production mode must pass:

- Tesla field comparison using the UppASD `alat` and `mu_B` values actually
  read by the run;
- `Dip` energy in mRy/atom;
- total energy equals the sum of printed columns;
- dipole-only and exchange-plus-dipole runs;
- `M=1` and `M=4`;
- multiple Hamiltonian evaluations with changed moments;
- no stale macro moment between predictor/corrector stages;
- non-32-aligned atom counts;
- CUDA memcheck/compute-sanitizer;
- HIP build and execution;
- release/re-init lifecycle;
- memory preflight just below and above a controlled budget; and
- default `OFF` behavior unchanged.

The local audit environment could compile the standalone CUDA driver but has
no CUDA-capable device, so runtime GPU validation must run on a GPU worker.

## 10. Implementation work packages and gates

### Work-package checklist

Tick a top-level box only after that work package's stated gate has passed,
not merely when its implementation has started.

- [x] **WP0 — Contract frozen:** `EWALD3D_FFT`, narrow fp64 SD scope,
  `do_dip=0`, and tolerance-driven automatic construction.
- [x] **WP1 — Independent oracle and red tests** (Luna; depends on WP0).
- [x] **WP2 — GPU layout ownership, input, units, and memory contract**
  (Terra; depends on WP0 and may run in parallel with WP1).
- [x] **WP3 — Complete fp64 host reference kernel with automatic
  tolerance convergence** (Terra; depends on WP1, consumes WP2 input
  contract).

  Status: complete.  CPU-only `dipole-ewald-host-builder` CTest covers
  cubic, non-cubic, and skew `NA=1` Luna-oracle values plus alpha invariance
  and reciprocity; Luna's independent Python oracle self-check also passes.
- [x] **WP4 — GPU convolution and normalization correctness** (Terra
  implementation, Luna acceptance; depends on WP3).
  - [x] **WP4a — Standalone FFT plumbing:** CUDA fp64 delta and periodic
    Builder A CTests pass, including compute-sanitizer memcheck with zero
    errors.
  - [x] **WP4b — Independent GPU acceptance closure:** direct independent
    oracle comparisons, spectral validation, metamorphic GPU tests, and the
    reciprocal-test helper correction described below.

  Status: complete for the accepted CUDA fp64 `NA=1` WP4b scope. Terra's
  CUDA CTest and compute-sanitizer run passed; the independent Python oracle
  self-check passed on commit `ab4c5518a7b22078f69900ba028073e280710468`.
  Maximum GPU-vs-independent-oracle errors were `3.89e-16` field and
  `1.81e-16` energy; maximum GPU-vs-host-convolution errors were
  `1.78e-14` field and `1.14e-13` energy. Alpha, translation, sign-flip,
  zero-moment, changed-moment, M=4 isolation, reciprocity, Hermitian, and
  compute-sanitizer gates passed. HIP compilation/runtime remains deferred by
  maintainer decision while ROCm access is unavailable. Independent `NA>1`
  physics coverage remains WP6, and production wiring remains WP5.
- [x] **WP5 — Production fp64 `NA=1` spin-dynamics field and energy wiring**
  (Terra; depends on WP4b acceptance).

  Closed 2026-07-26. CUDA fp64 `ctest --test-dir build_ab --output-on-failure`
  passed all 6 tests, including `dipole-gpu-wp5-e2e`; `compute-sanitizer
  --tool memcheck --error-exitcode 1 build_ab/bin/dipole_gpu_fft_tests` reported
  `ERROR SUMMARY: 0 errors`. WP5 E2E maximum errors were
  `4.3719150199989277e-40 T` (field), `3.8376175788916486e-43 mRy/atom`
  (energy), and `2.6189900457383946e-40` (exchange-plus-dipole delta).
- [x] **WP6 — Multi-basis block-one support** (Terra implementation, Luna
  acceptance; depends on WP5).
- [x] **WP7 — Optimized reciprocal-alias builder and automatic-construction
  performance** (Terra; depends on cross-builder tests from WP3/WP4).

  Closed 2026-07-27.  Builder B constructs real+self in displacement space,
  adds the reciprocal alias sum directly to the normalized R2C spectrum, and
  remains cross-checked against Builder A on the cubic, non-cubic, skew, and
  two-basis validation matrix.  CPU cross-builder spectrum differences were
  at most `8.38e-17`; the host benchmark records Builder-A/B setup work and
  construction storage for every case.  CUDA fp64 CTest passed
  `dipole-ewald-host-builder`, `dipole-gpu-fft-convolution`, and
  `dipole-gpu-wp5-e2e`; compute-sanitizer memcheck reported `ERROR SUMMARY: 0
  errors`.
- [x] **WP8 — Coarse projected Hamiltonian** (Terra implementation, Luna
  acceptance; depends on WP6 and preferably WP7).

  Closed 2026-07-27.  Uniform divisible blocks now use the exact projection
  of the complete block-one Ewald operator with `M_block=sum_i m_i`; no ad hoc
  point-macrocell self-demagnetizing correction is added.  Projected Builder B
  applies the matching reciprocal block form factor and agrees with Builder A
  to `2.78e-17` in the host spectrum check.  The deterministic non-uniform
  convergence fixture reports field/energy-per-atom errors of
  `5.654e-02/1.623e-03` for block width four, `1.973e-02/3.431e-04` for width
  two, and zero for block one.  CUDA fp64 CTest passed all six tests, including
  the host builder, GPU convolution, and `NA=1/2` SD coarse-block equivalence;
  partial edge blocks are rejected before allocation.
- [x] **WP9 — Performance, fp32, and CUDA/HIP parity** (depends on the
  relevant accepted physics scope).
- [ ] **WP10 — Open-island/racetrack `OPEN_FFT` mode.**  The detailed,
  independently delegable execution blueprint is
  `docs/WP10_OPEN_FFT_BLUEPRINT.md`.
- [ ] **WP11 — Legacy CPU dipole audit and repair** (begin after WP10.1
  freezes an independent finite/open reference; retain the accepted periodic
  oracle for regression coverage).
- [ ] **WP12 — True 2D-periodic dipole mode** (separate derivation and
  acceptance; proposed selector `EWALD2D_FFT`).

When a package is active, record the assignee, branch/commit, test command,
and acceptance report immediately below its checklist item rather than
changing the meaning of `[x]`.

### WP0 — Freeze the contract — complete

Owner: maintainer decision.

Frozen decisions:

- Rename the staged `PME3D` selector to `EWALD3D_FFT`.
- Use the GPU-specific `gpu_dipole_mode` keyword, default `OFF`.
- First production scope is GPU SD, P/P/P, tin foil, block one, `NA=1`, fp64.
- `EWALD3D_FFT` requires legacy `do_dip=0`.
- Use `gpu_dipole_tol` with default `1.0e-10`.
- Select alpha and real/reciprocal construction extents automatically; they
  are diagnostics/test API values rather than required production inputs.

Gate:

- Maintainer decisions are recorded in this plan. Synchronizing the older
  design/status documentation and implemented input parser belongs to
  WP1/WP2.

### WP1 — Harden the independent oracle and create red tests

Suggested owner: Luna.

Deliverables:

- Oracle derivative, convergence, geometry, and reciprocity tests.
- Complete displacement-kernel generator for small grids.
- Delta-convolution expected data.
- Correct literature DOI.
- Synchronize the frozen `EWALD3D_FFT` name, first-scope statement, and
  tolerance contract into `GPU_FFT_DIPOLE_DESIGN.md`.
- Tests that demonstrate why the current reciprocal-only GPU builder fails.

Gate:

- All oracle self-tests pass without building UppASD.
- At least one test fails against the current GPU-kernel formulation by
  construction.

### WP2 — Decouple GPU layout ownership and add unit data

Suggested owner: Terra.

Status: complete.  CUDA acceptance run (2026-07-23) printed the 40x40x40
descriptor geometry, `gpu_dipole_tol=1.0e-10`, the fp64 physical prefactor,
and the persistent/peak budget before reaching the intentional
kernel-not-ready guard.  CPU OFF regression and invalid input rejection tests
also passed; HIP runtime remains a later parity gate.

Deliverables:

- GPU layout creation based on GPU mode.
- Mutual exclusion with legacy `do_dip`.
- GPU-only cleanup.
- `alat` bridge and physical prefactor storage.
- block-divisibility and BC validation.
- shared descriptor construction for budget and runtime.
- corrected memory estimate.

Gate:

- `gpu_dipole_mode OFF` regression unchanged.
- A requested valid mode reaches the existing kernel-not-ready guard with a
  complete descriptor.
- Invalid combinations fail before device allocation.

### WP3 — Host complete-kernel builder

Suggested owner: Terra.

Deliverables:

- Pure fp64 host Builder A.
- Complete `K(d)` for `NA=1`.
- `gpu_dipole_tol`-driven automatic alpha/shell selection and convergence
  diagnostics.
- reference tensor upload format matching the current kernel batches.
- CPU-only unit tests against the Python oracle.

Gate:

- All small grids match the oracle at the agreed tolerance.
- Alpha invariance and reciprocity pass.

### WP4 — GPU convolution correctness

WP4 is an umbrella gate.  WP4a records the accepted standalone FFT plumbing;
WP4b closes the independent-acceptance gaps before production wiring begins.

#### WP4a — Standalone FFT plumbing — complete

Suggested implementation owner: Terra.
Suggested independent reviewer: Luna.

Deliverables:

- Real kernel staging allocation.
- kernel R2C execution.
- explicit `1/Ngrid` spectral normalization.
- delta-convolution CTest executable.
- CMake-integrated periodic Ewald GPU driver.
- replacement of the temporary external compile script.

Gate:

- Delta-convolution matrix passes first.
- Same-kernel direct convolution passes for `1x1x1`, non-cubic, skew, and
  `M=4` cases on CUDA fp64.
- HIP builds; HIP runtime passes on available hardware.

These tests prove the FFT normalization, displacement indexing, batch
ordering, energy contraction, and memory safety.  They remain useful but do
not by themselves prove that Builder A and the GPU together implement the
independent periodic Hamiltonian.

#### WP4b — Independent GPU acceptance closure

Suggested repair owner: Terra.
Suggested independent acceptance owner: Luna.

Scope:

- Keep the WP4a same-kernel direct-convolution tests as plumbing tests.
- Fix the reciprocal-cell helper in `test_gpu_dipole_fft.cpp`: the affected
  cofactor must use `h[1]`, not `h[5]`.  Add a skew regression geometry for
  which those entries differ and the old expression fails.
- Add versioned golden fixtures generated by the independent Python oracle.
  Fixtures must contain geometry, moments, dimensionless fields, and total
  dimensionless energy.  C++ expected values must not be derived from
  Builder A's returned kernel, `directConvolution(complete.kernel, ...)`,
  GPU output, or a legacy CPU `do_dip` implementation.
- Compare every GPU field component and total energy directly with those
  independent values for `1x1x1`, `2x1x1`, `2x3x1`, a genuinely skew cell,
  and `M=4`, all in the accepted `NA=1` fp64 scope.
- Exercise the complete construction-to-GPU path with two valid explicit
  alpha choices through the test-only builder API and demonstrate invariant
  fields and energy.
- Add GPU properties for lattice translation, global moment sign flip, zero
  moments, ensemble isolation, and two consecutive evaluations with changed
  moments.  The last test covers stale standalone buffers; WP5 must repeat
  it through production moment aggregation.
- Before setting `kernel_ready`, reject non-finite spectra and validate the
  real-transform conjugacy and Cartesian/basis block Hermiticity implied by
  `K(d,a,b)=K(-d,b,a)^T`.  Record the maximum residual.  Initialization-only
  diagnostics/readback are acceptable; normal field evaluation must remain
  synchronization-free.
- Rerun the CUDA fp64 suite and compute-sanitizer.  HIP runtime remains
  deferred until suitable hardware is available.

Gate:

- The independent-oracle field and energy matrix passes the section 9.5
  tolerances without any expected value being generated from the uploaded
  kernel.
- Cross-alpha GPU fields and energies agree within the construction/error
  budget.
- Translation, sign-flip, zero-moment, changed-moment, and `M=4` isolation
  properties pass.
- Finite, reciprocity, and spectral Hermitian checks pass before
  `kernel_ready`, with reciprocity/Hermitian residual `<= 1e-12` in host
  fp64 diagnostics.
- The reciprocal helper regression fails with the old index and passes with
  the corrected expression.
- Compute-sanitizer reports zero errors.
- The production guard remains enabled.  No `beff`, `eneff`, `Dip`, or total
  energy wiring belongs to WP4b.

Macro-centre layout validation is a prerequisite of WP5 production
initialization.  Independent `NA>1` basis/source oracle coverage remains the
WP6 gate; internal basis-two delta tests do not complete WP6.

### WP5 — Production SD field and energy wiring

Suggested owner: Terra.

Deliverables:

- kernel-ready runtime state.
- SD-only dispatch.
- physical Tesla scale.
- additive `beff`/`eneff` scatter.
- per-atom dipole and total energy integration.
- end-to-end fixtures and output checks.

Gate:

- End-to-end field and energy match the independent oracle.
- Exchange-plus-dipole total has neither omission nor double counting.
- Repeated changed-moment evaluations show no stale data.

Only at this gate may the narrow first production mode stop failing fast.

### WP6 — Multi-basis block-one support

Suggested implementation owner: Terra.
Suggested validation owner: Luna.

Deliverables:

- basis phases in Builder A.
- `NA>1` scatter and energy reduction.
- basis permutation tests.
- multi-basis orthogonal/skew fixtures.

Gate:

- Full field/energy matrix passes for `NA=2`, `M=1/4`.

### WP7 — Optimized alias builder and automatic-construction performance

Suggested owner: Terra after Luna freezes cross-builder tests.

Deliverables:

- Builder B reciprocal alias sums.
- equivalence tests against Builder A.
- optimized automatic alpha/shell search preserving the already accepted
  `gpu_dipole_tol` semantics.
- setup-time and memory benchmarks.
- optional geometry/kernel cache with a complete cache key.

Gate:

- Builder A and B agree over the full validation matrix.
- Initialization cost is acceptable for representative systems.

### WP8 — Coarse projected Hamiltonian

Suggested implementation owner: Terra.
Suggested validation owner: Luna.

Deliverables:

- uniform divisible block projection.
- intra-block/self-demagnetizing tensor.
- block-size convergence suite.
- explicit rejection of partial edges.

Gate:

- block one is recovered exactly.
- coarse results match direct projection for cell-uniform moment states.
- approximation error versus unrestricted atomistic states is documented.

Status: complete 2026-07-27; the full CUDA fp64 CTest gate passed.

### WP9 — Performance and precision

Status: complete for the accepted CUDA scope (2026-07-27).  The implementation now has portable
event-backed phase timing for macro reduction, pack, forward FFT, spectral
contraction, inverse FFT, scatter, and energy reduction.  The standalone
`dipole_gpu_fft_benchmark` target reports one-time setup, steady-state device
time, and FFT-buffer memory for parameterized grid/basis/ensemble sweeps.

The standalone convolution/oracle suite also builds in `SINGLE` precision.
It compares fp32 storage against the existing fp64 host/oracle authority with
a `5e-5` mixed operator/spectrum budget.  CUDA fp32 acceptance was verified
on 2026-07-27 with
`ctest --test-dir build_gpu_wp9_fp32 --output-on-failure -R
'^(dipole-gpu-fft-convolution|dipole-gpu-wp5-e2e)$'`.

The gate includes the actual block-reduced `energyM` path at a dimensionless
O(1) scale and the real UppASD SD executable: Tesla `beff`/`eneff` scatter,
mRy/atom measurement output, predictor/corrector evaluation, M=4 shape,
non-warp-aligned N=3, additive exchange-plus-dipole accounting, and input
rejections.  The E2E writer comparison uses a `5e-5` relative budget with a
`1e-42` physical absolute floor in fp32; fp64 retains its `2e-8` writer budget
and its tighter standalone thresholds.  `EWALD3D_FFT` is therefore enabled
for the accepted CUDA fp32 periodic-SD scope.

CUDA memcheck completed on 2026-07-27 for both `dipole_gpu_fft_tests` fp32
and fp64 binaries with `ERROR SUMMARY: 0 errors`.  The focused CUDA fp32
host-builder/convolution/UppASD-E2E CTest gate passed 3/3; the complete fp64
CTest suite passed 6/6.  The separate fp32 Kagome tensor regression baseline
is not a dipole case and remains outside this FFT-dipole gate.

Builder-B widens the stored cell matrix and recomputes its reciprocal matrix
in host fp64 before Ewald construction, so the fp64 geometry check validates
a mathematically consistent pair rather than a rounded fp32 reciprocal
vector.  The generic Fortran-double-to-fp32 host staging path also preserves
null optional Fortran buffers, matching the fp64 pointer contract.  Neither
change relaxes the fp64 builder or accuracy gates.  HIP runtime parity remains
a deferred backend follow-up; it was not available for this CUDA acceptance
closure.

Builder-B host tensors are cached only in-process and only for an exact
immutable construction key (builder revision, storage precision, tolerance,
cell/grid/block geometry, basis, and offsets).  Device spectra and FFT plans
remain owned by each convolution instance.  There is no persistent/on-disk
cache.  The GPU test verifies a cache miss followed by an equivalent cache
hit; cache policy expansion is contingent on the setup-time sweep.

Suggested measurement commands:

```text
cmake --build build_gpu --target dipole_gpu_fft_benchmark
build_gpu/bin/dipole_gpu_fft_benchmark --grid 64 64 64 --basis 1 --ensembles 1
build_gpu/bin/dipole_gpu_fft_benchmark --grid 63 64 65 --basis 2 --ensembles 4

# CUDA precision gates.  Keep the FFT-dipole gate separate from unrelated
# stochastic fp32 GPU-regression baselines.
ctest --test-dir build_gpu_wp9_fp32 --output-on-failure \
  -R '^(dipole-ewald-host-builder|dipole-gpu-fft-convolution|dipole-gpu-wp5-e2e)$'
ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
  -R '^(dipole-ewald-host-builder|dipole-gpu-fft-convolution|dipole-gpu-wp5-e2e)$'
compute-sanitizer --tool memcheck --error-exitcode 1 build_gpu_wp9_fp32/bin/dipole_gpu_fft_tests
compute-sanitizer --tool memcheck --error-exitcode 1 build_gpu_wp9_fp64/bin/dipole_gpu_fft_tests

cmake -S . -B build_hip -DUPPASD_GPU_BACKEND=HIP -DUPPASD_PRECISION=DOUBLE -DBUILD_TESTING=ON
cmake --build build_hip --target dipole_gpu_fft_tests dipole_gpu_fft_benchmark
ctest --test-dir build_hip --output-on-failure -R 'dipole-(ewald-host-builder|gpu-fft-convolution)'
```

Record CUDA and HIP results separately.  Parity is satisfied by both backends
meeting the same documented field/energy budget, not by bitwise-identical FFT
output; include backend, compiler/runtime, device, and precision in every
published sweep row.

Deliverables:

- fp32 acceptance/error budget.
- FFT size performance sweep.
- macro reduction/scatter profiling.
- kernel cache policy.
- CUDA/HIP parity.

Gate:

- No precision mode is enabled without a documented field and energy error.

### WP10 — Open finite mode

Implement `OPEN_FFT` as a zero-padded linear convolution with correct
active-grid embedding and finite/open analytic references.  The sole source
representation is the existing basis-resolved macrocell map; the atomic
baseline means block `(1,1,1)` with `NA=1`.

The detailed task gates and separate Sol, Terra, and Luna prompts are in
`docs/WP10_OPEN_FFT_BLUEPRINT.md`.

### WP11 — Legacy CPU audit and repair

Begin after WP10.1 has frozen an independent finite/open oracle.  The accepted
GPU periodic oracle remains required for regression coverage but is not the
reference for the finite `do_dip=1/2/3` Hamiltonians.

Order:

1. Freeze `do_dip=1` finite analytic tests.
2. Define the intended physics of `do_dip=2` explicitly as a finite projected
   macrocell Hamiltonian.
3. Add a backend-neutral direct finite convolution oracle for `do_dip=3`.
4. Add non-cubic grids that expose the FFTW allocation/unpack bugs.
5. Make FFTW and MKL use common pack/unpack/index helpers.
6. Compare FFTW/MKL to the analytic finite oracle, not only to each other.
7. Decide whether a CPU periodic Ewald FFT mode should share the new host
   kernel builder rather than changing the semantics of legacy `do_dip=3`.

### WP12 — True 2D-periodic mode

Derive and validate a true 2D-periodic dipolar Green function or a complete
dipolar layer-corrected formulation.  Plain vacuum padding of the 3D Ewald
mode is not sufficient.  The proposed selector is `EWALD2D_FFT` if the
regular-lattice implementation precomputes a complete 2D-periodic operator
rather than implementing a general particle mesh.

The retained physics analysis and preliminary checklist are in
`docs/WP10_OPEN_FFT_BLUEPRINT.md`.

## 11. Delegation prompts

The prompts below are intentionally self-contained. Give each agent one prompt
at a time and preserve the work-package dependencies. After WP0 is frozen,
Prompt A and Prompt B may run in parallel because they own independent
test/documentation and bridge/lifecycle scopes; later prompts consume their
results. Require a fresh `git status` before edits because the current
workspace contains many unrelated untracked files.

### Prompt A — Luna: oracle and red-test package

```text
You are validating UppASD's new GPU periodic FFT dipole solver. Work only on
the independent/test side; do not enable production dispatch.

Read completely:
- docs/FFT-dipole_status.md
- docs/FFT-dipole_implementation_plan.md, especially sections 4, 9, and WP1
- tests/dipole_validation/periodic_ewald_reference.py
- tests/dipole_validation/test_periodic_ewald_reference.py

Tasks:
1. Harden the independent fp64 Ewald oracle with automatic shell-convergence
   checks, H^T*Brec=2*pi*I, reciprocity, and a finite-difference check that
   -dE/dM equals the returned field.
2. Add an independent small-grid displacement-kernel generator using
   K(d,a,b)=complete real+reciprocal+self Ewald.
3. Add cases for 1x1x1, 2x1x1, a non-cubic grid, a skew cell, M=4, and at
   least one two-basis geometry. Keep the oracle independent of UppASD C++.
4. Add an explicit regression that shows the current single-representative
   reciprocal FFT kernel is incomplete: for a 1x1x1 grid the physical FFT bin
   q=0 still contains all nonzero reciprocal aliases.
5. Add energy/field, global-sign-flip, and lattice-translation properties.
6. Correct the Cerdà/Ballenegger/Holm DOI in the design documentation to
   10.1063/1.3657407.
7. Synchronize GPU_FFT_DIPOLE_DESIGN.md with the frozen EWALD3D_FFT name,
   narrow first scope, do_dip=0 rule, and gpu_dipole_tol=1.0e-10 default.

Constraints:
- Preserve the dimensionless B=Hessian(1/r)M convention.
- Tin-foil means only physical k=0 is omitted.
- Do not use CPU UppASD output as expected data.
- Do not loosen tolerances to accommodate an unexplained mismatch.
- Do not edit production GPU code.

Deliver:
- changed files and a concise rationale;
- exact commands and results for all oracle tests;
- a short list of any unresolved physics ambiguity;
- no production-mode activation.
```

### Prompt B — Terra: bridge, ownership, and memory contract

```text
Prepare UppASD's bridge and lifecycle for the first periodic FFT dipole slice,
but keep the physical runtime disabled.

Read completely:
- docs/FFT-dipole_status.md
- docs/FFT-dipole_implementation_plan.md, sections 3, 6, 7, 8 and WP2
- source/CoarseGraining/macrocells.f90
- source/uppasd.f90
- source/chelper.f90
- source/gpu_files/fortranData.{hpp,cpp}
- source/gpu_files/gpuSimulation.cpp
- source/gpu_files/gpuHamiltonianCalculations.cpp
- source/gpu_files/gpuDipoleConvolution.{hpp,cpp}

Implement:
1. Create/export the basis-resolved GPU macrocell layout when the GPU dipole
   mode is requested, independently of legacy do_dip and do_macro_cells.
2. Replace the staged PME3D spelling with EWALD3D_FFT, keep OFF as default,
   and reject gpu_dipole_mode!=OFF together with do_dip!=0.
3. Add a GPU-layout-only cleanup path.
4. Export alat through the C ABI and store a fp64 physical field prefactor
   1e-7*mu_B/alat^3 in the descriptor/owner.
5. Add gpu_dipole_tol throughout the persistent Fortran/C++ bridge with
   default 1.0e-10; reject non-finite/non-positive values.
6. For the first slice, require P/P/P, TINFOIL, block divisibility, and NA=1.
   Reject nonzero legacy alpha/rcut overrides, nonzero gpu_dipole_mesh, and
   GPU MC.
7. Refactor descriptor creation so runtime allocation and memory preflight use
   exactly the same function.
8. Correct memory accounting: both real grids, all spectra, geometry, energy,
   FFT workspace, and peak construction storage with overflow checks.
9. Preserve the final kernel-not-ready fail-fast guard.

Constraints:
- Do not allocate the legacy dense Qdip_macro for GPU use.
- Do not change legacy CPU dipole physics.
- Do not enable non-OFF field evaluation.
- Preserve unrelated working-tree files.

Verify:
- relevant CPU and CUDA/HIP compilation;
- OFF-mode regression;
- invalid input combinations fail before allocation;
- a valid requested mode reaches the intentional kernel-not-ready guard with
  correct printed geometry and budget.

Report changed files, test commands/results, and remaining blockers.
```

### Prompt C — Terra: pure host reference kernel

```text
Implement the first complete periodic Ewald displacement-kernel builder for
UppASD. This package is host-only physics construction; do not wire production
simulation dispatch.

Read:
- docs/FFT-dipole_implementation_plan.md sections 4, 5.1, 8.3, 9 and WP3
- Luna's updated independent oracle/tests
- source/gpu_files/gpuDipoleConvolution.{hpp,cpp}

Implement a pure fp64, CUDA/HIP-independent Builder A that:
1. accepts H, Brec, volume, grid, basis offsets, and gpu_dipole_tol;
2. constructs complete K(d,a,b) in target-minus-source convention;
3. includes all converged real images;
4. includes all nonzero reciprocal vectors, not only FFT representatives;
5. adds +4*alpha^3/(3*sqrt(pi))*I only at d=0,a=b;
6. implements tin-foil with no surface term;
7. returns [cell + Ngrid*kernel_batch] with the documented batch order;
8. automatically chooses internal alpha and real/reciprocal shell extents,
   requiring independently tested real and reciprocal tail budgets for two
   consecutive expansions so cancellation cannot fake convergence;
9. verifies the reference result with a second alpha candidate;
10. reports selected alpha, extents, residual, setup work, and reciprocity
    diagnostics; and
11. performs all construction arithmetic in double even if device real is
   configured as float.

Add CPU-only unit tests against Luna's independent values for cubic,
non-cubic, and skew NA=1 cases. Include alpha invariance and reciprocity.

Do not:
- call buildReciprocalEwaldKernel();
- use GPU results or CPU do_dip as expected values;
- apply the Tesla prefactor inside the dimensionless builder;
- expose alpha/cutoffs as required production input;
- optimize into reciprocal alias grouping yet;
- enable runtime dispatch.

Return the API, data-layout explanation, tests, results, and complexity.
```

### Prompt D — Terra: GPU convolution and normalization

```text
Turn the accepted complete host displacement tensor into a validated GPU FFT
operator, still outside production simulation dispatch.

Read:
- docs/FFT-dipole_implementation_plan.md sections 4.5, 4.7, 8.4, 9.3-9.5
- the accepted Builder A
- gpuDipoleConvolution, gpuFftWrapper, and CMake GPU test conventions

Implement:
1. transient real kernel storage [real_cell,kernel_batch];
2. upload of Builder A's complete real tensor;
3. batched kernel R2C using the owned kernel plan;
4. an explicit spectrum scaling kernel by 1/Ngrid exactly once;
5. release of construction storage after the kernel is ready;
6. a kernel_ready state checked by evaluate();
7. a standalone CMake/CTest GPU test rather than nvcc compilation from Python;
8. diagnostic field/energy readback only in the test API.

First implement a non-physical delta-convolution suite over 1x1x1, 2x1x1,
2x3x1, 3x2x2, basis 1/2, and ensembles 1/4. Only after it passes, run the
periodic Ewald cases.

Constraints:
- raw C2R remains unnormalised; normalization lives in stored Khat.
- preserve the documented field and kernel batch ordering.
- do not retain runtime addRealSpaceField/addPointSelfField composition.
- do not enable gpu_dipole_mode in UppASD.
- test fp64 first.

Run CUDA tests on a GPU worker and at least compile HIP. Report exact maximum
field/energy errors per case, sanitizer results, and all test commands.
```

### Prompt E — Luna: independent GPU acceptance review

```text
Independently review the complete-kernel GPU vertical slice. Do not repair
implementation failures yourself in this task; produce precise failing
evidence for Terra.

Read:
- docs/FFT-dipole_implementation_plan.md sections 4 and 9
- all Builder A and GPU test changes
- periodic_ewald_reference.py

Audit and test:
1. prove from code where the single 1/Ngrid normalization occurs;
2. inspect displacement sign and n1-fast indexing;
3. inspect target/source basis ordering;
4. compare every GPU field component and total energy to the independent
   oracle for the mandatory matrix;
5. run alpha, translation, sign-flip, and M=4 isolation tests;
6. check reciprocity and Hermitian residuals before device upload;
7. run compute-sanitizer;
8. confirm the 1x1x1 q=0 kernel is nonzero and correct.

Acceptance targets are those in section 9.5. Any relaxation needs numerical
convergence evidence and maintainer approval.

Deliver a pass/fail table with exact errors and a prioritized defect list.
Do not approve production wiring unless every mandatory case passes.
```

### Prompt E2 — Terra: WP4b independent GPU acceptance closure

```text
Close WP4b after Luna's independent review. Repair only the standalone
complete-kernel construction and GPU acceptance slice; do not wire dipoles
into production simulation dispatch.

Read:
- docs/FFT-dipole_implementation_plan.md sections 4, 9.5, and WP4b
- Luna's WP4 acceptance report and prioritized defects: docs/luna_acceptance_review.md
- periodic_ewald_reference.py and its independent golden cases
- dipoleEwaldKernel, gpuDipoleConvolution, and test_gpu_dipole_fft.cpp

Implement:
1. Correct the reciprocal-cell cofactor typo in the GPU test helper
   (`h[1]`, not `h[5]`) and add a skew fixture that exposes the old typo.
2. Preserve the existing directConvolution tests as FFT-plumbing tests, but
   add independent versioned oracle fixtures containing geometry, moments,
   dimensionless fields, and total dimensionless energy.
3. Compare every GPU field component and energy directly to independent
   oracle values for 1x1x1, 2x1x1, 2x3x1, a genuinely skew cell, and M=4 in
   the accepted NA=1 fp64 scope. Expected values must not be calculated from
   Builder A's returned/uploaded kernel.
4. Through a test-only explicit-alpha API, build two converged kernels and
   compare their GPU fields and energies.
5. Add GPU translation, global-sign-flip, zero-moment, ensemble-isolation,
   and consecutive changed-moment tests.
6. Add finite and spectral validation before kernel_ready: real-transform
   conjugacy plus Cartesian/basis block Hermiticity derived from
   K(d,a,b)=K(-d,b,a)^T. Report the maximum residual and require <=1e-12 in
   host fp64 diagnostics.
7. Run the full CUDA fp64 CTest matrix and compute-sanitizer.

Constraints:
- Do not remove the useful same-kernel delta/direct-convolution tests.
- Do not use the GPU, Builder A's returned kernel, directConvolution, or
  legacy CPU do_dip to generate independent expected values.
- Do not apply the physical Tesla prefactor in the dimensionless WP4 tests.
- Do not add beff/eneff scatter, production Dip/total-energy integration, or
  remove the intentional production guard; those belong to WP5.
- Do not expand independent physics acceptance to NA>1; that belongs to WP6.
- Normal runtime evaluation must remain free of host synchronization.

Return:
- changed files and a short explanation of each acceptance layer;
- provenance/regeneration instructions for every golden fixture;
- exact test and sanitizer commands;
- maximum field, energy, alpha-invariance, reciprocity, and Hermitian errors;
- a gate table against WP4b.

After implementation, hand the branch/commit and report to Luna to rerun
Prompt E. Tick WP4b and the WP4 umbrella only after Luna independently
approves every WP4b gate.
```

### Prompt F — Terra: production spin-dynamics wiring

```text
Wire the already accepted fp64 NA=1 block-one periodic operator into GPU spin
dynamics. Do not expand scope.

Read:
- docs/FFT-dipole_implementation_plan.md sections 6, 7, 8.5-8.6, 9.6, WP5
- accepted GPU validation report
- gpuHamiltonianCalculations, gpuSDSimulation, gpuSimulation, gpuStructures
- measurement energy kernels and writers

Implement:
1. EWALD3D_FFT initialization from the complete automatically converged
   gpu_dipole_tol kernel and kernel_ready dispatch;
2. macro moment refresh from current emomM at every Hamiltonian evaluation;
3. FFT evaluation and additive Tesla-scaled scatter to beff and required
   eneff path;
4. measurement-time -0.5 sum_macro M.B using the exact same macro field;
5. per-atom convention and Dip/total GPU energy columns;
6. end-to-end fixtures that compare Tesla fields and mRy/atom energy with the
   independent oracle;
7. explicit rejection for MC, NA>1, coarse blocks, non-P/P/P, non-tin-foil,
   fp32, open, and slab modes.

Constraints:
- do_dip must remain 0 for GPU dipoles.
- default gpu_dipole_tol=1.0e-10 must require no alpha/cutoff tuning.
- apply 1e-7*mu_B/alat^3 exactly once.
- do not synchronize the host in the normal field path.
- do not enable any mode beyond the accepted slice.
- preserve OFF-mode behavior.

Verify dipole-only and exchange-plus-dipole totals, M=1/4, changed moments
between evaluations, predictor/corrector stages, non-32-aligned N, sanitizer,
and memory preflight. Report exact output comparisons.
```

### Prompt G — Terra/Luna: multi-basis expansion

```text
Expand the accepted block-one periodic FFT dipole solver from NA=1 to regular
NA>1. Terra implements; Luna owns expected values and acceptance.

Requirements:
- r(g,a)=R(g)+Bas(:,a);
- complete K(d,a,b) with target-minus-source phase;
- existing kernel batch order remains unchanged;
- macro map remains basis-fast;
- scatter decodes every macro channel correctly;
- energy sums every (cell,basis) once;
- basis permutation changes only labels, not physical results.

Test orthogonal and skew NA=2 cells, source impulses in every basis/component,
M=1 and M=4, alpha invariance, reciprocity
K(d,a,b)=K(-d,b,a)^T, and energy derivatives.

Do not begin coarse blocks, true PME deposition, open boundaries, slabs, or
MC in this package. Production NA>1 remains rejected until Luna's full matrix
passes.
```

### Prompt H — Luna then Terra: legacy CPU audit

```text
Audit and repair legacy CPU dipoles only after the GPU periodic solver is
accepted.

Luna first:
- define separate analytic contracts for do_dip=1 finite direct,
  do_dip=2 finite projected macrocell, and do_dip=3 finite zero-padded FFT;
- add non-cubic and NA=2 fixtures;
- demonstrate current FFTW failures without using MKL as the sole oracle;
- inventory all indexing, allocation, normalization, and BC-semantic defects.

Terra second:
- fix the FFTW N3 temporary dimension;
- unify FFTW/MKL moment pack, kernel displacement, and field unpack indexing;
- retain exactly one inverse normalization;
- make unsupported BC semantics explicit;
- pass both backends against the independent finite analytic convolution;
- do not silently change legacy do_dip=3 into periodic Ewald.

If a CPU periodic Ewald FFT mode is desired, implement it as a new mode that
shares the accepted host periodic kernel builder, with its own tests and
documentation.
```

## 12. Definition of done for the first production slice

The first GPU periodic dipole mode is done only when all statements below are
true:

- The algorithm name and boundary convention are unambiguous.
- `gpu_dipole_mode EWALD3D_FFT` is the implemented selector and `OFF` remains
  the default.
- Legacy `do_dip` is not required and cannot double-apply the interaction.
- `gpu_dipole_tol` defaults to `1.0e-10`, automatically selects all internal
  construction parameters, and needs no user alpha/cutoff tuning.
- The GPU-only layout has complete allocation and cleanup.
- `alat` and the Tesla prefactor are wired and tested.
- The complete periodic kernel includes reciprocal aliases.
- The FFT normalization is proven by delta-convolution tests.
- `1x1x1`, non-cubic, skew, and `M=4` fp64 GPU cases match the independent
  oracle.
- Field and energy satisfy `-dE/dM=B`.
- Production `beff`, `eneff`, `Dip`, and total energy are consistent.
- Runtime is enabled only for GPU SD, P/P/P, tin foil, block one, `NA=1`,
  fp64.
- CUDA runtime tests and sanitizer pass.
- HIP at least builds, and HIP runtime is required before claiming HIP
  production support.
- Default `OFF` calculations remain unchanged.
- Open, slab, MC, coarse, fp32, and unsupported surface modes still fail with
  precise messages.

That narrow result is a proper, independently validated FFT-based periodic
dipole implementation. All later extensions should preserve its kernel,
normalization, unit, and energy contracts rather than reopening them.
