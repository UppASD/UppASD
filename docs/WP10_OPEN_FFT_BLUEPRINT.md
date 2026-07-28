# WP10 open-boundary FFT dipole blueprint

## 1. Purpose, authority, and current decision

This document is the execution blueprint for the next GPU dipole work package
after the accepted and optimized three-dimensionally periodic solver.  It
expands the short WP10 entry in `FFT-dipole_implementation_plan.md` into
independently delegable tasks with explicit ownership, prerequisites,
deliverables, tests, and stop conditions.

The work packages are now separated as follows:

```text
WP10  GPU finite/open-boundary dipoles: OPEN_FFT
WP11  Legacy CPU dipole audit and repair
WP12  GPU two-dimensionally periodic dipoles: proposed EWALD2D_FFT
```

WP10 is the only implementation project specified in full detail here.  WP11
and WP12 retain the audit and physics analysis needed to prepare later,
separately budgeted blueprints.

The controlling decision is:

> Implement the existing basis-resolved macrocell representation as a finite
> linear convolution.  There is no second atomistic source representation.
> The reference-resolution case called "atomic" in WP10 means exactly
> `block_size_x=block_size_y=block_size_z=1` and `NA=1`.

For `NA>1`, block size one remains basis-resolved and preserves one macro
channel per basis site, but this blueprint calls it the **basis-resolved
block-one case**, not the atomic baseline.  Every GPU runtime evaluation still
reduces atom moments into the existing macrocell array, transforms macro
channels, and scatters the macro fields back through the existing map.

The first production contract is:

```text
public selector:       gpu_dipole_mode OPEN_FFT
legacy selector:       do_dip 0
boundary flags:        0 0 0 only
execution:             GPU spin dynamics only
source representation: existing basis-resolved macrocell map only
first acceptance gate: block 1 x 1 x 1, NA=1, fp64
physical sum:          finite/open point-dipole Hamiltonian
surface convention:    none
Ewald split:           none
self term:             zero for the exact same point
energy:                -1/2 sum_macro M_macro dot B_macro
field units:           Tesla after one production-boundary prefactor
```

`OPEN_FFT` is not a fallback for a failed `EWALD3D_FFT` initialization.  It is
a separately requested Hamiltonian with different physics.

## 2. Repository baseline

The implementation must extend, not replace, the accepted periodic
infrastructure:

- `source/gpu_files/gpuDipoleConvolution.{hpp,cpp}` owns the regular-grid
  descriptor, R2C/C2R plans, device buffers, spectral block contraction,
  field scatter, energy reduction, workspace, and memory accounting.
- `source/gpu_files/dipoleEwaldKernel.{hpp,cpp}` owns periodic Ewald
  construction.  Open construction must be kept in a separate host component
  so a finite tensor cannot accidentally acquire Ewald terms.
- `source/CoarseGraining/macrocells.f90` owns the basis-resolved macrocell
  layout.  C++ must not infer membership from floating-point coordinates.
- `source/gpu_files/gpuHamiltonianCalculations.{hpp,cpp}` owns GPU Hamiltonian
  selection and runtime dispatch.
- `source/uppasd.f90`, `source/chelper.f90`,
  `source/Input/inputdata.f90`, and `source/Input/inputhandler.f90` own the
  input and Fortran/C bridge.
- `tests/dipole_validation/` owns independent references and acceptance
  drivers.

The current descriptor can already calculate padded extents for an open
boundary, but its runtime is not open-ready:

- `layout.real_cells` currently means both physical source cells and padded
  FFT cells;
- `packMacroMoments()` requires the source macro tensor to contain
  `layout.real_cells`;
- the pack path does not clear padding and embed the smaller active grid;
- scatter and energy iterate using the padded count;
- production descriptor creation and dispatch accept only
  `EWALD3D_FFT`; and
- `buildPeriodicKernel()` correctly refuses an open descriptor.

These are expected staging limitations, not reasons to weaken the open
physics contract.

## 3. Frozen physical and numerical contract

### 3.1 Active macro lattice

Let the active macro grid be:

```text
G = (G1,G2,G3)
Nactive = G1*G2*G3
```

All equations and C++ storage formulas in this blueprint use zero-based
indices.  With `g1` fastest:

```text
active_cell(g) = g1 + G1*(g2 + G2*g3)
macro(g,a)     = a + NA*active_cell(g)
```

The Fortran map stores `1+macro(g,a)`.  `active_cells` counts Bravais
translation cells, not basis channels; the number of active macro channels is
`active_macros = active_cells*NA`.

The runtime field/moment batch remains:

```text
field_batch(alpha,a,e) = alpha + 3*(a + NA*e).

active_moment_index(alpha,g,a,e) =
    alpha + 3*(macro(g,a) + active_macros*e).
```

The first formula is the batch used after packing; the second is the existing
component-fast basis-resolved macro tensor consumed before packing.  Scatter
decodes the one-based Fortran map back to `(active_cell,a)` and reads only the
corresponding active coordinate from the padded field batch.

The first accepted case has block shape `(1,1,1)` and `NA=1`.  A later gate
admits basis-resolved block-one `NA>1`; a still later gate admits uniform
divisible coarse blocks.

### 3.2 Finite convolution

For target cell `g`, target basis `a`, source cell `h`, and source basis `b`,
the dimensionless field is:

```text
d = g-h
B^0_alpha(g,a,e) =
    sum_(h,b,beta) K_alpha_beta(d,a,b) M_beta(h,b,e).
```

For the block-one point operator:

```text
x(g,a)   = g1*C1 + g2*C2 + g3*C3 + Bas(:,a)
r(d,a,b) = x(g,a)-x(h,b)
           = d1*C1 + d2*C2 + d3*C3 + Bas(:,a) - Bas(:,b)

K_alpha_beta(d,a,b) =
    3*r_alpha*r_beta/|r|^5 - delta_alpha_beta/|r|^3.
```

`C1/C2/C3` and `Bas` are Cartesian vectors in lattice-coordinate units.  The
host matrix storage is column-major:
`primitive_vectors[3*axis+component] = C_(axis+1)[component]`.

Only the indexed identity `d=(0,0,0)` together with `a=b` is point self, and
all nine tensor components of that entry are set to exactly `0.0`.  A
different `(d,a,b)` that also produces `r=0` is invalid overlapping geometry
and must be rejected; it is not another self entry to zero.  Conversely,
`d=0,a!=b` is retained whenever the two basis positions are distinct.  There
is no Ewald point-self correction, reciprocal contribution, `k=0` removal,
surface term, or vacuum-box parameter.

The displacement range is:

```text
-(G1-1) <= d1 <= G1-1
-(G2-1) <= d2 <= G2-1
-(G3-1) <= d3 <= G3-1.
```

### 3.3 Padding and embedding

Choose FFT extents `P=(P1,P2,P3)` satisfying:

```text
Pi >= 2*Gi - 1
```

for every axis.  The first implementation should use the exact minimum.
Choosing `2*Gi` or another FFT-friendly extent is a later performance
optimization and is physically equivalent only when the kernel embedding uses
the selected `Pi`.

Active moments occupy:

```text
0 <= i1 < G1
0 <= i2 < G2
0 <= i3 < G3
```

inside an otherwise zero FFT array.  Define:

```text
fft_cell(q) = q1 + P1*(q2 + P2*q3)
q_i(d_i)   = d_i                    for 0 <= d_i <= Gi-1
             P_i + d_i              for -(Gi-1) <= d_i < 0
           = modulo(d_i,P_i).

kernel_batch(alpha,beta,a,b) =
    alpha + 3*(beta + 3*(a + NA*b))

kernel[fft_cell(q(d)) + fft_cells*kernel_batch(alpha,beta,a,b)] =
    K_alpha_beta(d,a,b).
```

Here `alpha` is the target field row and `beta` is the source moment column;
the order is component row fastest, component column, target basis, then
source basis.

The entire padded kernel is zero-initialised before these assignments.  Along
axis `i`, indices `Gi <= qi <= Pi-Gi` are the unused gap (empty when
`Pi=2*Gi-1`) and remain zero.  Since the signed displacement interval has
`2*Gi-1` members, `Pi>=2*Gi-1` makes the mapping injective; the definition is
therefore valid for the minimum or any larger legal `P`.

Padded moment and field buffers use:

```text
real_index(q,alpha,a,e) =
    fft_cell(q) + fft_cells*field_batch(alpha,a,e).
```

The R2C half-spectrum has `S=(P1/2+1,P2,P3)` with integer division and:

```text
spectral_cell(k) = k1 + S1*(k2 + P2*k3)

kernel_spectrum_index(k,alpha,beta,a,b) =
    spectral_cell(k) + spectral_cells*kernel_batch(alpha,beta,a,b)

field_spectrum_index(k,alpha,a,e) =
    spectral_cell(k) + spectral_cells*field_batch(alpha,a,e).
```

Thus `k1` is the retained and fastest R2C axis.  At each retained wavevector
the contraction is:

```text
Bhat(alpha,a,e) =
    sum_(b,beta) Khat(alpha,beta,a,b) Mhat(beta,b,e).
```

The target-minus-source convention must be fixed by impulse tests rather than
inferred from the evenness of the single-basis point tensor.  Multi-basis
tests must detect an accidental basis-phase reversal.

### 3.4 FFT normalization

cuFFT and hipFFT leave the inverse transform unnormalised.  Preserve the
accepted periodic convention:

```text
Nfft                   = P1*P2*P3
stored_kernel_spectrum = R2C(K) / Nfft
runtime C2R output      = B^0, the finite convolution
```

Equivalently, with the libraries' forward-minus/inverse-plus phase pair, raw
C2R of `R2C(K)*R2C(M)` would produce
`Nfft*sum_h K[(g-h) mod P]M[h]`.  Dividing the immutable kernel spectrum once
removes that factor.  No builder-side real-tensor scaling, moment scaling,
runtime C2R scaling, scatter scaling, or energy scaling by `1/Nfft` may also
be added.

### 3.5 Units and energy

The host kernel and `B^0` remain dimensionless.  The numerical macro moments
are in units of `mu_B`.  Apply exactly once:

```text
field_prefactor = 1.0e-7 * mu_B / alat^3
B_Tesla         = field_prefactor * B^0
```

when the device field is added to `beff` and `eneff`.

For each ensemble, the energy reduction uses the same active macro moments
and the same active entries of the unscaled grid field:

```text
E_grid[T*mu_B-number] =
    -0.5 * field_prefactor
    * sum_(active cell,a,alpha) M_alpha(g,a)*B^0_alpha(g,a)

E_dip[J] =
    mu_B * E_grid

E_dip[mRy/atom] =
    (mu_B/mRy) * E_grid / Natom.
```

The first line is the existing pre-conversion numerical quantity, not already
joules.  Padding cells never contribute to energy.  Production measurement
performs the division by `Natom` and the existing `mu_B/mRy` conversion
exactly once, and adds the result to the dipole and total columns exactly
once.  The builder and its cache contain neither `alat` nor the physical
prefactor.

### 3.6 Coarse representation

Uniform divisible coarse blocks use the exact finite projection of the
block-one point operator under the existing macro-moment approximation:

```text
M_A = sum_(i in A) m_i
m_i = M_A/n_A                         for the projected subspace

K_coarse(A,B) =
    1/(n_A*n_B) * sum_(i in A,j in B) K_point(i,j)

B_A = sum_B K_coarse(A,B) M_B
E_coarse = -0.5 * sum_A M_A dot B_A.
```

In `K_point`, only `i=j` is zero.  Therefore the projected diagonal includes
all distinct intra-block atomistic interactions.  Do not add an ad hoc
point-macrospin self term or a Newell finite-volume tensor.  At block one the
formula reduces exactly to the basis-resolved point operator above.

Newell's cuboid tensor describes a different continuum finite-cell model.  It
may be valuable in a future explicitly named solver, but it must not silently
change WP10.

Partial edge blocks have different populations and shapes and are not a
translation-invariant block convolution.  They remain rejected in WP10.

### 3.7 Supported and rejected geometry

The first production slice requires a complete regular Bravais embedding:

```text
s_i = block_size_i > 0
Natom = N1*N2*N3*NA                         (checked arithmetic)
N_i % s_i = 0
G_i = N_i/s_i
pme_macro_grid = G
pme_Num_macro = NA*G1*G2*G3
population(macro(g,a)) = s1*s2*s3
```

Every atom must have exactly one one-based map value in
`[1,pme_Num_macro]`; the observed histogram must equal
`pme_macro_nlistsize`, and its sum must equal `Natom`.  Population equality
alone is not a sufficient regular-grid proof.  For an atom whose integer
lattice identity is `(i1,i2,i3,a)`, the map must equal:

```text
macro = a + NA*(floor(i1/s1)
              + G1*(floor(i2/s2) + G2*floor(i3/s3))).
```

Atom-array order may be geoblocked and is not itself part of the contract.
The Fortran layout owner must construct or attest this identity from integer
lattice membership; C++ must not reverse-engineer it from floating-point
coordinates.  The exported one-based map remains the runtime authority, while
averaged macro centres are diagnostics only.

This distinction matters in the current source:
`create_pme_macrocell_layout()` advances `iatom` in its own PME block-major
walk, while global coordinate creation can use a separate scalar geoblocking
walk and the PME input permits anisotropic `block_size_x/y/z`.  Those walks
must be proven identical for an enabled input or the Fortran builder must form
the map from the atom's integer lattice identity.  Uniform populations cannot
detect a permutation between equal-sized macrocells.  Block one with the
ordinary unit geometry walk is the first accepted case; anisotropic/coarse or
otherwise differently geoblocked orders remain rejected until an
integer-identity map test passes.  Independently, the current Fortran routine
rounds grid counts up with ceiling division, so the GPU mode must apply the
divisibility rejection before FFT allocation rather than accepting partial
edge blocks.

Before allocation the production gate must also require `do_dip=0`,
`BC=0 0 0`, fp64 GPU spin dynamics, finite nonsingular primitive vectors,
finite basis offsets, distinct physical point positions, and the exact scope's
accepted block/basis values.  It must reject MC, mixed or periodic BC,
non-default Ewald/surface overrides, partial edges, unequal populations,
out-of-range or unproven integer-identity maps, occupancy masks/deleted sites, and
coordinate-derived irregular layouts.  Singleton axes, non-cubic grids, and
skew primitive cells are supported regular shapes, not rejection cases.

A geometrical island or racetrack is usable in the first slice only when it is
already represented as a complete regular lattice whose absent magnetic
material has zero moment without violating UppASD's moment and mapping
semantics.  General occupancy masks, deleted sites, partial cells, or
coordinate-derived irregular layouts are later work and must not be inferred
from the mode name.

### 3.8 Hand-worked convolution orientation

First use a deliberately non-even scalar sentinel, because the one-basis
physical point tensor is even and cannot expose a reversal.  Let `G=2`,
`P=4`, `K(-1)=7`, `K(0)=0`, `K(+1)=11`, and active moments `(M0,M1)=(2,3)`.
The padded arrays are:

```text
K[q] = [0,11,0,7]
M[q] = [2, 3,0,0].
```

The active cyclic-convolution outputs are:

```text
B(0) = K(0)M(0) + K(-1)M(1) = 21
B(1) = K(+1)M(0) + K(0)M(1) = 22.
```

Thus `d=g-h` and `q=modulo(d,P)` give the desired linear convolution.  A
source-minus-target embedding would instead give `(33,14)`.  The unused
oversized-padding index `q=2` is explicitly zero.

Now use the physical tensor to settle the basis phase without appealing to
evenness.  Take `G=(2,1,1)`, `P=(4,1,1)`, `C1=(1,0,0)`,
`Bas(0)=(0,0,0)`, and `Bas(1)=(1/4,0,0)`.  Put a unit `x` moment at
`h=0,b=1` and inspect target `g=1,a=0`.  Then:

```text
d = +1
r = C1 + Bas(0)-Bas(1) = (3/4,0,0)
K_xx = 2/(3/4)^3 = 128/27.
```

This value is stored at `q1=1`,
`kernel_batch(x,x,a=0,b=1)=18`, hence at
`kernel[1 + 4*18]`.  Reversing the basis offset in the same batch would use
`r=(5/4,0,0)` and incorrectly produce `128/125`; tensor evenness cannot make
those values equal.  Reciprocity places the matching
`K_xx(d=-1,a=1,b=0)=128/27` at `q1=3`, batch `9`.  If only these two physical
points carry unit `x` moments, their dimensionless interaction energy is
`-0.5*(128/27+128/27)=-128/27`, confirming the field/energy factor of one
half.

## 4. Target architecture

### 4.1 Layout vocabulary

Replace ambiguous cell counts with:

```text
active_grid
active_cells
active_macros
fft_grid
fft_cells
spectral_grid
spectral_cells
field_batches
kernel_batches
```

These are the canonical field names.  `active_grid=G` and
`active_cells=G1*G2*G3` describe Bravais translation cells;
`active_macros=active_cells*NA` includes basis channels.  `fft_grid=P` and
`fft_cells=P1*P2*P3` describe the real FFT allocation.  “Padded cells” means
the `fft_grid` locations outside the embedded active box; there is no separate
`padded_cells` count.  `fft_grid` equals `G` for `EWALD3D_FFT` and a legal
padded `P` for `OPEN_FFT`.  `real_grid` and `real_cells` are deprecated
because they do not distinguish active real data from the real-valued FFT
allocation.  Periodic behavior must remain byte-for-byte or
tolerance-equivalent after the refactor.

`active_macros` is the sole runtime source representation.  An
`atomistic_grid` value retained for later projection describes fine geometry;
it does not authorize a second atom-array FFT path.  In particular, the
atomic baseline still enters through the macro map with unit block and
`NA=1`.

### 4.2 Device path

Initialization:

```text
validate input and regular macro map
-> choose active and FFT layouts
-> construct complete finite host kernel in fp64
-> embed it on P
-> batched R2C
-> normalize spectrum once
-> validate spectrum
-> release construction storage
```

Every field evaluation:

```text
reduce current atom moments to active macro moments
-> clear padded moment grid
-> embed active macro channels
-> batched R2C
-> spectral 3*NA by 3*NA contraction
-> raw C2R
-> scatter active fields additively to atoms
```

Measurement:

```text
reduce -1/2 M_macro dot B_macro over active cells only
-> apply physical prefactor
-> compose the existing dipole and total energy columns once
```

### 4.3 Host builder API

The implementation should introduce a backend-neutral component such as:

```text
source/gpu_files/dipoleOpenKernel.hpp
source/gpu_files/dipoleOpenKernel.cpp
```

with value types similar to:

```cpp
struct DipoleOpenGeometry {
    std::array<std::size_t,3> active_grid{};
    std::array<std::size_t,3> fft_grid{};
    std::array<std::size_t,3> block_shape{{1,1,1}};
    // Column-major [C1 C2 C3]:
    // primitive_vectors[3*axis+component].
    std::array<double,9> primitive_vectors{};
    unsigned int basis = 0;
    std::vector<std::array<double,3>> basis_offsets;
};

struct DipoleOpenKernelDiagnostics {
    std::array<std::size_t,3> active_grid{};
    std::array<std::size_t,3> fft_grid{};
    std::size_t active_cells = 0;
    std::size_t fft_cells = 0;
    std::size_t kernel_batches = 0;
    std::size_t nonfinite_values = 0;
    bool all_finite = false;
    double minimum_nonself_r2 = 0.0;
    double max_reciprocity_error = 0.0;
    double max_point_self_abs = 0.0;
    double max_padding_gap_abs = 0.0;
};

struct DipoleOpenKernelResult {
    std::vector<double> kernel;
    DipoleOpenKernelDiagnostics diagnostics;
};

DipoleOpenKernelResult buildOpenDipoleDisplacementKernel(
    const DipoleOpenGeometry& geometry);
```

The public names may be adjusted mechanically during implementation, but
these value semantics and types are frozen.  In particular there is no
settings/tolerance object.  The builder must:

- use host fp64 even in an fp32 device build;
- have no CUDA/HIP dependency;
- have no Ewald settings;
- reject zero dimensions, overflow, `Pi<2*Gi-1`, non-finite/singular
  primitive geometry (`det([C1 C2 C3])` must be finite and positive),
  basis-count/offset mismatch, and any nonself `r=0`;
- initially require `block_shape=(1,1,1)` and reject other values until the
  separately derived WP10.7 projection is supplied;
- return `fft_cells*9*basis*basis` doubles in the storage order of section
  3.3;
- report the active and FFT dimensions/counts, minimum nonself distance,
  non-finite count, reciprocity error, exact-self error, and unused-gap error;
- validate finite values, reciprocity, the exact self rule, and zero padding;
  and
- be independently testable without creating a GPU context.

Contract violations must fail without returning a partial tensor:
`std::invalid_argument` for invalid geometry/model fields and
`std::overflow_error` for an unrepresentable extent, batch, element, or byte
count.  The diagnostic message must name the failed field and its relevant
values.  A successful result has `all_finite=true`,
`nonfinite_values=max_point_self_abs=max_padding_gap_abs=0`, and a strictly
positive `minimum_nonself_r2` whenever a nonself pair exists.

### 4.4 Cache and memory

Open padding can approach eight times the active real-grid volume in a fully
three-dimensional sample.  The preflight must include:

- padded moment and field real buffers;
- padded R2C moment and field spectra;
- padded kernel spectrum;
- shared FFT workspace;
- construction-only real kernel storage;
- any transient host-to-device conversion buffer; and
- active macro moments and map storage.

The minimum immutable cache key is:

- `OPEN_FFT` plus an exact open point/projection model and artifact-format
  revision;
- cached artifact/storage precision (`host-fp64` for the builder result;
  device real/complex precision as well if a converted spectrum is cached);
- active `G` and selected FFT `P`;
- block shape;
- all nine primitive-vector doubles;
- ordered basis count and all ordered basis-offset doubles.

Integral fields use their exact values and floating geometry uses the exact
post-staging IEEE-754 bit patterns, not rounded text.  Ensemble count,
`Natom`, `alat`, field prefactor, moment values, `do_dip`, and BC do not affect
the immutable dimensionless kernel and are not key fields.  FFT backend/device
identity is needed only if the cached artifact format or backend ABI is
backend-specific; it is not part of a backend-neutral host-double tensor key.
No Ewald alpha, cutoff, mesh, tolerance, surface value, or vacuum gap belongs
in an open-kernel cache key.

## 5. Global validation matrix

Every accepted implementation must cover:

| Case | Purpose |
|---|---|
| `1x1x1`, `NA=1`, block one | exact zero self field and energy |
| `2x1x1`, axial moment | tensor sign and factor two |
| `2x1x1`, transverse moment | transverse sign |
| `2x3x1` | non-cubic active/padded indexing |
| a case with `G3>G2` | expose the legacy FFTW-style dimension error |
| skew `C1/C2/C3` | Cartesian displacement construction |
| block one, `NA=2` | basis phases and channel order |
| `M=4` | batch/ensemble independence |
| nonuniform moments | prevent symmetry from hiding permutations |
| uniform thin film, in-plane/out-of-plane | finite-sample shape response |
| exchange plus `OPEN_FFT` | additive field and energy composition |

Required mathematical checks:

- same-kernel direct convolution versus FFT application;
- independent finite pair sum versus host builder;
- independent finite pair sum versus GPU field;
- reciprocity
  `K(d,a,b)=K(-d,b,a)^T`;
- R2C conjugacy and block Hermiticity;
- delta sources in every component, basis, corner, and ensemble;
- no wraparound between opposite sample faces;
- energy equals `-1/2 sum M dot B`;
- finite-difference energy derivatives reproduce the field;
- changing unused padding values from uninitialised memory is impossible
  because padding is explicitly zeroed; and
- repeated evaluation with changed moments does not reuse stale padding or
  macro moments.

Suggested fp64 budgets:

```text
host direct/operator tests:  relative 1e-12, absolute case-specific
GPU internal diagnostics:    relative 1e-11 to 1e-12
Fortran text E2E:             limited by printed precision, normally ~2e-8
energy derivative:            step-size convergence must be documented
```

CUDA and HIP must satisfy the same physics budget; bitwise equality is not
required.

## 6. Delegation rules

Each prompt below is intended for a separate session.  Give an agent only one
prompt at a time.  The prompt contains its own scope and required reading.

All sessions must:

1. Run `git status --short` before editing.
2. Preserve unrelated tracked and untracked files.
3. Use one branch/commit or one clearly identified working-tree change set per
   prompt.
4. Record the commit, build commands, test commands, and numerical results
   under the corresponding checkbox.
5. Stop when a prerequisite gate is not satisfied rather than implementing a
   different physical mode.
6. Never use CPU output or another GPU builder as the sole expected value.

File ownership is deliberately separated:

```text
Sol:    physical contract, host open builder, coarse projection
Terra:  descriptor/layout, device runtime, input/bridge/dispatch, performance
Luna:   independent oracle, red tests, acceptance and audit reports
```

Sol and Luna may work concurrently on Tasks 10.0 and 10.1.  Terra may begin
Task 10.2 after the layout names in Task 10.0 are frozen.  Sol's Task 10.3 and
Terra's Task 10.2 may then proceed concurrently because their primary files do
not overlap.  Task 10.7 is deliberately a three-way handoff: Sol owns the
coarse projection, Terra owns production integration, and Luna owns the
independent acceptance decision.

Recommended session order:

```text
WP10.0 Sol ─┬─> WP10.2 Terra ─┐
            │                 ├─> WP10.4 Terra -> WP10.5 Luna
WP10.1 Luna ┴─> WP10.3 Sol ───┘

WP10.5 GO -> WP10.6a Terra -> WP10.6b Luna
WP10.6b GO -> WP10.7a Sol -> WP10.7b Terra -> WP10.7c Luna
last physical-scope GO -> WP10.8 Terra -> WP10.9 Luna
```

## 7. Explicit WP10 task checklist and prompts

### [x] WP10.0 — Sol — freeze the open operator and host interface

Dependencies: accepted periodic WP9; this blueprint.

Deliverables:

- a short contract confirmation or corrections committed to this document;
- the final active/padded layout vocabulary;
- the target host-builder API and kernel storage convention;
- explicit decisions for open input validation and cache keys; and
- no production dispatch.

Gate:

- Luna and Terra can quote one unambiguous displacement, embedding,
  normalization, unit, energy, and self convention.

Delegation prompt:

```text
You are Sol.  Freeze the physics and host-side interface for UppASD GPU WP10
OPEN_FFT.  This is a design/contract task, not production enablement.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md in full;
- docs/FFT-dipole_implementation_plan.md sections 3.3, WP8, WP10, and WP11;
- docs/GPU_FFT_DIPOLE_DESIGN.md physical contracts;
- source/gpu_files/gpuDipoleConvolution.hpp;
- source/gpu_files/dipoleEwaldKernel.hpp;
- source/CoarseGraining/macrocells.f90 create_pme_macrocell_layout().

Frozen intent:
- There is only the existing basis-resolved macrocell source representation.
- "Atomic baseline" means block_size=(1,1,1) and NA=1.
- NA>1 block one is called basis-resolved block one.
- OPEN_FFT is the exact finite point-dipole Hamiltonian at block one.
- do_dip must be zero; BC must be 0 0 0; first scope is fp64 GPU SD.
- No Ewald, surface, vacuum-gap, or continuum finite-cell term is allowed.

Audit every equation and index convention in the blueprint.  Resolve:
1. target-minus-source displacement and basis offset convention;
2. embedding of d in an arbitrary P >= 2G-1;
3. kernel batch order;
4. sole inverse normalization;
5. physical prefactor and energy;
6. exact point self;
7. active versus padded layout vocabulary;
8. regular-grid, population, divisibility, and unsupported-shape gates;
9. the host-fp64 builder value types and diagnostics; and
10. the minimum complete cache key.

Use small hand-worked 1D and two-basis examples to verify the convolution
orientation.  Do not rely on the tensor being even to settle the multi-basis
phase.  If this blueprint is wrong or ambiguous, update only the relevant
design text and explain the correction.  Do not edit GPU dispatch, Fortran
input, production kernels, or Luna's independent reference.

Finish with:
- files changed;
- frozen equations and storage order;
- unresolved decisions, if any;
- commands/checks run; and
- an explicit GO/NO-GO for WP10.1, WP10.2, and WP10.3.
```

Acceptance record:

```text
Assignee: Sol
Branch/commit: working-tree design change; no commit requested
Checks: full required-source audit; hand sentinel and two-basis arithmetic;
        Fortran/PME ordering audit; Markdown fence and whitespace checks
Decision: WP10.0 GO; WP10.1 GO; WP10.2 GO;
          WP10.3 NO-GO until WP10.1 supplies its independent oracle/goldens
```

### [x] WP10.1 — Luna — independent finite oracle and red-test package

Dependencies: WP10.0 equations may be reviewed concurrently, but any
disagreement must be resolved before expected values are frozen.

Deliverables:

- a backend-neutral finite direct oracle;
- deterministic golden cases;
- active/padded impulse tests that fail on the current implementation; and
- a concise red-test report.

Gate:

- the tests distinguish finite convolution from periodic wraparound and catch
  every known active/padded indexing failure before Terra implements the fix.

Delegation prompt:

```text
You are Luna.  Own the independent oracle and red tests for WP10 OPEN_FFT.
Do not implement production GPU construction or dispatch.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md sections 1-5;
- tests/dipole_validation/README.md;
- tests/dipole_validation/check_fields.py;
- tests/dipole_validation/generate_cases.py;
- tests/dipole_validation/test_gpu_dipole_fft.cpp;
- tests/dipole_validation/test_host_periodic_ewald_kernel.cpp only for test
  harness style, not expected values;
- source/Hamiltonian/dipolecommon.f90 only to understand the legacy sign and
  units, not as an oracle.

Create an independent fp64 finite point-dipole evaluator:
- positions are generated directly from C1/C2/C3 and basis offsets;
- sum every distinct finite source-target pair exactly once;
- use the analytic point tensor;
- exclude only exact self;
- return per-site fields, total energy, and energy per atom;
- apply units in a separate layer from the dimensionless tensor;
- never call UppASD, the future open builder, periodic Ewald code, do_dip=3,
  or GPU code to generate expected values.

Add deterministic cases for:
- 1x1x1;
- axial and transverse 2x1x1 pairs;
- 2x3x1 with nonuniform moments;
- a case with G3>G2;
- one skew cell;
- NA=2 with unequal basis offsets and nonuniform basis moments;
- M=4 with distinct ensembles;
- a finite thin film magnetized in plane and normal to the film.

Add red tests for the GPU layout test API:
- active macro storage has Nactive, not Nfft, entries;
- padded buffers are explicitly zero outside the embedded active box;
- delta sources in every active corner do not wrap to the opposite face;
- every component/basis/ensemble uses the documented batch;
- direct same-kernel convolution agrees with expected values before any
  physical prefactor;
- energy and a finite-difference derivative agree.

It is acceptable for GPU-facing red tests not to compile until a small test
seam from WP10.2 exists.  In that case, land the oracle and host tests first
and document the exact test seam Terra must provide.  Do not weaken a test to
make the current periodic-only class pass.

Finish with:
- files changed;
- independent derivation and constants used;
- golden cases and why each is non-redundant;
- expected failing tests and failure signatures;
- commands run; and
- an explicit statement that no expected value came from CPU or GPU output.
```

Acceptance record:

```text
Assignee: Luna
Branch/commit: WP10.1 independent oracle and red-test package
Red tests: host oracle/embedding tests pass; GPU layout red source compiles and
           is intentionally pending Terra's `luna_open_fft_test::run` seam
Failure signatures: current GPU red link fails with undefined
                   `luna_open_fft_test::run(...)`; runtime gates cover
                   active/padded counts, zero padding, no-wrap corners, batch
                   order, dimensionless convolution, energy, and derivative
Oracle independence review: GO; expected values come only from the finite
                            analytic point-dipole evaluator
```

### [x] WP10.2 — Terra — active/padded layout and device embedding

Dependencies: WP10.0 layout vocabulary; Luna's requested test seam from
WP10.1 when available.

Deliverables:

- explicit active and FFT grid/count fields;
- zero-and-embed pack path;
- active-only scatter and energy indexing;
- updated memory/workspace estimates; and
- no production `OPEN_FFT` input enablement.

Gate:

- all periodic tests remain accepted;
- nonphysical open delta convolutions pass for NA=1/2 and M=1/4.

Delegation prompt:

```text
You are Terra.  Refactor GpuDipoleConvolution so active source cells and padded
FFT cells are distinct.  This is an FFT plumbing task.  Do not enable the
OPEN_FFT user mode and do not invent its physical kernel.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md sections 1-5;
- the accepted output of WP10.0;
- Luna's WP10.1 test-seam request, if present;
- source/gpu_files/gpuDipoleConvolution.{hpp,cpp};
- source/gpu_files/gpuSimulation.cpp memory preflight;
- tests/dipole_validation/test_gpu_dipole_fft.cpp.

Required changes:
1. Represent active_grid/active_cells separately from fft_grid/fft_cells.
2. Preserve periodic layout behavior where active and FFT grids are equal.
3. Size macro_moments against active_cells*basis, not fft_cells*basis.
4. Before every open forward transform, zero every padded moment buffer entry.
5. Embed active cells using n1-fast active and padded strides.
6. Keep channel/ensemble batches unchanged.
7. Make scatter decode the existing one-based basis-fast macro map but fetch
   fields using padded strides and active cell coordinates.
8. Make diagnostic and production energy iterate active cells only.
9. Update all shape validation, overflow checks, persistent bytes,
   construction bytes, and workspace estimates.
10. Provide a test-only way to upload an arbitrary padded real kernel and run
    open-layout delta convolutions without making OPEN_FFT a user mode.

Do not:
- change the periodic Ewald builder;
- change the physical prefactor;
- introduce an atomistic device array;
- infer macrocells from coordinates;
- accept partial edge blocks;
- use stale or uninitialised padding;
- add a runtime inverse normalization; or
- enable Fortran input.

Test at minimum:
- periodic existing 1x1x1, non-cubic, skew, NA=2, and M=4;
- open-layout G=2x3x1 and G=2x1x3;
- P=2G-1 and one deliberately larger legal P;
- source impulses in every corner/component/basis/ensemble;
- repeated evaluation after changing all source values;
- memory estimates equal the allocation inventory.

Finish with:
- files changed;
- before/after layout table;
- test API added;
- commands and numerical maxima;
- periodic regression result;
- sanitizer result if a GPU is available; and
- confirmation that production OPEN_FFT is still unreachable.
```

Acceptance record:

```text
Assignee: Terra
Branch/commit: active/padded FFT plumbing commit
CUDA tests: fp64 focused CTest 4/4 passed; open-layout impulse matrix
            field max 1.1368683772161603e-13, energy max 4.5796699765787707e-16
HIP tests: not run (CUDA-only acceptance host)
Sanitizer: compute-sanitizer memcheck ERROR SUMMARY: 0 errors
Periodic regression: existing 1x1x1, non-cubic, skew, NA=2, and M=4 suites passed
Production OPEN_FFT: still unreachable; no selector, Fortran input, or host
                     open-kernel builder was added
```

### [x] WP10.3 — Sol — host-fp64 finite open-kernel builder

Dependencies: WP10.0 GO; Luna's oracle cases.  May proceed in parallel with
WP10.2 if the frozen host storage contract is not changed.

Deliverables:

- `dipoleOpenKernel` host component;
- construction diagnostics;
- host tests against Luna's oracle;
- no GPU production dispatch.

Gate:

- the builder passes every dimensionless host field, reciprocity, basis phase,
  self, and energy check independently of FFT libraries.

Delegation prompt:

```text
You are Sol.  Implement the backend-neutral host-fp64 finite open dipole
kernel for WP10.  Do not enable production dispatch.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md sections 1-5;
- the frozen WP10.0 decision;
- Luna's independent oracle/golden cases from WP10.1;
- source/gpu_files/dipoleEwaldKernel.{hpp,cpp} for value-type and CMake style
  only;
- source/gpu_files/gpuDipoleConvolution.hpp for the accepted kernel batch;
- source/Hamiltonian/dipolecommon.f90 only as historical context.

Implement a separate dipoleOpenKernel host component that:
- accepts active G, selected padded P, C1/C2/C3, basis offsets, basis count,
  and later a uniform block shape;
- validates Pi >= 2*Gi-1 and all arithmetic;
- constructs the complete padded displacement tensor in fp64;
- uses target-minus-source displacement consistently;
- sets only exact same-point/same-basis self to zero;
- returns layout [fft_cell + fft_cells*kernel_batch];
- contains no Ewald alpha, tolerance, reciprocal aliases, surface term, or
  physical-unit prefactor;
- reports all section 4.3 diagnostics, including max reciprocity,
  exact-point-self, unused-gap, minimum-nonself-distance, and finite-value
  status;
- supports skew primitive cells and NA>1 from the start.

First implement and accept block shape 1x1x1.  If a coarse API field is
included, reject non-unit blocks until WP10.7; do not guess the projection.

Tests must compare fields formed from the built tensor with Luna's independent
finite oracle for all WP10.1 cases.  Add direct convolution tests that exercise
all components and basis pairs.  Include a hand-worked two-basis phase test
which would fail if Bas(a)-Bas(b) were reversed.  Test an arbitrary legal
padding larger than 2G-1.

Do not use periodic Ewald output, CPU do_dip=3, or GPU FFT output as expected
values.  Do not copy the periodic cache or tolerance semantics.

Finish with:
- files changed and public API;
- equations actually implemented;
- diagnostics;
- host test commands and maximum errors;
- build integration; and
- confirmation that OPEN_FFT remains production-disabled.
```

Acceptance record:

```text
Assignee: Sol
Branch/commit: this WP10.3 implementation commit
Host tests: dipole-open-fft-oracle, dipole-open-host-builder,
            dipole-open-host-goldens, and dipole-ewald-host-builder passed;
            strict-warning and UBSan host runs passed
Maximum field error: 5.3290705182007514e-15 versus Luna's goldens;
                     7.8159700933611020e-14 for the complete direct impulse matrix
Maximum energy error: 3.5527136788005009e-15 dimensionless
Maximum reciprocity error: 7.1054273576010019e-15
Production OPEN_FFT: still unreachable; no selector, Fortran input, C bridge,
                     GPU upload, or runtime dispatch was added
```

### [x] WP10.4 — Terra — input, lifecycle, kernel upload, and production dispatch

Dependencies: WP10.1 oracle, WP10.2 plumbing, WP10.3 host builder.

Deliverables:

- `gpu_dipole_mode OPEN_FFT`;
- full input and bridge validation;
- host builder upload and production runtime;
- additive field and energy integration;
- startup diagnostics and memory reporting.

Gate:

- first production slice runs with `do_dip=0`, block one, `NA=1`, fp64 and
  all-open BC, while every unsupported combination fails before allocation.

Delegation prompt:

```text
You are Terra.  Integrate the accepted WP10 open builder and active/padded FFT
plumbing into the first production OPEN_FFT mode.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md in full;
- WP10.0 decision;
- WP10.1 oracle/test report;
- WP10.2 and WP10.3 implementation reports;
- source/Input/inputdata.f90 and inputhandler.f90;
- source/uppasd.f90 GPU dipole validation;
- source/chelper.f90 and source/gpu_files/fortranData.{hpp,cpp};
- source/gpu_files/gpuSimulation.cpp;
- source/gpu_files/gpuHamiltonianCalculations.{hpp,cpp};
- source/gpu_files/gpuDipoleConvolution.{hpp,cpp}.

Implement:
1. Parse and export OPEN_FFT as a distinct mode ID.
2. Keep OFF and EWALD3D_FFT semantics unchanged.
3. Require do_dip=0 and BC=0 0 0.
4. Require GPU spin dynamics; reject GPU MC and other unsupported drivers.
5. Reject Ewald alpha/rcut/mesh overrides and non-default surface requests
   because they have no meaning for OPEN_FFT.
6. Build/export the existing basis-resolved macrocell map based on the GPU
   mode, independently of legacy do_dip and do_macro_cells.
7. First require block=(1,1,1), NA=1, complete regular population, and fp64.
8. Build the host open tensor once, upload/transform/normalize it once, then
   release construction storage.
9. At each field call reduce current atom moments, clear/embed the active grid,
   execute the FFT operator, and add to beff and eneff.
10. On measurement steps accumulate the exact matching active-macro energy
    into the established dipole and total columns once.
11. Include padded extents, half-spectrum size, persistent bytes, peak bytes,
    block, basis, ensembles, and prefactor in startup diagnostics.
12. Keep release idempotent and plan/workspace destruction ordered.

Failure checks must occur before device allocation for:
- wrong BC pattern;
- do_dip conflict;
- invalid or partial layout;
- non-unit block in this first gate;
- NA != 1 in this first gate;
- fp32 before its later acceptance;
- MC;
- overflow or unsafe memory;
- any attempt to treat OPEN_FFT as a periodic fallback.

Run host, GPU convolution, and E2E tests.  Verify exchange plus OPEN_FFT is
additive and that OFF/periodic cases are unchanged.

Finish with:
- files changed;
- complete input matrix;
- startup output example;
- build/test/sanitizer commands;
- field and energy errors;
- memory figures; and
- explicit list of still-rejected WP10 extensions.
```

Acceptance record:

```text
Assignee: Terra
Branch/commit: recorded with the WP10.4 implementation change
Input rejection matrix: OPEN_FFT accepts only do_dip=0, BC=0 0 0, GPU SD,
                        fp64, block=1x1x1, NA=1, default Ewald/surface inputs;
                        rejects PBC/mixed BC, MC, fp32, partial/non-unit layouts,
                        NA>1, legacy dipoles, and Ewald/surface overrides.
CUDA E2E: focused convolution/layout/WP5 periodic regressions passed. A real
           2x1x1 OPEN_FFT SD run gave the finite axial pair field
           (6.86963703e-32, 6.86963703e-33, -1.37392741e-32) T; OFF-to-OPEN
           exchange coexistence gave the same field delta and one dipole-only
           addition to total energy.
HIP build/run: unavailable: neither `hipcc` nor `rocminfo` is present and no
               AMD device is available on this runner.
Sanitizer: final CUDA Compute Sanitizer memchecks passed for both
           `dipole_open_fft_layout_tests` and a valid 2x1x1 production
           `sd.f95.cuda` OPEN_FFT SD launch: `ERROR SUMMARY: 0 errors`.
Memory re-address (CUDA 13.3, RTX A4000, fp64): the preflight now covers the
        complete run rather than only base matrices plus dipole FFT buffers.
        The 2x1x1 total is 576 B base matrices + 23,184 B GPU measurement
        phase + 256 B integrator/thermfield + 1,040 B OPEN_FFT phase =
        25,056 B. The OPEN_FFT phase is 776 B persistent + 216 B finite
        real-kernel construction + 0 B FFT workspace + 48 B macro moments.
        The prior 504 B construction figure incorrectly included a periodic
        reciprocal-alias staging buffer that OPEN_FFT never allocates.
```

WP10.4 re-address record (2026-07-28):

```text
Build: cmake --build build_gpu_wp9_fp64 --target sd.f95.cuda \
       dipole_open_fft_layout_tests dipole_gpu_fft_tests -j1
Functional regressions:
  ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
    -R 'dipole-open-fft-oracle|dipole-open-host-builder|dipole-open-host-goldens|dipole-ewald-host-builder'
  PASS 4/4
  ctest --test-dir build_gpu_wp9_fp64 --output-on-failure \
    -R 'dipole-gpu-fft-convolution|dipole-open-fft-layout|dipole-gpu-wp5-e2e'
  PASS 3/3
Memcheck:
  /usr/local/cuda-13.3/bin/compute-sanitizer --tool memcheck --error-exitcode=99 \
    build_gpu_wp9_fp64/bin/dipole_open_fft_layout_tests
  ERROR SUMMARY: 0 errors
  OMP_NUM_THREADS=1 /usr/local/cuda-13.3/bin/compute-sanitizer --tool memcheck \
    --error-exitcode=99 build_gpu_wp9_fp64/bin/sd.f95.cuda  # valid OPEN_FFT SD fixture
  ERROR SUMMARY: 0 errors
```

| CUDA production input | Projected/full-run B | Persistent B | Construction B | Workspace B | OPEN_FFT phase peak B | Tracker peak B | Release B | Comparison |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `2x1x1`, `NA=M=1`, padded `3x1x1` | 25,056 | 776 | 216 | 0 | 1,040 | 25,056 | 0 | +0.0% |
| `2x2x1`, `NA=1`, `M=4`, padded `3x3x1` | 34,400 | 5,048 | 648 | 0 | 6,080 | 34,400 | 0 | +0.0% |
| thin film `4x3x1`, `NA=M=1`, padded `7x5x1` | 37,336 | 6,632 | 2,520 | 0 | 9,440 | 37,336 | 0 | +0.0% |

The tracker is a process-lifetime peak; the preflight is now explicitly the
same full-run allocation scope, and the release line is the live remaining
allocation after phase-plan/workspace destruction and base-matrix release.
The 25,056 B tracker value is therefore expected, not an OPEN_FFT-only
workspace figure.  Releasing the thermfield's owned field and sigma tensors
also makes the final live inventory exactly zero and keeps repeated release
idempotent.

### [ ] WP10.5 — Luna — independent acceptance of the first production slice

Dependencies: WP10.4.

Deliverables:

- an acceptance report independent of the implementation;
- numerical field, energy, derivative, coexistence, rejection, and memory
  results;
- GO/NO-GO for enabling the first public slice.

Gate:

- no known physics, indexing, lifecycle, dispatch, or reporting defect remains
  in block-one `NA=1` fp64 GPU SD.

Delegation prompt:

```text
You are Luna.  Perform independent acceptance of the first production
OPEN_FFT slice.  Do not repair implementation failures in the acceptance
change set; report them and return the package to its owner.

Read:
- docs/WP10_OPEN_FFT_BLUEPRINT.md sections 1-7;
- all WP10.0-WP10.4 handoff reports;
- your independent oracle and red tests;
- the production input/dispatch diff only to identify test surfaces.

Acceptance matrix:
- 1x1x1 zero self;
- axial and transverse pair;
- non-cubic 2x3x1 with nonuniform moments;
- G3>G2;
- skew cell;
- distinct M=4 ensembles;
- finite thin film in-plane and normal magnetization;
- exchange plus dipole;
- OFF regression;
- EWALD3D_FFT regression.

Check:
1. every field component against the independent direct finite oracle;
2. dipole and total energy against -1/2 sum M dot B;
3. finite-difference energy derivatives;
4. no opposite-face wraparound;
5. changed moments refresh both macro moments and zero padding;
6. beff and eneff additive composition;
7. no double application through do_dip;
8. startup mode/geometry/memory diagnostics;
9. rejection before allocation for PBC, mixed BC, MC, do_dip conflict,
   non-unit block, NA>1, fp32, partial blocks, and irrelevant Ewald inputs;
10. estimated versus actual allocation inventory;
11. CUDA memcheck and HIP numerical execution when available.

Use internal diagnostic precision for the tight fp64 gate.  Treat formatted
Fortran text as a separate lower-precision E2E layer.  CPU do_dip=1 may be
reported as an additional regression comparison, but your analytic oracle is
the authority.  Do not use do_dip=3 as expected data.

Write docs/luna_wp10_open_acceptance.md with commands, hardware/backend,
maximum errors, memory, rejection results, and any limitations.  Mark WP10.5
complete only with an explicit GO.  A NO-GO report must identify the owning
task and leave production enablement gated.
```

Acceptance record:

```text
Assignee:
Branch/commit:
Report:
Decision:
```

WP10.4 re-addressing required after Luna's 2026-07-28 independent
acceptance (NO-GO):

- [x] Reconcile the full-run GPU memory budget with the live allocation
      inventory. In the valid CUDA 2x1x1 `NA=M=1` fp64 run, preflight
      projected `2.2 kB`, the OPEN_FFT-specific diagnostic reported
      `0.001 MiB` persistent and `0.001 MiB` peak, but the allocator tracker
      reported `Max GPU 25.056 KB` and release printed `estimate -91.4%`.
      The 2x2x1 M=4 and 4x3x1 thin-film runs likewise printed `-64.9%` and
      `-54.4%`. Determine and document whether these are different scopes
      (full simulation versus dipole-only workspace); if so, print the scopes
      and like-for-like figures, otherwise correct the estimate or inventory.
- [x] Re-run valid CUDA production OPEN_FFT startup cases after the
      memory-accounting change and retain projected bytes, persistent and
      peak bytes, live tracker peak, and release comparison in the WP10.4
      handoff. Do not treat formatted `0.000 GB` output as sufficient
      evidence; retain byte-level values.
- [x] Re-run CUDA memcheck on both the standalone layout seam and an actual
      production `gpu_dipole_mode OPEN_FFT` SD launch, recording the complete
      `ERROR SUMMARY` lines. Luna found zero sanitizer errors before this
      re-address, so this is a regression guard rather than a reported CUDA
      memory fault.
- [x] Preserve the already-passing functional scope: independent direct-
      oracle fields/energies, M=4 isolation, thin films, exchange additivity,
      OFF/EWALD3D_FFT regressions, rejection gates, and no double application
      through `do_dip`. No new physics mode or `do_dip=3` expected data is
      needed for this re-address.
- [x] If an AMD runner becomes available, run the equivalent HIP numerical
      production check; the current NO-GO is not based on a HIP numerical
      failure because no HIP device/toolchain was available.

Terra must attach the updated commands, hardware/backend, byte-level memory
table, and pass/fail interpretation to the WP10.4 handoff. WP10.5 remains
unchecked until Luna independently reruns the affected acceptance and gives
an explicit GO.

### WP10.6 — Basis-resolved block-one `NA>1`

Dependencies: WP10.5 GO.

Deliverables:

- remove only the `NA=1` production gate;
- basis-resolved GPU field and energy acceptance;
- preserve the single macrocell source representation.

Gate:

- every basis/component/ensemble passes the independent oracle and basis
  permutation tests.

#### [ ] WP10.6a — Terra — implement basis-resolved block-one `NA>1`

Delegation prompt:

```text
You are Terra.  Expand accepted OPEN_FFT from the NA=1 atomic baseline to
basis-resolved block-one NA>1.  This is not a new atomistic backend.

Read docs/WP10_OPEN_FFT_BLUEPRINT.md, the WP10.5 acceptance report, and the
accepted multi-basis periodic implementation/tests.

Retain:
- macro=basis+NA*cell;
- kernel batch row+3*(column+3*(target_basis+NA*source_basis));
- target-minus-source basis phase;
- active/padded embedding;
- one normalization and one physical prefactor;
- active-macro energy.

Remove only the NA=1 rejection after adding E2E coverage.  Do not begin coarse
blocks, partial edges, fp32, MC, or slabs.  Test NA=2 in orthogonal and skew
cells, source impulses in every basis/component, M=1/4, nonuniform moments,
basis permutations, reciprocity, and energy derivatives.

Hand the result to Luna without changing Luna's expected values.
```

#### [ ] WP10.6b — Luna — accept basis-resolved block-one `NA>1`

Delegation prompt:

```text
You are Luna.  Independently accept WP10 basis-resolved block-one NA>1.
Use your finite direct oracle, not implementation output.  Verify every basis
phase, source/target pair, component, ensemble, basis permutation, field,
energy, and derivative on orthogonal and skew NA=2 cases.  Confirm the runtime
still uses the existing macrocell representation and that NA=1 is unchanged.
Return GO/NO-GO with maximum errors and sanitizer/backend results.
```

Acceptance record:

```text
Terra branch/commit:
Luna report:
Decision:
```

### WP10.7 — Uniform divisible coarse blocks

Dependencies: WP10.6b GO; accepted periodic WP8 projection conventions.

Deliverables:

- exact finite projection builder;
- production coarse-block enablement for uniform divisible blocks;
- convergence and rejection report.

Gate:

- projected direct and FFT operators agree, block one is recovered exactly,
  and partial edges remain rejected.

#### [ ] WP10.7a — Sol — implement the finite coarse projection

Delegation prompt:

```text
You are Sol.  Extend the open host builder to uniform divisible coarse blocks
using the exact finite projection of the block-one point-dipole Hamiltonian.

Read docs/WP10_OPEN_FFT_BLUEPRINT.md section 3.6, the accepted periodic WP8
projection implementation/tests, and Luna's WP10.6b report.

Derive and implement the block tensor under M_block=sum atom moments.  Include
the finite diagonal produced by intra-block atom pairs.  Do not add a point
macrospin self-demag term and do not use Newell's continuum-cell tensor.
Require identical full block shape/population in every active macro channel.
Reject partial edges.

Cross-check the projected host tensor against an explicit projection of
Luna's block-one finite operator on small non-cubic and skew cases.  Demonstrate
exact recovery at block one and document approximation/convergence for
nonuniform atom moments as block width decreases.
```

#### [ ] WP10.7b — Terra — integrate accepted coarse blocks

Delegation prompt:

```text
You are Terra.  Integrate Sol's independently tested finite coarse projection
into the accepted OPEN_FFT production path.  Do not change the projection
formula or regenerate Luna's expected values.

Read docs/WP10_OPEN_FFT_BLUEPRINT.md sections 1-7, Sol's WP10.7a report and
tests, Luna's WP10.6b acceptance report, and the accepted periodic WP8
coarse-layout validation.

Remove only the uniform-divisible-coarse rejection.  Before any FFT
allocation, validate that every active macro channel has the expected full
block population and that all active lattice dimensions are exactly divisible
by their block sizes.  Keep partial edge blocks, irregular memberships, and
mixed block populations rejected.  Reuse the accepted active/padded pack,
scatter, normalization, prefactor, and active-only energy paths; do not add a
second atomistic or point-macrospin runtime.

Add GPU and E2E tests for NA=1/2, M=1/4, non-cubic and skew cells, block-one
recovery, nonuniform convergence, exchange coexistence, and every rejection
case.  Record CUDA/HIP status, memory accounting, and the exact input
combinations enabled.  Hand the result to Luna for acceptance.
```

#### [ ] WP10.7c — Luna — accept uniform divisible coarse blocks

Delegation prompt:

```text
You are Luna.  Independently accept WP10 uniform coarse blocks.  Compare the
coarse field and energy to an explicit projection of the independent
block-one finite Hamiltonian, not CPU do_dip=2.  Verify the projected diagonal,
uniform-moment exactness, nonuniform convergence, block-one recovery, NA=1/2,
M=1/4, non-cubic/skew geometry, and deterministic rejection of every partial
edge/population case.  Return GO/NO-GO and numerical tables.
```

Acceptance record:

```text
Sol branch/commit:
Terra integration branch/commit:
Luna report:
Decision:
```

### [ ] WP10.8 — Terra — padding, caching, fp32, CUDA/HIP, and performance

Dependencies: WP10.7c GO, or a documented decision that WP10 ends at an
earlier accepted physical scope.  Optimizations may not weaken that scope's
error budget.

Deliverables:

- padding-size benchmark and policy;
- open-kernel cache;
- setup and steady-state timing/memory tables;
- fp32 accuracy decision;
- CUDA/HIP parity report.

Gate:

- every enabled optimization meets the accepted field and energy budget.

Delegation prompt:

```text
You are Terra.  Optimize accepted OPEN_FFT without changing its Hamiltonian.

Read docs/WP10_OPEN_FFT_BLUEPRINT.md, Luna's latest acceptance report, the WP9
periodic performance work, and the FFT wrapper/memory tracker.

Benchmark exact-minimum P=2G-1, P=2G, and nearby backend-friendly legal sizes.
Any selected P must satisfy P>=2G-1 and reproduce the same finite operator.
Measure thin films, racetrack-like long grids, and fully 3D samples.  Record
setup, pack/clear, forward FFT, contraction, inverse FFT, scatter, energy,
memory, and total steady-state time.

Implement/cache only after measuring:
- immutable host open spectra keyed by the complete open geometry/model;
- singleton-axis plan optimizations if CUDA and HIP both support them safely;
- reduced construction allocations;
- fp32 device storage/execution after a documented error sweep.

Do not introduce periodic images, truncated real-space sums, stale padding,
different coarse physics, or a hidden fallback.  Compare CUDA and HIP against
the same independent oracle.  Run sanitizer and allocation-accounting checks.

Finish with a docs performance report and an explicit enable/reject decision
for fp32 and every padding/cache optimization.
```

Acceptance record:

```text
Assignee:
Branch/commit:
Performance report:
fp32 decision:
CUDA/HIP result:
```

### [ ] WP10.9 — Luna — final WP10 acceptance and closure

Dependencies: all intended WP10 production scopes have their individual GO.

Deliverables:

- final acceptance document;
- capability and rejection matrix;
- synchronized user/design documentation;
- explicit deferred-work list for WP11 and WP12.

Gate:

- WP10 can be marked complete without implying CPU repair, general irregular
  islands, MC, or 2D periodicity.

Delegation prompt:

```text
You are Luna.  Perform the final WP10 closure review.  This is an acceptance
and documentation-consistency task, not a repair task.

Read this blueprint, every WP10 acceptance/performance report, the final input
parser and startup messages, FFT-dipole_implementation_plan.md,
GPU_FFT_DIPOLE_DESIGN.md, and tests/dipole_validation/README.md.

Re-run the accepted numerical matrix on available CUDA and HIP backends,
including fp32 only if it was explicitly accepted.  Verify the capability
matrix exactly distinguishes:
- OPEN_FFT versus EWALD3D_FFT;
- atomic baseline versus basis-resolved block one versus coarse blocks;
- all-open versus mixed/periodic BC;
- SD versus MC;
- uniform divisible grids versus partial/irregular layouts.

Confirm all documentation uses the existing macrocell representation and does
not claim a second atomistic backend.  Confirm "atomic baseline" means only
block=(1,1,1), NA=1.  Confirm WP11 CPU repair and WP12 2D periodicity remain
open and independently gated.

Write the final acceptance report with commits, commands, hardware, numerical
maxima, memory/performance tables, sanitizer results, and known limitations.
Mark WP10 complete only if no required work remains.
```

Acceptance record:

```text
Assignee:
Branch/commit:
Final report:
Decision:
```

## 8. Definition of done for WP10

WP10 is complete only when:

- `gpu_dipole_mode OPEN_FFT` is explicit and `OFF` remains the default;
- legacy `do_dip=0` is mandatory and double application is impossible;
- BC must be exactly `0 0 0`;
- the existing basis-resolved macrocell representation is the sole source
  representation;
- the atomic baseline is explicitly block `(1,1,1)`, `NA=1`;
- active and FFT grids/counts are distinct throughout pack, FFT, scatter, energy,
  diagnostics, cache, and memory accounting;
- padding is cleared on every evaluation;
- the complete finite point tensor has no Ewald or surface terms;
- field and energy come from the same operator;
- the physical prefactor and inverse normalization each occur once;
- independent finite references accept every enabled basis, ensemble,
  geometry, precision, and coarse scope;
- CUDA and HIP meet documented budgets;
- periodic GPU dipoles and `OFF` regressions remain accepted;
- unsupported MC, mixed/periodic BC, partial edges, irregular layouts, and
  irrelevant Ewald inputs fail before allocation; and
- documentation does not claim that WP11 or WP12 is complete.

## 9. WP11 retained analysis — legacy CPU audit and repair

WP11 is a separate project and should receive its own detailed blueprint after
WP10's independent finite oracle is established.

### 9.1 Frozen CPU semantics

```text
do_dip=1  finite/open direct point-dipole sum
do_dip=2  finite projected legacy macrocell Hamiltonian
do_dip=3  finite/open zero-padded FFT convolution
```

None is a periodic Ewald mode.  `do_dip=3` currently ignores BC because the
boundary flags are not passed into its setup or evaluation.  WP11 should
require all-open BC explicitly rather than silently retaining finite physics
under periodic input.

`do_dip=3` is the existing CPU implementation of the same broad numerical
idea as WP10: build a finite point-dipole tensor on the
`NA x N1 x N2 x N3` lattice, zero-pad each axis to `2*Ni-1`, and apply the
linear convolution by FFT.  It is therefore useful as a diagnostic comparison
once repaired.  It is not the production acceptance oracle.

Its source representation must not be confused with either CPU `do_dip=2` or
the GPU contract:

- `dipole_setup()` passes `NA`, `N1/N2/N3`, `Bas`, and the cell vectors to the
  FFTW/MKL builders; it does not pass the legacy macrocell membership;
- `dipole_field_calculation()` passes the full atom-moment array `emomM`
  directly to the FFT implementation; and
- `Num_macro`, `cell_index`, `emomM_macro`, and `Qdip_macro` are used by
  `do_dip=2`, not by the `do_dip=3` FFT call.

Thus CPU `do_dip=3` and GPU `OPEN_FFT` should agree in the atomic baseline
`block=(1,1,1), NA=1`, but WP10 must still use only the existing
basis-resolved GPU macrocell representation.  No CPU atom-array backend is to
be copied into the GPU path.

### 9.2 Confirmed defects

- FFTW allocates its kernel temporary as
  `(N1_pad,N2_pad,N2_pad)` instead of
  `(N1_pad,N2_pad,N3_pad)`.
- FFTW input packing uses an `I1`-contiguous flattened index, while field
  unpacking uses an `I3`-contiguous expression.
- Symmetric cubic/square fixtures can hide this permutation.
- A local `3x2x1`, uniform out-of-plane test on 2026-07-28 showed a maximum
  affected-site field error of `30.5903%` and a dipole-energy error of
  `-8.88328%` versus `do_dip=1`.
- `do_dip=3` does not add dipole energy in
  `dipole_field_calculation()`, although modes 1 and 2 do.
- FFTW and MKL create/destroy transform descriptors repeatedly instead of
  owning reusable batched plans.
- Both retain full complex grids and large automatic temporaries.
- overflow, allocation, plan, and backend status handling is incomplete.
- driver/local-update paths outside ordinary spin dynamics generally describe
  only modes 1 and 2; especially MC must not be assumed correct for mode 3.
- `do_ewald` is parsed but not wired into ordinary Hamiltonian dispatch.
  `ewaldmom.f90` is not an acceptance oracle and requires a separate audit.

### 9.3 Preliminary WP11 checklist

- [ ] **WP11.0 — Luna:** freeze independent finite CPU fixtures and reproduce
  every defect with non-cubic, `G3>G2`, skew, `NA=2`, nonuniform, and `M=4`
  cases.
- [ ] **WP11.1 — Sol:** specify one backend-neutral indexing, normalization,
  field-addition, energy, lifecycle, and BC contract.
- [ ] **WP11.2 — Terra or CPU owner:** repair FFTW allocation/unpack and make
  FFTW/MKL share the common helpers.
- [ ] **WP11.3 — CPU owner:** reuse plans, consider R2C/C2R, remove unsafe
  automatic temporaries, and add checked arithmetic/status handling.
- [ ] **WP11.4 — Luna:** accept FFTW and MKL separately against the independent
  direct finite oracle.
- [ ] **WP11.5 — maintainer:** explicitly gate supported drivers and decide
  whether a new CPU periodic mode should reuse the accepted periodic host
  builder.  Never change `do_dip=3` semantics into periodic Ewald.

## 10. WP12 retained analysis — true 2D-periodic dipoles

WP12 begins only after WP10 is accepted.  It is not open padding in all axes
and not ordinary 3D Ewald with an arbitrary vacuum gap.

### 10.1 Proposed contract

Prefer the public name:

```text
gpu_dipole_mode EWALD2D_FFT
```

over the older staged name `SLAB_PME` if the implementation, like the 3D
solver, precomputes a complete regular-lattice Green function rather than
using a general particle mesh.

For open axis `o` and periodic axes `p,q`:

```text
active grid:   Gp x Gq x Go
FFT grid:      Gp x Gq x Po
Po >= 2*Go-1
```

The periodic directions are cyclic; the open direction is a linear
convolution.  For UppASD's fixed regular lattice, the preferred architecture
is:

```text
initialization:
    converged true 2D-periodic dipole tensor
    K(dp mod Gp, dq mod Gq, do, target basis, source basis)
    -> embed open-axis displacement
    -> FFT and normalize once

runtime:
    existing macro moments
    -> periodic/open mixed FFT convolution
    -> active scatter and energy
```

The special in-plane `k=0` term, slab surface convention, self correction,
and open-axis asymptotics must be derived and tested explicitly.

### 10.2 Candidate mathematical routes

1. **Direct true 2D Ewald/periodic Green function — preferred for the first
   exact regular-lattice builder.**

   Construct the complete 2D-periodic displacement tensor directly and use
   the FFT only to apply the already converged discrete operator.  This most
   closely mirrors accepted `EWALD3D_FFT`.

2. **3D Ewald plus electrostatic layer correction.**

   A dipolar EW3DLC formulation exists and includes correction/error
   estimates.  The gap and layer correction are inseparable parts of the
   method; plain vacuum padding is not valid.

3. **Spectral 2-periodic Ewald.**

   Mixed Fourier-series/Fourier-integral methods offer spectrally accurate
   `O(N log N)` electrostatic solvers.  For dipoles, the differentiated Green
   tensor and special zero mode still require an independent derivation and
   validation.

### 10.3 Required WP12 validation

- an independent true 2D-periodic direct/Ewald oracle;
- invariance to the Ewald split and construction cutoffs;
- lateral lattice-translation invariance;
- convergence as the direct lateral image extent increases;
- explicit zero-mode and slab-surface tests;
- comparison between true 2D and any layer-corrected 3D construction;
- skew in-plane cells and all three choices of open axis;
- `NA=1/2`, `M=1/4`, field, energy, reciprocity, and derivatives;
- proof that increasing vacuum without the correction is not the implemented
  acceptance definition; and
- CUDA/HIP execution using the mixed periodic/open layout accepted by WP10.

### 10.4 Preliminary WP12 checklist

- [ ] **WP12.0 — Sol:** derive and freeze the 2D-periodic dipolar Hamiltonian,
  zero mode, surface convention, and construction-error contract.
- [ ] **WP12.1 — Luna:** implement an independent slow true-2D oracle and
  analytic limiting cases.
- [ ] **WP12.2 — maintainer:** choose direct 2D Green function versus dipolar
  EW3DLC for the first production builder and freeze the selector name.
- [ ] **WP12.3 — Sol/Terra:** implement the host builder and mixed
  periodic/open GPU application only after WP12.0-WP12.2 pass.
- [ ] **WP12.4 — Luna:** perform independent acceptance and compare any
  layer-corrected route with the true 2D oracle.

### 10.5 Literature starting set

- A. J. Newell, W. Williams, and D. J. Dunlop, "A generalization of the
  demagnetizing tensor for nonuniform magnetization," *J. Geophys. Res.*
  **98**, 9551 (1993), DOI
  <https://doi.org/10.1029/93JB00694>.  This supports continuum cuboid
  demagnetizing tensors, not a silent replacement of WP10's point model.
- A. Bródka and A. Grzybowski, "Electrostatic interactions in computer
  simulations of a three-dimensional system periodic in two directions:
  Ewald-type summation," *J. Chem. Phys.* **117**, 8208 (2002), DOI
  <https://doi.org/10.1063/1.1513151>.
- A. Bródka, "Ewald summation method with electrostatic layer correction for
  interactions of point dipoles in slab geometry," *Chem. Phys. Lett.*
  **400**, 62 (2004), DOI
  <https://doi.org/10.1016/j.cplett.2004.10.086>.
- D. Lindbo and A.-K. Tornberg, "Fast and spectrally accurate Ewald
  summation for 2-periodic electrostatic systems," *J. Chem. Phys.*
  **136**, 164111 (2012), DOI
  <https://doi.org/10.1063/1.4704177>.

These papers establish candidate mathematics.  None is permission to copy an
electrostatic charge formula into the magnetic dipole path without deriving
the tensor, self, zero mode, field, and energy conventions.
