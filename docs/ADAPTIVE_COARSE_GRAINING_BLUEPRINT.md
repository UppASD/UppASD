# Adaptive two-scale coarse graining for UppASD

## Physics, numerical, software-engineering, validation, and execution blueprint

**Status:** Working blueprint for discussion and staged implementation  
**Date:** 2026-07-29  
**Primary scope:** Optional low-temperature two-scale spin dynamics on CPU and GPU  
**Explicitly out of scope:** Reimplementation, integration, or redesign of the existing μASD package

---

## 1. Executive decision

UppASD should replace the obsolete constellation optimization with an optional,
two-scale adaptive spin-dynamics package built around the same regular spatial
blocks used by the FFT dipole implementation.

Each spatial block is in one of three operational states:

1. **Atomistic:** all atomic moments in the block are active ASD degrees of
   freedom.
2. **Coarse:** the block is evolved through one or more coarse magnetic
   channels.
3. **Buffer/interface:** atomistic resolution is retained to mediate coupling
   between atomistic and coarse interiors.

The first production target is deliberately narrow:

- low-temperature, fixed-length moment dynamics;
- low-energy localized textures such as domain walls, skyrmions, hopfions, and
  vortices in an otherwise smooth host;
- regular supercell geometries for which the FFT block layout is valid;
- a static coarse/fine mask first, followed by a simple adaptive mask;
- one ferromagnetic dynamical channel first, while the topology and APIs are
  multi-channel from the beginning;
- uniformly coarse FFT dipole fields, including in atomistic blocks;
- CPU reference implementation followed by a GPU production backend.

The initial implementation is not a complete finite-temperature
atomistic–continuum method. It is also not a replacement for μASD. It is a
thin, optional, regular-block alternative intended to obtain substantial
savings when most of the simulated volume is magnetically smooth.

The coarse interior should support two interchangeable operator families:

- a **micromagnetic tensor operator**, using exchange-stiffness,
  spiralization, anisotropy, moment, and inter-channel parameters extracted
  from the atomistic Hamiltonian; and
- a **smooth energy-projection operator**, using a smooth prolongation from
  block variables to atom positions and an adjoint restriction of atomistic
  fields.

The tensor operator is the intended fast production model. The smooth
projection is the principal correctness/reference model and a useful
interface construction. A piecewise-constant rigid macrospin projection is
permitted only as a diagnostic and block-size-one limit, not as the default
coarse exchange model.

---

## 2. Non-negotiable design principles

### 2.1 Optional means behaviorally isolated

When adaptive coarse graining is disabled:

- the normal atomistic CPU and GPU paths must be unchanged;
- atom ordering and existing Hamiltonian storage must be unchanged;
- no coarse-graining allocation, preprocessing, selector work, or device
  synchronization may occur;
- field and energy results must match the pre-feature baseline to the existing
  numerical tolerance;
- μASD behavior must be unchanged.

### 2.2 One canonical spatial partition

The spatial blocks used by adaptive dynamics are exactly the spatial blocks
used by the FFT dipole path. The block sizes
`block_size_x`, `block_size_y`, and `block_size_z` therefore determine:

- FFT dipole resolution;
- coarse-dynamics resolution;
- coarse/fine selection granularity;
- the basic update and compaction unit on the GPU.

The canonical topology must exist independently of whether an FFT dipole
backend is enabled. A CPU run may use the same spatial partition without
performing an FFT dipole calculation.

### 2.3 Spatial blocks are not magnetic channels

The following counts and mappings must remain distinct:

- number of spatial blocks;
- number of crystallographic basis sites, `NA`;
- number of FFT source channels per spatial block;
- number of coarse dynamical channels per spatial block;
- number of chemical species or Hamiltonian types.

For a regular replicated supercell, `NA` is extremely useful: each complete
spatial block contains the same number of copies of every basis site. It must
nevertheless not become a permanent synonym for the number of dynamical
channels.

### 2.4 Energy and field must form an adjoint pair

Every coarse and mixed interaction must have a stated energy functional and a
field obtained from its derivative. If a prolongation operator \(P\) maps
coarse variables to atomistic variables, the corresponding field restriction
must use the adjoint \(P^T\), with the correct magnetic-moment and volume
weights.

This contract is more important than whether a particular implementation is
called atomistic, projected, finite-volume, or micromagnetic.

### 2.5 No silent unsupported physics

Unsupported combinations must be rejected during setup. In particular, the
first implementation should reject rather than approximate silently:

- `NA=Natom` explicit-device geometries;
- arbitrary nonuniform block populations;
- finite-temperature coarse dynamics;
- LSF/induced moments;
- spin-lattice dynamics;
- Monte Carlo and GNEB evolution through the new coarse layer;
- unsupported multi-channel material models;
- partial edge FFT blocks;
- coarse blocks larger than a configured or derived resolution safety limit.

### 2.6 Validation precedes optimization

The CPU reference, analytic tests, and finite-difference energy derivatives
must be accepted before compact GPU kernels replace the reference loops.

---

## 3. Current repository baseline

### 3.1 Existing macrocell layouts

`source/CoarseGraining/macrocells.f90` currently contains two related but
different layouts:

- the legacy `Num_macro/cell_index` layout used by CPU macrocell dipole code,
  which groups all basis atoms in a spatial block; and
- the basis-resolved `pme_Num_macro/pme_cell_index` layout used by the GPU FFT
  dipole code, which keeps one FFT channel per basis site.

The basis-resolved layout is the better starting point for a canonical regular
topology. It already provides the useful identity

\[
  \text{FFT channel} =
  \text{basis id} + NA \times \text{spatial block id}.
\]

The topology work should consolidate the shared spatial information without
forcing the CPU legacy dipole code to change behavior immediately.

### 3.2 Hamiltonian block metadata

`source/CoarseGraining/hamiltonianmacroblocks.f90` creates:

- atom-to-block and block-local mappings;
- CSR block membership;
- directed block-pair adjacency;
- atomistic Hamiltonian entries grouped by source and destination block.

It deliberately preserves atom and neighbour-list order. This is a good
property and should be retained. The module is currently metadata only and
has no physics consumer. It should either evolve into the canonical topology
builder or become an adapter to it.

### 3.3 GPU FFT dipole path

The GPU Hamiltonian code already implements:

1. atom-to-macro moment reduction on the device;
2. FFT evaluation on the macro grid;
3. macro-field scatter to member atoms;
4. optional energy accumulation.

This supplies useful device-side primitives, but all atomic degrees of freedom
are still evolved and all short-range fields are still atomistic. The new
package must reduce the active dynamical and short-range Hamiltonian work, not
only reuse the dipole reduction.

### 3.4 Existing stiffness and spiralization code

`source/SpinWaves/stiffness.f90` already implements:

- Pajda-style regularized exchange lattice sums;
- scalar and \(3\times3\) exchange-stiffness tensors;
- internal `NA×NA` and `3×3×NA×NA` intermediate matrices;
- DMI spiralization through `DMI_stiffness`;
- ordered and random-alloy paths;
- rational and least-squares extrapolations.

This is valuable existing physics. It is not yet a suitable runtime material
service because:

- its public results collapse multi-basis matrices through selected
  eigenvalues;
- raw per-basis/per-channel tensors are internal work arrays;
- the code is documented primarily for ferromagnets;
- signs are inferred from input moment directions;
- the calculation selects a central unit cell using replicated-supercell atom
  ordering;
- calculation, fitting, output, units, and memory management are tightly
  coupled;
- its formulas and prefactors require tests against the Hamiltonian convention
  used by the runtime solver;
- multi-sublattice acoustic and optical modes cannot be represented by simply
  selecting one maximum eigenvalue.

The extraction code should be refactored and validated, not replaced
unnecessarily and not consumed blindly.

The regularization follows the approach introduced by
[Pajda et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.64.174402).
Modern micromagnetic tensor derivations and spin-spiral expansions provide
additional validation routes; see
[Mankovsky, Polesya, and Ebert](https://arxiv.org/abs/1810.13175).

### 3.5 Obsolete constellation implementation

`source/CoarseGraining/optimizationroutines.f90` is an unrelated region
optimization based on classifying similar unit-cell constellations and
reusing representative configurations.

It must be removed before the new package is implemented because:

- it is incomplete and unsupported on the GPU;
- it has propagated obsolete data and arguments through many drivers;
- it is not a spatial atomistic–continuum model;
- at least one active optimized Heisenberg branch loops over neighbours while
  repeatedly adding an already accumulated constellation field, apparently
  multiplying the result by neighbour count;
- retaining the same names would make the new design harder to understand and
  test.

Removal is a separate prerequisite task, not part of the new feature patch.

The normal no-touch rule for μASD has two explicitly approved exceptions:

- remove obsolete constellation arguments and calls from `source/ms_driver.f90`;
- remove obsolete constellation arguments, imports, and branches from
  `source/Multiscale/evolution_ms.f90`.

The already commented constellation arguments in
`source/Multiscale/midpoint_ms.f90` may be deleted as dead comments. No other
μASD physics, interpolation, damping-band, setup, or file-format behavior is
in scope.

### 3.6 Explicit/device geometries

Some device-like inputs use `NA=Natom` instead of a small unit cell replicated
by `N1×N2×N3`. The current basis-resolved PME layout then treats every atom as
a separate basis channel and cannot form useful spatial coarse blocks.

Initial support may reject this geometry, but the new data model must not make
future support impossible. This requires explicit geometry capability flags
and mapping arrays rather than pervasive assumptions that channel count equals
`NA`.

---

## 4. Physical model

### 4.1 State variables

Let \(b\) denote a spatial FFT block and \(\alpha\) a coarse dynamical channel.
For each ensemble, a coarse block holds:

- a unit direction \(\mathbf m_{b\alpha}\);
- a total magnetic moment \(\mu_{b\alpha}\);
- an effective gyromagnetic ratio \(\gamma_{b\alpha}\);
- an effective damping parameter \(\lambda_{b\alpha}\);
- optionally an unresolved angular variance used only for reconstruction and
  diagnostics.

Atomistic blocks retain the usual atomic directions and moment magnitudes.
Dormant atomic arrays are retained in coarse blocks during the first
implementation. This preserves the existing memory model and makes refinement
possible, but it means initial savings are computational rather than
atom-storage savings.

### 4.2 Exploiting `NA` in regular supercells

For a regular supercell and a complete block of
\(B_xB_yB_z\) unit cells:

- the block has \(NA B_xB_yB_z\) atoms;
- each basis site occurs exactly \(B_xB_yB_z\) times;
- restriction can accumulate basis-resolved moments cheaply and exactly;
- a basis-to-dynamical-channel mapping can merge equivalent basis sites;
- a basis-to-FFT-channel mapping can preserve the present dipole convention.

The topology should therefore store:

```text
atom -> spatial_block
atom -> basis_id
atom -> dynamic_channel_id
atom -> fft_channel_id
basis_id -> default_dynamic_channel_id
```

Default regular behavior may use one dynamical channel per magnetic basis
site. Users or material setup may merge equivalent sites into a smaller
number of channels. Nonmagnetic sites must not create dynamical channels.

Automatic inference from the sign of an initial moment is not a sufficiently
robust permanent interface. A future channel-map input should allow users to
identify magnetic sublattices explicitly.

### 4.3 Ferromagnetic coarse energy

A general continuum energy density for a single normalized magnetization
field may be written schematically as

\[
\mathcal E =
  A_{pq}\,\partial_p\mathbf m\cdot\partial_q\mathbf m
  + D_{kp}\,(\mathbf m\times\partial_p\mathbf m)_k
  + \mathcal E_\mathrm{ani}(\mathbf m)
  - \mu_0 M_s\mathbf H_\mathrm{ext}\cdot\mathbf m
  + \mathcal E_\mathrm{dip}.
\]

Indices \(p,q\) refer to Cartesian spatial directions and \(k\) to spin
components. The block grid is indexed along crystallographic repetition
directions, which may not be orthogonal. The discretization must transform
gradients using the physical block vectors and their metric; it must not
assume that block index directions are Cartesian.

The initial finite-volume/finite-difference operator should support:

- a symmetric exchange-stiffness tensor;
- a full spiralization tensor with an explicitly tested sign convention;
- existing uniaxial/multiaxial anisotropy only where a well-defined block
  energy can be constructed;
- external field;
- the existing FFT dipole field.

### 4.4 Multi-channel ferri- and antiferromagnets

For multiple dynamical channels, the coarse energy becomes matrix-valued:

\[
\begin{aligned}
\mathcal E ={}&
 \sum_{\alpha\beta}
 A^{\alpha\beta}_{pq}
 \,\partial_p\mathbf m_\alpha\cdot\partial_q\mathbf m_\beta \\
&+\sum_{\alpha\beta}
 D^{\alpha\beta}_{kp}
 \,(\mathbf m_\alpha\times\partial_p\mathbf m_\beta)_k \\
&+\sum_{\alpha<\beta}
 \Lambda_{\alpha\beta}
 \,\mathbf m_\alpha\cdot\mathbf m_\beta
 + \sum_\alpha \mathcal E_{\mathrm{ani},\alpha}
 + \cdots .
\end{aligned}
\]

Here \(\Lambda_{\alpha\beta}\) retains strong local intersublattice exchange,
while \(A^{\alpha\beta}\) and \(D^{\alpha\beta}\) govern long-wavelength
spatial modes.

The first implementation need not supply all these terms, but all topology,
storage, field APIs, restart design, and GPU indexing must include a channel
dimension. The code must never normalize the net block moment to represent a
compensated material.

For a two-sublattice antiferromagnet, evolving the two sublattice directions
directly is the safest first formulation. Néel/net variables may be added as
an alternative later. Continuum antiferromagnetic theory naturally contains
both soft staggered and hard net-magnetization modes; see the two-sublattice
discussion in
[Phys. Rev. B 104, 064415](https://link.aps.org/accepted/10.1103/PhysRevB.104.064415).

### 4.5 Effective block dynamics

Each coarse channel should obey the same selected deterministic LLG
convention as the atomistic simulation, using a field derived from the coarse
energy:

\[
\mathbf H^\mathrm{eff}_{b\alpha}
= -\frac{\kappa}{\mu_{b\alpha}}
\frac{\partial E}{\partial\mathbf m_{b\alpha}},
\]

where \(\kappa\) is the required UppASD energy-to-field unit conversion and
\(\mu_{b\alpha}\) is the total moment associated with that channel. If the
runtime variables use a different weighted coordinate, the corresponding
dual weight must replace this schematic prefactor. The exact expression must
be derived from UppASD's existing Hamiltonian and field conventions and
verified by finite-difference energy derivatives.

For a channel made of atoms with different Landé factors, the effective
gyromagnetic ratio must follow total angular momentum rather than an
unqualified arithmetic average. Schematically,

\[
\gamma^\mathrm{eff}_{b\alpha}
=
\frac{\sum_i \mu_i}
     {\sum_i \mu_i/\gamma_i}.
\]

The effective damping requires a corresponding angular-momentum-weighted
derivation. Human physics review is mandatory before heterogeneous
\(\gamma_i\) or damping is enabled.

At the first milestone, require identical gyromagnetic ratio and damping
within each dynamical channel and reject violations.

### 4.6 Tensor extraction

The improved stiffness service should produce a typed material descriptor:

```text
coarse_material
  n_channels
  channel_moment
  channel_gamma
  channel_damping
  local_exchange(channel, channel)
  exchange_stiffness(space, space, channel, channel)
  spiralization(spin, space, channel, channel)
  anisotropy descriptors
  units and convention metadata
  extraction diagnostics
```

Extraction must offer two independent checks:

1. real-space lattice sums, including the existing regularized/extrapolated
   method; and
2. direct small-\(q\) comparison against atomistic spin-spiral energies or the
   spin-wave dynamical matrix.

The second check is essential for long-range, frustrated, chiral, or
multi-basis materials. Even accepted relations between atomistic parameters,
spiralization, and helical period may be subtle in real chiral magnets; see
[Grytsiuk et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.214406).

### 4.7 Smooth energy projection

Let \(P\) prolong coarse channel directions to atom positions. A projected
coarse energy is

\[
E_P(\mathbf M)=E_\mathrm{atom}(P\mathbf M),
\qquad
\mathbf H_P=
W_c^{-1}P'(\mathbf M)^T W_a
\mathbf H_\mathrm{atom}(P\mathbf M).
\]

The prolongation should use compact smooth shape functions on the regular
block grid, initially trilinear interpolation in fractional block
coordinates, followed by carefully defined normalization.

Here \(P'(\mathbf M)\) is the derivative of the possibly normalized
prolongation, while \(W_a\) and \(W_c\) contain the atomistic and coarse
moment/unit weights. This weighted adjoint form is schematic; its exact sign
and conversion factors are fixed by the field convention from Section 4.5.

A piecewise-constant prolongation makes a smooth coarse rotation occur as
atomic-scale jumps at block faces. Its exchange energy has the wrong
block-size dependence and it should not define the production coarse model.

The smooth operator has three intended uses:

- an all-coarse reference for validating the tensor discretization;
- a fallback for materials whose continuum tensor extraction is not yet
  supported;
- an energy-consistent interface operator connecting atomistic and tensor
  interiors.

Because normalization makes \(P\) nonlinear, its derivative—not merely the
linear interpolation weights—must be considered in strict energy/field
checks. A tangent-plane linearized version may be accepted initially if its
domain of validity and error are quantified.

### 4.8 Atomistic–coarse interface

The initial interface should be conservative and deliberately wider than the
minimum possible:

1. Dilate every atomistic block by at least the maximum atomistic interaction
   radius, rounded up in blocks.
2. Add one optional safety block beyond that radius.
3. Evaluate normal atomistic interactions throughout the atomistic and buffer
   regions.
4. Generate coarse-side ghost/dormant atomic directions using the smooth
   prolongation.
5. Apply atomistic fields to real fine atoms.
6. Accumulate the reaction field on coarse variables through the adjoint
   prolongation.
7. Suppress the corresponding tensor-coarse energy in the overlap/interface
   band to prevent double counting.

This is an overlapping buffer construction, not a claim of a mathematically
perfect sharp interface. Uniform states and constant long-wavelength spirals
must pass interface patch tests.

Published atomistic–continuum magnetization methods show that ghost torques,
padding/interpolation, and unresolved wave reflection are central interface
issues, not implementation details. Relevant references include:

- [Poluektov, Eriksson, and Kreiss](https://www.sciencedirect.com/science/article/abs/pii/S0045782517302463);
- [De Lucia et al.](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.184415);
- [Andreas, Kákay, and Hertel](https://juser.fz-juelich.de/record/201878/files/PhysRevB.89.134403-1.pdf).

Waves above the coarse-grid cutoff cannot be transmitted faithfully. The
first package should report this limitation rather than reproduce the μASD
damping-band machinery.

### 4.9 Dipole treatment

The first model keeps a single regular FFT dipole grid for the entire system:

- every spatial block contributes through the existing FFT channel mapping;
- all atoms in a spatial block receive the corresponding coarse dipole field;
- a block becoming atomistic does not change dipole resolution;
- no adaptive/nonuniform FFT is attempted;
- no atomistic local dipole correction is attempted initially.

This makes the method a pure two-scale model for short-range dynamics but not
for magnetostatics. The block size must therefore resolve all texture features
whose internal magnetostatic structure matters.

The spatial topology may be shared while FFT and dynamical channel maps remain
different. For example, several basis-resolved FFT channels may be summed or
merged into fewer equivalent ferromagnetic dynamical channels.

### 4.10 Restriction and refinement

#### Restriction

Restriction must be channel resolved:

\[
\mathbf M_{b\alpha}
=\sum_{i\in(b,\alpha)}\mu_i\mathbf s_i.
\]

The package should store both the vector resultant and the sum of moment
magnitudes. Their ratio is a polarization diagnostic. A block with
insufficient polarization must not be coarsened unless its multi-channel
model explains the cancellation.

#### Initial refinement

The accepted initial reconstruction is constrained cone-angle sampling:

1. interpolate a mean direction for each channel from the refined block and
   neighbouring active variables;
2. generate reproducible transverse random perturbations;
3. balance or remove the mean transverse component;
4. renormalize individual moments;
5. iteratively correct the resultant to match the parent coarse moment;
6. taper the cone angle toward an already atomistic neighbour;
7. measure and report the energy change.

The default cone angle at nominal zero temperature should be zero. Random-cone
reconstruction must be opt-in or explicitly configured.

The random seed should be a deterministic function of:

- global seed;
- ensemble;
- spatial block id;
- dynamical channel;
- refinement epoch.

Future reconstruction schemes may perform local ASD or MC thermalization
conditioned on the coarse state. They remain outlook items because correct
finite-temperature magnetization length and conditional microscopic sampling
are nontrivial. See
[Arjmand, Poluektov, and Kreiss](https://link.springer.com/article/10.1007/s10444-017-9575-3).

---

## 5. Adaptive decision system

### 5.1 Separation of responsibilities

The adaptive system has four independent layers:

1. **Criteria:** calculate block scores or hard resolution requests.
2. **Combination policy:** combine multiple criteria.
3. **State machine:** apply thresholds, hysteresis, dwell time, and update
   cadence.
4. **Spatial safety policy:** dilate atomistic regions into buffer regions.

No criterion may directly allocate arrays, change block state, or reconstruct
atoms.

### 5.2 Initial criterion

The first criterion should be maximum local neighbour misalignment:

\[
q_b =
\max_{\substack{i\in b\\j\in\mathcal N(i)}}
\left(1-\mathbf s_i\cdot\mathbf s_j\right).
\]

It is inexpensive, rotationally invariant, and directly related to the
failure of a long-wavelength description.

The initial state policy should use:

- `refine_threshold`;
- a lower `coarsen_threshold`;
- `minimum_dwell_updates`;
- `selector_update_interval`;
- atomistic/buffer dilation measured in blocks.

### 5.3 Planned criteria

The registry should later support:

- channel-resolved polarization;
- maximum or RMS angular variance;
- projected-versus-atomistic field residual;
- atomistic exchange/DMI energy density;
- estimated continuum truncation error;
- topological-density or emergent-field indicators;
- user-defined static masks;
- moving geometric regions of interest;
- material, defect, surface, and interface exclusion masks.

Criteria should return diagnostics even when they are not part of the active
combination policy.

### 5.4 Combination

The initial combination is logical OR for refinement and logical AND for
coarsening:

- any hard or scored criterion may request refinement;
- every active criterion must permit coarsening;
- spatial safety dilation is applied after the raw decision;
- a block with an unsupported material/channel model is permanently
  atomistic.

Weighted-score combinations may be introduced later, but should not be the
first implementation.

### 5.5 State transitions

Allowed transitions are:

```text
COARSE -> PENDING_REFINE -> ATOMISTIC
ATOMISTIC -> PENDING_COARSEN -> COARSE
ATOMISTIC/COARSE -> BUFFER  (derived spatial state)
BUFFER -> ATOMISTIC/COARSE  (after dilation is recomputed)
```

Transitions occur only at accepted synchronization points between complete
integrator steps. They must not occur between predictor and corrector stages.

Every transition log should include:

- step and physical time;
- block and channel;
- reason and criterion scores;
- previous dwell age;
- energy before and after;
- resultant moment before and after;
- reconstruction scheme and seed.

---

## 6. Numerical contracts

### 6.1 Required invariants

The implementation must maintain:

- unit norm of every active spin/channel direction within solver tolerance;
- exact atom membership: each atom belongs to one spatial block;
- exact regular-block basis populations in supported geometries;
- valid channel membership for every magnetic atom;
- no empty active dynamical channel;
- consistent energy/field signs and factors of two;
- no double counting across fine, buffer, and coarse terms;
- deterministic state transitions for a fixed seed and backend;
- block-size-one identity for supported atomistic Hamiltonian terms.

### 6.2 Energy derivative test

For every operator and interface term, perturb a normalized active degree of
freedom in two tangent directions and compare:

\[
\frac{E(\mathbf m+\epsilon\mathbf t)-E(\mathbf m-\epsilon\mathbf t)}
     {2\epsilon}
\]

with the implemented effective field projection. Test a range of
\(\epsilon\) to distinguish truncation error from implementation error.

This is a mandatory unit test, not a one-time notebook calculation.

### 6.3 Geometry and metric

Block vectors are

\[
\mathbf B_1=B_x\,a\,\mathbf C_1,\quad
\mathbf B_2=B_y\,a\,\mathbf C_2,\quad
\mathbf B_3=B_z\,a\,\mathbf C_3.
\]

All gradient and face operators must be derived using these physical vectors.
Tests must include:

- orthogonal cubic cells;
- anisotropic orthorhombic cells;
- at least one skew/nonorthogonal cell;
- periodic and open boundaries supported by the chosen FFT mode.

### 6.4 Resolution safety

The block grid must resolve the continuum solution. Where parameters are
available, estimate:

- exchange length;
- domain-wall width;
- expected spiral pitch;
- user-provided minimum feature size.

Setup should warn or reject when too few blocks cover the smallest requested
length scale. The exact threshold is a policy decision, but the diagnostic
must print the physical block dimensions and relevant ratios.

### 6.5 Time stepping

The first CPU reference should use a coarse integrator algebraically
consistent with an existing deterministic atomistic LLG integrator. Atomistic
and coarse variables should initially use the same global time step.

Subcycling can be investigated later, but would complicate:

- interface synchronization;
- selector timing;
- energy accounting;
- GPU stream ordering;
- reproducibility.

### 6.6 Precision

- Material extraction, topology geometry, energies, and validation references
  should use double precision.
- GPU field storage may follow the configured device precision.
- Block reductions in fp32 builds should document and test their accumulation
  error.
- CPU/GPU acceptance tolerances must be precision specific.

---

## 7. Proposed software architecture

### 7.1 Directory intent

After constellation removal, `source/CoarseGraining` should contain only the
new spatial-block package and legacy adapters required by active dipole code.
A suggested decomposition is:

```text
source/CoarseGraining/
  coarsegraining_types.f90
  blocktopology.f90
  blockchannels.f90
  blockrestriction.f90
  blockreconstruction.f90
  blockselectors.f90
  blockstatemachine.f90
  coarsematerial.f90
  coarseoperator_tensor.f90
  coarseoperator_projection.f90
  coarseinterface.f90
  coarseintegrator.f90
  macrocells.f90                 # temporary legacy/FFT adapters
  hamiltonianmacroblocks.f90     # retained or folded into topology
```

GPU implementation files should live with the GPU backend rather than placing
CUDA/HIP logic in the Fortran directory.

### 7.2 Canonical topology type

A conceptual topology type is:

```fortran
type :: block_topology_type
   logical :: ready
   integer :: geometry_mode
   integer :: n_spatial_blocks
   integer :: n_basis
   integer :: n_dynamic_channels
   integer :: n_fft_channels
   integer :: block_shape(3)
   integer :: block_grid(3)

   integer, allocatable :: atom_to_block(:)
   integer, allocatable :: atom_to_basis(:)
   integer, allocatable :: atom_to_dynamic_channel(:)
   integer, allocatable :: atom_to_fft_channel(:)

   integer, allocatable :: block_atom_offset(:)
   integer, allocatable :: block_atoms(:)
   integer, allocatable :: block_channel_population(:,:)

   real(dblprec), allocatable :: block_center(:,:)
   real(dblprec), allocatable :: block_vectors(:,:,:)
   real(dblprec), allocatable :: block_volume(:)

   ! CSR block adjacency and atomistic interaction metadata.
   ...
end type
```

The actual layout may use flattened arrays for C interoperability, but the
conceptual distinctions must remain.

### 7.3 Geometry modes

Define capability modes now:

- `REGULAR_REPLICATED_CELL`: supported initially;
- `EXPLICIT_DEVICE`: represented but rejected initially;
- possible future `IRREGULAR_SPATIAL_BINS`.

For regular mode, require and verify:

\[
N_\mathrm{atom}=NA\,N_1N_2N_3
\]

and the expected atom ordering or build the map from explicit unit-cell
indices if they are available.

For explicit-device mode, future support will require:

- coordinate-based spatial binning in lattice/fractional coordinates;
- explicit or inferred magnetic-channel labels independent of `NA`;
- FFT treatment for arbitrary populations or a separate deposition scheme;
- no allocation proportional to `NA` channels when `NA=Natom`.

### 7.4 State and runtime data

Runtime state should be separate from immutable topology:

```text
block_runtime
  resolution_state(block)
  pending_state(block)
  state_age(block)
  selector_scores(criterion, block)
  coarse_moment(3, channel, block, ensemble)
  coarse_direction(3, channel, block, ensemble)
  coarse_field(3, channel, block, ensemble)
  transition_epoch(block)
  active_atom_list
  active_coarse_list
  interface_atom_list
```

The GPU should consume compact lists. Initial CPU code may use masks for
clarity, but the public runtime contract should not require branch-heavy
per-atom GPU kernels.

### 7.5 Strategy interfaces

Provide registries or narrow abstract interfaces for:

- `coarse_operator`;
- `selector_criterion`;
- `restriction_scheme`;
- `reconstruction_scheme`;
- `selection_combiner`.

For GPU interoperability, strategy selection should resolve to compact enum
and descriptor data during setup. Do not require device-side virtual dispatch.

### 7.6 Fortran/GPU ownership

Recommended ownership:

- Fortran constructs and owns the canonical topology and material descriptor.
- Plain arrays and scalar descriptors are exported through the existing
  Fortran/C bridge.
- GPU code validates counts and capabilities before allocation.
- GPU runtime owns device masks, compact lists, block moments, block fields,
  and operator buffers.
- Dynamic selector updates should remain on the device where practical.
- If the block mask is copied to the host initially, copy only block-scale
  data and only at selector update intervals.

### 7.7 Inputs

Final names should follow UppASD conventions, but a candidate input surface is:

```text
do_adaptive_cg                 Y/N
cg_operator                    TENSOR/PROJECTED
cg_mask_mode                   STATIC/ADAPTIVE
cg_selector                    MAX_ANGLE[, ...]
cg_refine_threshold            <real>
cg_coarsen_threshold           <real>
cg_update_interval             <integer>
cg_minimum_dwell_updates       <integer>
cg_buffer_blocks               <integer>
cg_channel_mode                BASIS/FILE
cg_channel_file                <filename>
cg_reconstruction              ALIGNED/CONE
cg_cone_angle                  <real>
cg_static_mask_file            <filename>
cg_energy_jump_limit           <real>
cg_diagnostics                 <level>
```

The existing block-size inputs are reused. Enabling adaptive coarse graining
must trigger the same uniform/divisible block validation as the supported FFT
layout, even if the dipole FFT is disabled.

### 7.8 Restart and output

Dynamic coarse graining adds state not present in existing restart files.
Before adaptive runs are declared production ready, restart data must include:

- block state and dwell ages;
- transition epochs;
- active coarse directions and moments;
- dormant atom directions if they are not already present;
- selector configuration hash;
- reconstruction RNG state or deterministic epoch information;
- topology/material version.

Initial static-mask development may explicitly reject restart.

### 7.9 Diagnostics and observability

Provide machine-readable or easily parsed diagnostics for:

- topology counts and memory;
- channel populations;
- extracted material tensors and units;
- characteristic-length checks;
- state counts over time;
- transition reasons and energy jumps;
- time spent in restriction, selectors, compaction, atomistic fields, coarse
  fields, interface fields, FFT dipole, and integration;
- achieved active-DOF reduction;
- peak GPU memory.

---

## 8. Verification and validation matrix

### 8.1 Cleanup regression

- Feature-disabled CPU build and standard regression suite pass.
- Feature-disabled CUDA and HIP builds pass where available.
- No constellation symbol remains, except explicitly approved historical
  documentation if retained.
- Removed input keywords fail with an informative message.
- μASD reference inputs produce unchanged results within existing tolerances.

### 8.2 Topology

- `1×1×1`, `NA=1`.
- `NA=2` and `NA>2`.
- non-cubic block shapes.
- multiple ensembles.
- periodic and open boundaries.
- atom-to-block bijection.
- exact per-basis populations.
- block-size-one identity.
- invalid/nondivisible layout rejection.
- `NA=Natom` explicit-device rejection with a precise capability message.

### 8.3 Material extraction

- simple-cubic nearest-neighbour ferromagnet with analytic stiffness;
- anisotropic exchange yielding a diagonal anisotropic tensor;
- skew-cell coordinate transformation;
- controlled long-range exchange demonstrating regularization convergence;
- simple DMI chain/film with known handedness and spiralization sign;
- direct small-\(q\) spin-spiral energy fit;
- `NA=2` ferromagnetic acoustic and optical branches;
- `NA=2` antiferromagnetic branch structure before multi-channel production
  is enabled.

### 8.4 All-coarse operator

- uniform magnetization: zero exchange/DMI torque where physically expected;
- global rotation invariance for isotropic exchange;
- energy derivative checks;
- spin-spiral energy and dispersion versus atomistic reference;
- domain-wall width and energy convergence with block size;
- skyrmion equilibrium radius and energy versus full ASD;
- smooth-projection versus tensor-operator comparison.

### 8.5 Static hybrid

- all-fine mask equals baseline ASD;
- all-coarse mask equals accepted all-coarse reference;
- uniform-state interface patch test;
- constant long-wavelength spiral interface patch test;
- domain wall translated across the interface;
- skyrmion moved across the interface;
- no double-counted energy;
- total reaction torque across the interface is consistent with the energy
  derivative;
- interface error decreases as the coarse grid is refined.

### 8.6 Adaptive transitions

- hysteresis prevents chatter;
- dwell time is respected;
- buffer dilation covers the interaction radius;
- deterministic replay with the same seed;
- aligned reconstruction conserves the block direction;
- cone reconstruction matches the requested resultant;
- transition energy jumps remain below configured limits;
- topology and state remain valid after repeated refine/coarsen cycles;
- selector update does not occur inside an integrator stage.

### 8.7 Multi-channel

- compensated AFM block does not produce normalization failure;
- each sublattice norm remains valid;
- intra-block intersublattice torque has the correct sign;
- acoustic and optical small-\(q\) modes match atomistic references;
- ferrimagnetic net and angular-momentum compensation behavior is sensible;
- FFT source moment agrees with the sum over all dynamical channels.

### 8.8 CPU/GPU parity

- topology and channel maps match exactly;
- static masks match;
- block reductions match within precision tolerance;
- tensor and interface fields match;
- energies match;
- state-transition decisions match, or documented deterministic thresholds
  avoid backend-dependent ties;
- CUDA and HIP use the same indexing and sign conventions;
- fp32 and fp64 acceptance budgets are separate.

### 8.9 Performance acceptance

Measure:

- wall time per step;
- fraction of atomistic blocks;
- active atoms and active block channels;
- selector/compaction overhead;
- memory overhead;
- crossover point where coarse graining becomes beneficial;
- scaling with block size and ensembles;
- CUDA/HIP portability.

No speedup claim should be made from skipping atomic integration alone. The
short-range Hamiltonian work must also be reduced.

---

## 9. Phased delivery and quality gates

### Phase 0: Remove constellation optimization

**Exit gate:** The obsolete implementation and inputs are gone, ordinary
workflows pass, and only the approved μASD template cleanup has occurred.

### Phase 1: Foundations

Deliver canonical topology, channel mapping, immutable/runtime type separation,
input capability validation, and an audited material-extraction API.

**Exit gate:** Topology and analytic stiffness/DMI tests pass without enabling
new runtime dynamics.

### Phase 2: Static CPU reference

Deliver all-coarse tensor and smooth-projection operators, then static
fine/coarse masks with an energy-consistent buffer interface.

**Exit gate:** Energy derivatives, spin spirals, domain walls, and static
interface tests pass.

### Phase 3: Adaptive CPU prototype

Deliver modular criteria, hysteretic state machine, restriction, aligned/cone
reconstruction, transitions, and diagnostics.

**Exit gate:** A domain wall or skyrmion can move through an adaptively refined
region without unacceptable energy jumps or loss of the texture.

### Phase 4: GPU production path

Deliver device state, compact work lists, block operator kernels, interface
kernels, selector/compaction workflow, and CPU/GPU parity.

**Exit gate:** Accepted physics tests pass on CUDA and HIP where available,
and a meaningful crossover/speedup is demonstrated.

### Phase 5: Multi-channel production

Deliver validated ferri-/antiferromagnetic material extraction and coupled
block dynamics.

**Exit gate:** Acoustic/optical modes and at least one ferrimagnetic or
antiferromagnetic texture benchmark match atomistic references.

### Phase 6: Future geometry and thermal work

Investigate explicit-device geometry, arbitrary deposition, local
thermalization, and temperature-dependent coarse models.

This phase is not required for the initial feature.

---

## 10. Delegation guide

The labels below are capability recommendations, not provider-specific
requirements:

- **Human:** physics ownership, convention approval, scope decisions, and
  acceptance of model error.
- **Opus/Terra:** high-reasoning architecture, mathematical derivation,
  codebase-wide review, and adversarial design review.
- **Sol:** substantial implementation, cross-language integration, debugging,
  GPU work, and performance engineering.
- **Luna/Sonnet:** bounded refactors, focused unit tests, mechanical cleanup,
  documentation, and independent verification.

Every task that changes physics conventions requires human approval, even when
implemented by a strong model. Separate implementation and review agents are
recommended for stiffness signs/units, interface energy, and GPU parity.

Agents should work in isolated branches or worktrees for independent tasks.
Two agents should not concurrently edit central driver signatures,
`macrocells.f90`, or GPU Hamiltonian dispatch.

---

## 11. Task prompts and acceptance checklists

The following prompts are intended to be reusable task briefs. Each task
should begin by reading this blueprint and inspecting the current source
rather than assuming filenames and signatures have remained unchanged.

### Task CG-00: Remove the obsolete constellation approach

**Dependencies:** None  
**Suggested primary:** Sol  
**Suggested review:** Luna/Sonnet for symbol audit; Human for μASD boundary  
**Risk:** High mechanical breadth, low intended physics change

#### Prompt

> Remove the obsolete constellation/region-optimization implementation from
> UppASD end to end. Delete its runtime behavior, inputs, storage, branches,
> and argument plumbing while preserving all normal paths. You may edit
> `source/ms_driver.f90` and `source/Multiscale/evolution_ms.f90` solely to
> remove leftover constellation imports, arguments, declarations, calls, and
> branches. You may delete already-commented constellation argument remnants
> in `source/Multiscale/midpoint_ms.f90`. Do not otherwise change μASD
> algorithms, damping bands, interpolation, setup, inputs, or results.
>
> Start from a complete `rg` inventory of `optimizationroutines`,
> `OPT_flag`, `OPT_ON`, `constellation*`, `constlNCoup`, `unitCellType`,
> `cos_thr`, `cellPosNumNeigh`, `buildOptRegions`, `pre_optimize`, and related
> input keywords. Remove `optimizationroutines.f90` from the build and delete
> it after consumers are cleaned. Remove or explicitly reject the old input
> keywords, including `do_region_optimization`,
> `region_optimization_rebuild_step`, and `glob_cos_thr`. Do not silently
> accept an input that requests removed behavior.
>
> Keep the patch focused. Do not begin implementing the new coarse-graining
> package. Build representative CPU and GPU configurations and run existing
> regression and μASD tests. Report any pre-existing failures separately.

#### Checklist

- [x] Complete before/after symbol inventories are attached.
- [x] The faulty optimized Hamiltonian branches are gone.
- [x] `optimizationroutines.f90` is deleted and removed from CMake.
- [x] Old inputs fail clearly or are removed according to input policy.
- [x] Procedure signatures no longer carry dead constellation arguments.
- [x] Allocations and cleanup calls are gone.
- [x] `ms_driver.f90` and `evolution_ms.f90` contain only the approved
      mechanical removals.
- [x] No other μASD source behavior changed.
- [x] CPU build succeeds.
- [x] CUDA and HIP builds succeed where available.
- [x] Ordinary ASD regression tests pass.
- [x] μASD reference tests pass within their existing tolerance.
- [x] No unrelated formatting or refactoring obscures the removal.

---

### Task CG-01: Freeze physics conventions and capability matrix

**Dependencies:** CG-00 may run separately, but must complete before merge  
**Suggested primary:** Opus/Terra  
**Required participant:** Human  
**Suggested review:** Sol for implementability  
**Risk:** High; decisions affect every later task

#### Prompt

> Produce a concise implementation specification that freezes the initial
> adaptive coarse-graining physics conventions. Derive and document:
>
> - the exact coarse energy terms and field derivatives;
> - UppASD unit conversion, sign, and factor-of-two conventions;
> - block moment, effective field, gyromagnetic ratio, and damping
>   definitions;
> - tensor index ordering for exchange stiffness and spiralization;
> - the coordinate transformation for nonorthogonal block vectors;
> - the energy partition between atomistic, interface, and coarse regions;
> - the initially supported Hamiltonian terms;
> - the initial geometry, solver, temperature, and boundary capabilities;
> - the multi-channel storage contract even where physics support is deferred.
>
> Use small analytic systems and the current atomistic Hamiltonian convention
> to check every sign and prefactor. Identify decisions requiring human
> approval rather than selecting them silently. Do not implement production
> kernels in this task.

#### Checklist

- [x] Energy and field equations use one unambiguous index convention.
- [x] Cartesian versus lattice/fractional coordinates are explicit.
- [x] DMI handedness has an analytic fixture.
- [x] Moment and field units match existing UppASD conventions.
- [x] Heterogeneous Landé/damping support is either derived or rejected.
- [x] Fine/interface/coarse energy ownership is explicit.
- [x] Initial capability matrix lists accept/reject behavior.
- [ ] Human approves the convention document.
- [ ] Sol review confirms descriptors can represent every approved term.

**CG-01 evidence:** `docs/CG-01_PHYSICS_CONVENTIONS.md` records the proposed
contract, source-derived UppASD unit/sign conventions, analytic fixtures, and
the two outstanding approval gates.  No production kernels are included in
this task.

---

### Task CG-02: Build the canonical regular block topology

**Dependencies:** CG-00, CG-01  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra for abstraction; Luna/Sonnet for tests  
**Risk:** Medium-high

#### Prompt

> Implement an immutable canonical spatial-block topology for
> `REGULAR_REPLICATED_CELL` geometry. Reuse the spatial partition of the FFT
> dipole blocks, preserve atom ordering, and represent basis, FFT-channel, and
> dynamical-channel mappings separately. Consolidate or adapt the current
> `macrocells.f90` and `hamiltonianmacroblocks.f90` metadata without changing
> active legacy dipole results.
>
> The builder must validate `Natom=NA*N1*N2*N3`, positive divisible block
> dimensions, exact atom membership, uniform per-basis populations, and
> nonempty magnetic channels. Add a geometry capability enum including
> `EXPLICIT_DEVICE`, but reject `NA=Natom` device-like setups for adaptive
> coarse graining with a precise diagnostic. Avoid allocations whose channel
> dimension becomes `NA=Natom` in the rejected path.
>
> Expose plain arrays suitable for Fortran/C/GPU transfer. Add focused unit
> tests independent of the dipole numerical kernels.

#### Checklist

- [ ] Spatial block ids match the FFT block grid.
- [ ] `atom_to_basis`, `atom_to_dynamic_channel`, and
      `atom_to_fft_channel` are distinct.
- [ ] CSR membership covers every atom exactly once.
- [ ] `NA=1`, `NA=2`, and `NA>2` tests pass.
- [ ] Non-cubic block tests pass.
- [ ] Block-size-one mapping is exact.
- [ ] Skew-cell geometry is stored correctly.
- [ ] Invalid/nondivisible cases fail before allocation.
- [ ] Explicit-device geometry fails with a capability message.
- [ ] Existing CPU/GPU dipole tests remain unchanged.
- [ ] Feature-off path allocates no new topology unless another active feature
      requests it.

---

### Task CG-03: Refactor and validate stiffness/spiralization extraction

**Dependencies:** CG-01, topology descriptors from CG-02  
**Suggested primary:** Sol with Opus/Terra physics pairing  
**Required review:** Human  
**Suggested test support:** Luna/Sonnet  
**Risk:** Very high physics sensitivity

#### Prompt

> Refactor `source/SpinWaves/stiffness.f90` so its existing exchange-stiffness
> and DMI spiralization calculations can supply a typed, side-effect-light
> coarse-material descriptor. Preserve the existing user-facing stiffness
> output path unless a verified bug requires a separately documented fix.
>
> Separate calculation from printing, fitting, file units, and module-global
> storage. Expose raw per-basis/per-channel local exchange,
> `A(space,space,channel,channel)`, and
> `D(spin,space,channel,channel)` data before eigenvalue collapse. Record
> units, signs, Hamiltonian convention, convergence parameters, and fit
> diagnostics in the result.
>
> Do not assume the present maximum-eigenvalue reduction is a correct
> multi-sublattice runtime model. Validate real-space sums against analytic
> nearest-neighbour systems and direct small-q atomistic spin-spiral energies.
> Include anisotropic exchange, DMI handedness, a skew cell, and at least one
> two-basis fixture. Gate ferri/AFM runtime use until acoustic and optical mode
> extraction is validated.

#### Checklist

- [ ] Existing stiffness output remains available.
- [ ] Pure calculation API has no mandatory file output.
- [ ] Raw channel matrices/tensors are exposed.
- [ ] Units and tensor ordering are self-describing.
- [ ] Long-range regularization convergence is reported.
- [ ] Analytic FM stiffness test passes.
- [ ] Anisotropic tensor test passes.
- [ ] DMI sign/handedness test passes.
- [ ] Small-\(q\) atomistic energy fit agrees.
- [ ] Nonorthogonal geometry test passes.
- [ ] Multi-basis matrices are not silently reduced to one eigenvalue.
- [ ] Unsupported AFM/ferri consumption fails explicitly.
- [ ] Human signs off signs, prefactors, and units.

---

### Task CG-04: Implement the all-coarse CPU tensor operator

**Dependencies:** CG-01, CG-02, CG-03  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra for numerical stencil; Human for physics  
**Risk:** High

#### Prompt

> Implement a CPU reference tensor coarse operator for a single ferromagnetic
> channel on the canonical regular block grid. Use the accepted material
> descriptor and physical block metric to evaluate exchange stiffness,
> spiralization, supported anisotropy, external field, and the existing
> uniformly coarse dipole field. Derive fields from an explicit discrete
> energy and provide per-term energies for testing.
>
> Start with all blocks coarse and a static configuration; do not add adaptive
> switching or GPU code. Use clear reference loops and double precision.
> Reuse an existing deterministic LLG integrator algebra where practical, but
> keep atomistic and coarse data separate. Reject unsupported material,
> geometry, solver, and temperature combinations during setup.

#### Checklist

- [ ] Uniform-state torque test passes.
- [ ] Global exchange-rotation invariance passes.
- [ ] Energy derivative tests pass for every term.
- [ ] Cubic and skew-cell stencils pass.
- [ ] Small-\(q\) dispersion matches atomistic reference.
- [ ] DMI chirality matches the accepted convention.
- [ ] Domain-wall energy and width converge with block refinement.
- [ ] Block-size-one limit is understood and tested where applicable.
- [ ] Per-term energy reporting is available.
- [ ] Unsupported modes reject cleanly.
- [ ] No adaptive or interface behavior is hidden in the operator.

---

### Task CG-05: Implement the smooth projected CPU operator

**Dependencies:** CG-01, CG-02  
**Suggested primary:** Sol  
**Suggested design review:** Opus/Terra  
**Risk:** High numerical subtlety

#### Prompt

> Implement a double-precision CPU reference operator based on smooth
> prolongation from block/channel variables to atom positions and adjoint
> restriction of atomistic fields. Begin with trilinear shape functions in
> regular fractional block coordinates. Define how spin normalization enters
> the prolongation and its derivative. Evaluate supported atomistic bilinear
> Hamiltonian terms without modifying their source arrays.
>
> Compare the projected energy and field against finite-difference derivatives
> and against the tensor operator for long-wavelength spin spirals. Include a
> piecewise-constant projection only as a test demonstrating why it is not the
> production exchange model. Do not optimize or port to GPU in this task.

#### Checklist

- [ ] Shape weights form the required partition of unity.
- [ ] Boundary handling is explicit.
- [ ] Normalization derivative is implemented or approximation bounds are
      documented.
- [ ] Adjoint energy/field tests pass.
- [ ] Block-size-one behavior is tested.
- [ ] Smooth spin-spiral energy has correct mesh scaling.
- [ ] Rigid projection failure is captured by a regression/diagnostic test.
- [ ] Tensor and projected long-wavelength results agree within an explained
      discretization error.
- [ ] General atom ordering remains unchanged.

---

### Task CG-06: Implement static hybrid CPU coupling

**Dependencies:** CG-04, CG-05  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra and Human  
**Risk:** Very high; interface correctness

#### Prompt

> Implement static fine/coarse masks and a conservative buffer interface on
> the CPU. Atomistic and buffer blocks use normal atomistic Hamiltonian
> evaluation. Coarse interiors use the tensor operator. Coarse-side ghost
> atomic states are supplied by smooth prolongation, and their reaction fields
> are accumulated through the adjoint. Define a non-overlapping ownership of
> energy terms so the interface is not double counted.
>
> Compute buffer width from the maximum active atomistic interaction radius,
> with a configurable safety dilation. Begin with exchange and DMI plus the
> already accepted on-site/coarse terms. Do not add dynamic switching.
> Implement all-fine and all-coarse masks through the same dispatch so both
> limiting cases are testable.

#### Checklist

- [ ] All-fine equals baseline ASD.
- [ ] All-coarse equals accepted coarse reference.
- [ ] Uniform patch test passes.
- [ ] Constant long-wavelength spiral patch test passes.
- [ ] Interface energy derivative test passes.
- [ ] No field or energy term is double counted.
- [ ] Buffer covers every cross-interface atomistic interaction.
- [ ] Domain wall crosses a fixed interface without spurious pinning beyond
      the accepted error.
- [ ] Skyrmion translation test remains stable.
- [ ] Interface error decreases under block refinement.
- [ ] μASD code is not used or modified.

---

### Task CG-07: Implement selector registry and state machine

**Dependencies:** CG-02; may be developed against mock fields before CG-06  
**Suggested primary:** Luna/Sonnet or Sol  
**Suggested review:** Opus/Terra for API  
**Risk:** Medium

#### Prompt

> Implement a modular block selector registry, combination policy, hysteretic
> state machine, and buffer dilation independent of field operators and
> reconstruction. The first score is maximum neighbour misalignment. Support
> static masks, refine/coarsen thresholds, update intervals, dwell ages, hard
> atomistic exclusions, and logical-OR refine/logical-AND coarsen combination.
>
> Criteria may read state and return scores/requests but must not change block
> state. State transitions occur only through the state machine at explicit
> synchronization points. Provide deterministic transition logs and unit
> tests using synthetic score sequences.

#### Checklist

- [ ] Criterion and state mutation APIs are separate.
- [ ] Hysteresis and dwell tests pass.
- [ ] Buffer dilation tests cover boundaries and periodic wrapping.
- [ ] Unsupported blocks remain atomistic.
- [ ] Multiple criteria combine as specified.
- [ ] Static mask works without selector evaluation.
- [ ] Synthetic chatter sequences do not chatter.
- [ ] Transition decisions are reproducible.
- [ ] Logs contain reason and score data.
- [ ] No integrator or Hamiltonian code is embedded in the selector module.

---

### Task CG-08: Implement restriction, reconstruction, and adaptive CPU runs

**Dependencies:** CG-06, CG-07  
**Suggested primary:** Sol  
**Suggested review:** Human for reconstruction error; Luna/Sonnet for tests  
**Risk:** High

#### Prompt

> Connect the selector state machine to the static hybrid CPU solver. Implement
> channel-resolved moment-weighted restriction, aligned reconstruction, and
> deterministic constrained cone-angle reconstruction. Apply transitions only
> between complete integration steps. Rebuild active lists and interface
> ownership safely after each accepted transition.
>
> Cone sampling must use deterministic per-block/channel/ensemble/epoch seeds,
> remove net transverse bias, match the requested block resultant within
> tolerance, and report the transition energy jump. Add a configurable energy
> jump limit that rejects or rolls back a transition without corrupting the
> previous state. Do not add finite-temperature claims or local
> thermalization.

#### Checklist

- [ ] Restriction conserves channel moment vector.
- [ ] Compensated multi-channel data does not normalize a net zero vector.
- [ ] Aligned reconstruction is exact.
- [ ] Cone reconstruction matches the requested resultant.
- [ ] Seeds and replay are deterministic.
- [ ] Failed transitions leave the previous state intact.
- [ ] Energy jumps are measured and logged.
- [ ] Repeated transitions preserve topology invariants.
- [ ] Domain-wall adaptive benchmark passes.
- [ ] Skyrmion adaptive benchmark passes.
- [ ] No transition occurs inside predictor/corrector stages.

---

### Task CG-09: Design GPU runtime data and compact work scheduling

**Dependencies:** Accepted CPU adaptive prototype  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra for architecture  
**Risk:** High

#### Prompt

> Design and implement GPU-side immutable topology transfer and mutable
> adaptive runtime storage without yet replacing every CPU reference
> operator. Add device masks, block/channel moments and fields, compact active
> atom/block/interface lists, and memory accounting. Integrate allocation,
> validation, staging, cleanup, and stream ownership with the existing GPU
> simulation lifecycle.
>
> Avoid per-atom virtual dispatch and persistent host synchronization.
> Selector-update compaction may initially use a simple accepted device scan
> or a block-mask host round trip, but only block-scale data may be transferred
> and the synchronization cost must be measured. Feature-off allocation and
> execution must remain unchanged.

#### Checklist

- [ ] Fortran and GPU topology counts match exactly.
- [ ] Allocation preflight includes every new device buffer.
- [ ] Cleanup restores device-memory accounting.
- [ ] Feature-off device inventory is unchanged.
- [ ] Compact lists contain exactly the expected work.
- [ ] Mask changes rebuild lists safely.
- [ ] Multiple ensembles index correctly.
- [ ] CUDA and HIP compilation paths use the same descriptors.
- [ ] Stream ordering is documented.
- [ ] Host synchronization is measured and localized.

---

### Task CG-10: Implement GPU coarse, interface, and adaptive kernels

**Dependencies:** CG-09 and accepted CPU references  
**Suggested primary:** Sol  
**Suggested independent review:** Sol or Opus/Terra not responsible for the
initial kernel; Luna/Sonnet for parity harness  
**Risk:** Very high

#### Prompt

> Port the accepted tensor coarse operator, interface coupling, restriction,
> selector score, reconstruction, and integration workflow to CUDA/HIP using
> compact work lists. Preserve the CPU energy, sign, tensor-index, and state
> transition contracts. Reuse existing GPU macro reductions and FFT dipole
> dispatch where appropriate, but do not couple adaptive short-range state to
> a new dipole resolution.
>
> Add CPU/GPU parity fixtures before performance tuning. Keep construction and
> runtime allocations visible to the existing memory preflight. Profile
> atomistic, coarse, interface, selector, compaction, FFT, and integration
> phases separately. Optimize only after fp64 correctness is accepted, then
> define fp32 budgets.

#### Checklist

- [ ] FP64 CPU/GPU fields and energies agree.
- [ ] CUDA and HIP parity tests pass where hardware is available.
- [ ] State decisions match for non-threshold-tie fixtures.
- [ ] Restriction and reconstruction are deterministic as specified.
- [ ] Energy derivative fixtures agree with CPU reference.
- [ ] Existing FFT dipole tests still pass.
- [ ] Feature-off performance is unchanged within noise.
- [ ] Selector/compaction overhead is reported.
- [ ] A real active-DOF crossover is measured.
- [ ] FP32 error budgets are documented separately.
- [ ] No untracked device allocations bypass memory accounting.

---

### Task CG-11: Enable validated multi-channel ferri/AFM dynamics

**Dependencies:** CG-03, accepted single-channel CPU/GPU paths  
**Suggested primary:** Human + Opus/Terra physics design; Sol implementation  
**Suggested review:** Independent human/model physics review  
**Risk:** Research-grade

#### Prompt

> Extend the coarse material extraction and tensor operator from one
> ferromagnetic channel to coupled magnetic sublattices. Preserve one
> direction, moment, gyromagnetic ratio, and damping parameter per dynamical
> channel and block. Derive local intersublattice exchange and matrix-valued
> spatial stiffness/spiralization from the atomistic Hamiltonian or accepted
> small-q dynamical matrix. Do not represent a compensated block by its net
> moment.
>
> Begin with a controlled two-sublattice model. Validate uniform ground state,
> acoustic and optical spin-wave branches, sublattice-resolved torques,
> ferrimagnetic compensation behavior, and a domain-wall or skyrmion
> benchmark. Enable runtime use only for channel models with accepted
> extraction diagnostics. Keep general random-alloy and finite-temperature
> claims out of scope unless separately derived.

#### Checklist

- [ ] Channel map can merge multiple basis sites explicitly.
- [ ] Nonmagnetic basis sites are excluded.
- [ ] Intersublattice exchange sign and field derivative pass.
- [ ] AFM zero-net-moment blocks remain numerically valid.
- [ ] Acoustic and optical modes match atomistic results.
- [ ] Per-channel \(\gamma\) and damping are defined.
- [ ] FFT source moment equals channel sum under its mapping.
- [ ] CPU reference passes before GPU enablement.
- [ ] At least one ferri/AFM texture benchmark passes.
- [ ] Human physics acceptance is recorded.

---

### Task CG-12: Design future explicit-device geometry support

**Dependencies:** Stable regular topology and FFT backend  
**Suggested primary:** Opus/Terra  
**Suggested implementation:** Sol in a later phase  
**Required review:** Human and GPU/dipole owner  
**Risk:** High; not initial scope

#### Prompt

> Produce a design for adaptive coarse graining and FFT dipole deposition when
> the input is an explicit device geometry with `NA=Natom`. Do not force every
> atom to become a basis or dynamical channel. Define coordinate-based spatial
> binning, explicit magnetic-channel labels, occupancy handling, FFT source
> deposition/interpolation, boundary semantics, and capability validation.
>
> Compare at least: regular bins with nonuniform populations, particle-mesh
> deposition to a fixed grid, and rejection/fallback behavior. State which
> existing FFT kernel assumptions must change and how block topology remains
> independent of crystallographic `NA`. This task is a design and prototype
> study, not authorization to modify the production FFT solver.

#### Checklist

- [ ] Memory does not scale as `Natom` channels by construction.
- [ ] Spatial binning is defined for skew cells.
- [ ] Channel labels are independent of basis numbering.
- [ ] Empty and nonuniform cells have explicit semantics.
- [ ] Dipole deposition and adjoint interpolation are paired.
- [ ] Boundary behavior is specified.
- [ ] Compatibility/rejection matrix is explicit.
- [ ] No production FFT behavior changes without a separate accepted plan.

---

### Task CG-13: Release validation, performance, and documentation

**Dependencies:** Each phase independently  
**Suggested primary:** Luna/Sonnet for harness/docs; Sol for performance  
**Required review:** Human  
**Risk:** Medium

#### Prompt

> Assemble the accepted coarse-graining tests into repeatable CPU, CUDA, and
> HIP validation suites. Add user documentation describing model scope,
> inputs, characteristic-length requirements, supported interactions,
> diagnostics, restart limitations, and known sources of interface and dipole
> error. Provide performance benchmarks reporting active fractions and all
> overheads, not only total speedup.
>
> Include reference cases for spin spirals, domain walls, skyrmions, static
> interfaces, adaptive transitions, and multi-channel systems when enabled.
> Ensure every unsupported combination has a tested setup-time error. Do not
> describe low-temperature cone reconstruction as thermodynamic
> equilibration.

#### Checklist

- [ ] Validation suite is documented and automated.
- [ ] Reference data provenance is recorded.
- [ ] CPU/CUDA/HIP matrices are explicit.
- [ ] Precision-specific tolerances are explicit.
- [ ] User input documentation is complete.
- [ ] Model limitations are prominent.
- [ ] Restart support or rejection is documented.
- [ ] Performance reports include active fractions and overhead.
- [ ] No finite-temperature claim exceeds the implemented model.
- [ ] Human approves release wording.

---

## 12. Cross-task review checklists

### Physics review

- [ ] Is there a written energy for every field?
- [ ] Do signs and factors match UppASD Hamiltonian conventions?
- [ ] Is the block-size dependence physically correct?
- [ ] Are DMI spin and spatial tensor indices unambiguous?
- [ ] Are compensated systems represented by channels rather than a net
      macrospin?
- [ ] Are gyromagnetic ratio and damping reductions justified?
- [ ] Is the FFT dipole approximation stated correctly?
- [ ] Are high-frequency/interface limitations acknowledged?
- [ ] Is finite temperature explicitly outside the initial claim?

### Numerical review

- [ ] Does every operator pass tangent energy derivatives?
- [ ] Does the metric support nonorthogonal cells?
- [ ] Are boundaries and halos correct?
- [ ] Is energy ownership non-overlapping?
- [ ] Are transitions synchronized with the integrator?
- [ ] Are RNG results reproducible?
- [ ] Are error thresholds and precision budgets stated?
- [ ] Does error decrease under spatial refinement?

### Software review

- [ ] Is the feature-off path unchanged?
- [ ] Are immutable topology and mutable state separate?
- [ ] Are spatial, basis, FFT, and dynamical channel counts separate?
- [ ] Are unsupported geometries rejected early?
- [ ] Are public interfaces narrow and testable?
- [ ] Are CPU and GPU descriptors identical in meaning?
- [ ] Are all allocations tracked and released?
- [ ] Are compact GPU work lists rebuilt safely?
- [ ] Are unrelated workflows untouched?
- [ ] Is μASD separation preserved except for approved constellation cleanup?

### Pull-request review

- [ ] The PR implements one blueprint task or a clearly declared subset.
- [ ] The prompt, dependencies, and acceptance evidence are linked.
- [ ] New tests fail without the change and pass with it.
- [ ] Baseline failures are distinguished from introduced failures.
- [ ] No unrelated formatting hides physics changes.
- [ ] Sign/unit changes receive explicit physics review.
- [ ] Performance claims include methodology and raw measurements.
- [ ] Follow-up limitations are recorded rather than hidden in TODO comments.

---

## 13. Risk register

| Risk | Consequence | Mitigation |
|---|---|---|
| Existing stiffness outputs are consumed without audit | Incorrect coarse exchange, DMI, or units | CG-03 is a gated refactor with analytic and small-\(q\) validation |
| `NA` is treated as channel count everywhere | Device geometries become impossible; excessive channels | Separate basis, FFT, and dynamic channel mappings |
| Rigid block projection is used for exchange | Wrong block-size scaling and artificial interface energy | Tensor operator plus smooth projection reference |
| Atomistic/coarse energies overlap | Double-counted fields and spurious pinning | Explicit energy ownership and derivative tests |
| Sharp interface reflects spin waves | Corrupted dynamics | Atomistic buffer, smooth ghost projection, documented cutoff |
| Dynamic mask chatters | Noise, overhead, energy injection | Hysteresis, dwell time, filtered/modular criteria |
| Refinement injects uncontrolled modes | Texture destabilization | Aligned default, constrained cone, energy-jump gate |
| Net macrospin is used for AFM/ferri | Singular normalization and lost dynamics | Multi-channel topology from day one |
| FFT blocks are too large | Texture tails or magnetostatics under-resolved | Characteristic-length warnings and convergence tests |
| GPU branching removes speedup | Poor performance despite fewer logical DOFs | Compact active lists and measured crossover |
| Full atom arrays remain allocated | Little memory saving | State this openly; optimize compute first |
| Explicit device geometry is accidentally accepted | Invalid block/channel layout | Geometry capability gate and tested rejection |
| Constellation cleanup perturbs μASD | Regression in protected package | Only approved mechanical removals and reference tests |
| Scope expands to finite temperature prematurely | Scientifically weak claims | Explicit phase boundary and human approval |

---

## 14. Initial definition of success

The first scientifically useful release does not need every outlook item. It
is successful if it can demonstrate all of the following:

1. The obsolete constellation implementation has been removed cleanly.
2. A regular FFT block grid is the single spatial topology for dipole and
   adaptive dynamics.
3. The topology supports multiple channels even though the first runtime
   material is ferromagnetic.
4. The existing stiffness/spiralization physics has been refactored and
   independently validated.
5. All-fine runs reproduce normal ASD.
6. All-coarse smooth textures reproduce atomistic long-wavelength energies
   and dynamics within a documented discretization error.
7. A static atomistic region can coexist with a coarse host without
   unacceptable ghost torque or pinning.
8. A simple angular criterion can move that region while preserving a domain
   wall or skyrmion.
9. The GPU implementation matches the CPU reference and shows a measured
   speedup when most blocks are coarse.
10. Unsupported finite-temperature, explicit-device, and multi-channel cases
    reject cleanly until their own acceptance gates are met.

That result would be a useful, honest, and extensible two-scale capability. It
would preserve UppASD's atomistic strengths, reuse the new FFT block
infrastructure, and leave a credible path to ferrimagnetic,
antiferromagnetic, explicit-device, and thermal extensions without entangling
the new package with μASD.

---

## 15. Selected literature

1. M. Pajda et al., “Ab initio calculations of exchange interactions,
   spin-wave stiffness constants, and Curie temperatures of Fe, Co, and Ni,”
   [Phys. Rev. B 64, 174402 (2001)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.64.174402).
2. C. Andreas, A. Kákay, and R. Hertel, “Multiscale and multimodel simulation
   of Bloch-point dynamics,”
   [Phys. Rev. B 89, 134403 (2014)](https://juser.fz-juelich.de/record/201878/files/PhysRevB.89.134403-1.pdf).
3. A. De Lucia et al., “Multiscale model approach for magnetization dynamics
   simulations,”
   [Phys. Rev. B 94, 184415 (2016)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.184415).
4. M. Poluektov, O. Eriksson, and G. Kreiss, “Coupling atomistic and continuum
   modelling of magnetism,”
   [Comput. Methods Appl. Mech. Eng. 329, 219–253 (2018)](https://www.sciencedirect.com/science/article/abs/pii/S0045782517302463).
5. D. Arjmand, M. Poluektov, and G. Kreiss, “Atomistic-continuum multiscale
   modelling of magnetisation dynamics at non-zero temperature,”
   [Adv. Comput. Math. 44, 1119–1151 (2018)](https://link.springer.com/article/10.1007/s10444-017-9575-3).
6. S. Mankovsky, S. Polesya, and H. Ebert, “Spin-wave stiffness and
   micromagnetic exchange interactions expressed by means of the KKR Green
   function approach,”
   [arXiv:1810.13175](https://arxiv.org/abs/1810.13175).
7. S. Grytsiuk et al., “Ab initio analysis of magnetic properties of the
   prototype B20 chiral magnet FeGe,”
   [Phys. Rev. B 100, 214406 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.214406).
8. M. Gija et al., “Multiscale micromagnetic / atomistic modeling of heat
   assisted magnetic recording,”
   [arXiv:2502.02236 (2025)](https://arxiv.org/abs/2502.02236).
