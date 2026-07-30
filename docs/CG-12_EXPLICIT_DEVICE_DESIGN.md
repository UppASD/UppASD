# CG-12: explicit-device adaptive coarse graining and FFT dipole design

## Decision

This is a design/prototype-study result. It does **not** change a production
FFT solver, GPU dispatch, input selector, or the current explicit-device
rejection. That rejection is correct: an explicit device commonly has
`NA=Natom`; treating it as an FFT or dynamical-channel count would allocate an
`O(Natom^2)` kernel contract and accidentally make atom numbering physics.

The future model has three separate identities:

```text
atom i  --coordinate bin--> spatial block b
atom i  --declared label--> magnetic channel c
source entity e=(fine atom i | coarse (b,c)) --assignment--> mesh nodes g
```

`NA` remains geometry provenance only. It is not a device-topology channel
count or an FFT batch dimension. The physical dipole mesh stores one vector
source per mesh node: all channel contributions are summed before the dipole
solve. Magnetic channels remain a local material/dynamics concept.

## Independent device topology

Introduce a later `DEVICE_COORDINATE_BINS` topology rather than weakening
`BlockTopology`'s `REGULAR_REPLICATED_CELL` contract. Its immutable contents
are an origin `r0`, a nonsingular frame `A`, bin grid `B=(B1,B2,B3)`, boundary
mode, declared magnetic labels `1..C`, atom-to-block/channel maps, CSR block
membership and nonempty `(block,channel)` groups, and occupancy/weight
diagnostics. It contains no spin state, adaptive mask, field, or FFT buffers,
and it retains the input atom numbering.

Labels must be supplied by a device channel file or a reviewed atom-property
table. They can represent a species or a ferrimagnetic sublattice, but cannot
be inferred from atom index, position-file row, or `NA`. Label `0` is
nonmagnetic; every positive label needs declared material parameters. `C` is
small and physical, independent of `NA`. Allocate local channel state only
for nonempty magnetic pairs using CSR (or a checked dense `B*C` array). Thus
the design does not construct `Natom` channels. A user can nevertheless ask
for one spatial bin per atom; that requested resolution is caught by explicit
memory/capability limits, not disguised as a basis map.

### Coordinate binning, including skew cells

Use Cartesian coordinates in one consistent physical unit (after the existing
`alat*coord` conversion). Let the columns of `A` span the full device frame.
For every atom:

```text
u_i = A^-1 (r_i-r0)
```

This fractional coordinate is the only binning coordinate: component-wise
Cartesian bins are prohibited because they are wrong for skew cells. For
periodic geometry set `A=H`, the full supercell matrix, and wrap
`u <- u-floor(u)`. For open geometry require `0 <= u_j < 1`, except for a
documented seam tolerance; a canonical upper physical boundary is assigned to
the final bin. An out-of-frame point is an error, never silently clamped.
Half-open bins are:

```text
b_j = min(B_j-1, floor(B_j*u_j))
r_b = r0 + A ((b+(1/2,1/2,1/2))/B)
V_b = abs(det(A))/(B1*B2*B3).
```

Persist the exact origin and frame across restart. Rebuilding an axis-aligned
bounding box can change seam membership.

### Occupancy and coarse variables

For fixed scalar moment `mu_i` and direction `s_i`, define

```text
W_bc = sum_{i in (b,c)} mu_i
M_bc = sum_{i in (b,c)} mu_i s_i.
```

Population is a count, never a normalisation. `W_bc` is the restriction/field
weight. Empty bins and empty pairs have zero source, no dynamic state, and no
division. A nonempty pair with insufficient `W_bc` is rejected unless an
explicit material model accepts it. Nonuniform populations are valid metadata
and must not be replaced by a neighbouring population. A compensated vector
`M_bc=0` is valid but has no single-spin direction; multi-channel dynamics
must address that explicitly.

For a coarse source, select one immutable reconstruction policy per run:

1. Distributed members (recommended): put the coarse vector back on its
   member positions with static weights `mu_i/W_bc`. This preserves total
   moment and within-bin geometry, without creating more dynamical channels.
2. Centroid macrospin: place `M_bc` at a fixed moment-weighted geometric
   centroid. This is a controlled monopole approximation; it loses higher
   multipoles and requires a bin-size/error gate. The centroid cannot depend
   on the current spin direction.

Fine atoms and reconstructed coarse entities are disjoint source entities, so
mixed adaptive regions deposit each physical moment exactly once.

## Fixed-grid particle--mesh deposition

The general explicit-device FFT approach is a dipolar particle--mesh operator
on a separately selected mesh `G=(G1,G2,G3)`, frame `A_mesh`, and origin
`r0_mesh`. It need not equal the bin grid. The initial prototype assignment
can be CIC, later a cardinal B-spline of declared order. For source `e` at
`x_e`, form compact weights `D_ge=w_g(x_e)` with nonnegative weights and
`sum_g w_g=1` to tolerance. Then use the exact adjoint pair:

```text
M_g    = sum_e D_ge m_e
B_g    = K_mesh * M_g                 (FFT convolution)
B_e    = sum_g D_ge B_g = (D^T B)_e
E      = -1/2 sum_g M_g . B_g = -1/2 sum_e m_e . B_e.
```

The same weights, positions, frame, boundary wrapping, and precision policy
must be used in both directions. For distributed coarse sources, average
member fields with `mu_i/W_bc` to obtain the channel field. Empty mesh nodes
are ordinary zero FFT sources; atom/bin occupancy is never used to divide a
mesh moment.

A proper dipolar P3M/PME kernel is required: assignment-order-aware influence
function, mesh/self correction, and a declared near-field correction policy.
Applying the current point-source lattice Green function after arbitrary CIC
deposition is not valid: it neither removes assignment smearing nor mesh
aliases.

## Boundary semantics

| Boundary | Coordinates and stencils | Field definition | Status |
|---|---|---|---|
| `P P P` | `A_mesh=H`; wrap fractional coordinates and any seam-crossing stencil | Periodic P3M/Ewald; full `H`; tin-foil physical `k=0` omission unless another reviewed surface term is selected | Future prototype |
| `0 0 0` | Never wrap; source and its assignment support must fit the frame plus guard band | Finite/open field; use zero-padded linear convolution, at least `2G_j-1` each axis, with physical skew-frame displacements | Future prototype |
| mixed/slab | No inferred wrapping choice | Requires separately derived mixed/slab Green function and surface convention | Reject |
| moving positions/frame | Membership/weights become stale | Requires an explicit rebuild/restart policy | Reject initially |

An FFT is only an accelerator: open zero-padding must prevent periodic images.
There is no automatic boundary choice based on device shape or vacuum.

## Options compared

| Option | Benefits | Limits | Recommendation |
|---|---|---|---|
| Regular coordinate bins with nonuniform populations | Small source grid; holes and variable occupancy are simple | Point at bin centre discards off-centre multipoles. Channels must be summed before the dipole solve because occupancy-dependent channel centroids are not a translationally invariant basis | Controlled macrospin prototype only, with bin-size/error evidence; not atomistically exact |
| Particle--mesh to a fixed grid | Arbitrary positions, holes, labels, and population without `NA` FFT channels | Requires new P3M influence, self/near correction, adjoint operators, and independent oracles | Recommended general explicit-device FFT design |
| Rejection/fallback | No silent physics change | No acceleration until a new mode passes review | Required default; unsupported pairs reject and users explicitly select an existing atomistic/direct method or disable CG |

The regular-bin model is a distinct macrospin Hamiltonian, not a shortcut for
particle mesh. It may use a new single-channel regular-grid kernel, but cannot
reuse the current `basis=NA` descriptor unchanged.

## Existing FFT assumptions that must change only in a new path

| Current assumption | Existing anchor | Future device replacement |
|---|---|---|
| `macro_count = active_cells*basis`, with `basis=NA` | `GpuDipoleConvolutionDescriptor::valid()` | One physical mesh source; labels exist only in topology CSR |
| Coarse grid divides `N1,N2,N3`; every basis macro has equal population | GPU Hamiltonian setup and periodic builder | Coordinate membership plus assignment validation; keep present checks for present modes |
| Sources lie at a cell plus a global basis offset | periodic/open displacement builders | Arbitrary positions through `D`; centre positions only for the bin-macrospin prototype |
| `9*NA*NA` kernel batches and basis block contraction | `apply_spectral_kernel` | Nine Cartesian mesh tensor batches; sum labels before deposition |
| One-based basis-fast macro pack/scatter | `pack_macro_moments_kernel`, `add_dipole_fields` | CSR source list, deposit, and adjoint interpolation kernels |
| Complete lattice Green function is exact for point sources on its regular macro lattice | `buildPeriodicKernel()` | Assignment-aware P3M influence plus self/near correction, retaining `H`, units, `k=0`, and FFT normalization conventions |
| Open grid is a regular macro lattice at origin | OPEN_FFT pack/builder | Affine mesh embedding, guarded sources, paired interpolation, and retained no-wrap proof |

The present complete periodic lattice Green function must not be renamed PME.
Its `FFT(K)/Ngrid` plus raw-inverse normalisation can remain, but not its
point-source/basis-grid source model.

## Capability validation

Before allocation, a future explicit-device preflight must validate:

1. finite coordinates/moments/frame, nonzero `det(A)`, mesh-count overflow,
   FFT dimensions, and memory budget;
2. exactly one coordinate bin per atom, deterministic seams, and no invalid
   open-frame point;
3. label coverage, material data, nonempty/moment thresholds, and no implicit
   atom-to-channel conversion;
4. finite nonnegative assignment weights, partition of unity, legal node IDs,
   and wrapping only for fully periodic mode;
5. boundary/kernel/precision/prefactor/reconstruction compatibility, rejecting
   instead of silently switching periodic-to-open or mesh-to-bin-centre; and
6. an independent direct reference for each new boundary, mesh order, and
   precision combination before runtime enablement.

| Geometry | Boundary and request | Result until separately accepted |
|---|---|---|
| Regular replicated cell | current `EWALD3D_FFT` / `OPEN_FFT` | Existing paths, unchanged |
| Explicit device, including `NA=Natom` | current FFT selector | Reject before allocation; existing behavior remains |
| Explicit device, coordinate bins | `P P P` or `0 0 0` centre-macrospin | Prototype only after multipole/bin-size oracle acceptance |
| Explicit device, particle mesh | `P P P` | Disabled until periodic P3M, correction, skew-cell, derivative, and energy acceptance |
| Explicit device, particle mesh | `0 0 0` | Disabled until padded affine-mesh/open-reference acceptance |
| Any explicit device | mixed/slab, moving geometry, undefined labels, invalid occupancy/assignment | Reject |

Prototype acceptance tests: skew-frame membership and periodic seam wrapping;
holes/nonuniform-population accounting; label permutation invariance; source
moment conservation; `D^T` energy identity; open no-wrap padding; periodic
translation invariance; two-particle comparison; refinement versus independent
periodic Ewald and finite-open direct references; and energy-versus-rotation
finite differences before adaptive dispatch.

## Recommended progression

1. Keep the explicit-device rejection and regular `BlockTopology` unchanged.
2. Add a host-only `DEVICE_COORDINATE_BINS` prototype with skew/seam/occupancy/
   label tests and no FFT dispatch.
3. Add host/device deposit--adjoint-interpolate tests on synthetic meshes,
   proving source and energy identities independently of an Ewald kernel.
4. Develop periodic and open particle--mesh dipole kernels as a reviewed new
   work package with independent oracles and a documented near-field policy.
5. Only then design nonuniform-block material extraction and request human and
   GPU/dipole-owner acceptance for a new explicit-device selector.

This preserves a topology independent of crystallographic `NA`, bounds
channels by physical labels, and makes every mesh approximation explicit.
