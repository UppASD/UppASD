# CG-01: initial adaptive coarse-graining physics conventions

**Status:** approved implementation contract.  Human physics approval and the
descriptor implementability review are recorded at the end.  This document
freezes the conventions used by CG-02 onward; it is not a claim that the
operators exist or that every current atomistic backend already conforms.

## 1. Canonical atomic convention and units

Let `emomM(:,i)` be the atomistic magnetic moment in units of the Bohr
magneton, `mmom(i) s_i`, where \(\lvert\mathbf s_i\rvert=1\).  The active
Hamiltonian field path defines `beff` in tesla.  Its pair conventions are

\[
 E_{J}=-\frac{\mu_B}{2}\sum_{i,j}\mathbf M_i\!\cdot
 J_{ij}\mathbf M_j,
 \qquad
 E_{D}=\frac{\mu_B}{2}\sum_{i,j}\mathbf D_{ij}\!\cdot
 (\mathbf M_i\times\mathbf M_j),
\]

with directed lists, \(J_{ij}=J_{ji}\) and \(\mathbf D_{ji}=-\mathbf
D_{ij}\).  Thus

\[
 \mathbf B_i=-\frac{1}{\mu_B}\frac{\partial E}{\partial\mathbf M_i}
 =\sum_j J_{ij}\mathbf M_j+\sum_j\mathbf D_{ij}\times\mathbf M_j.
\]

This is the approved public convention.  The scalar `heisenberg_field`
implements its exchange part.  `Energy:update_ene` reports
\(-f\,\mathbf M\cdot\mathbf B\): use \(f=1/2\) for a reciprocal pair
term, dipole, or DMI, and \(f=1\) for Zeeman.  The surrounding energy code
subsequently multiplies by `mub/mry`, so an energy written in
field-times-\(\mu_B\) units is converted to mRy by
\(\mu_B/\mathrm{mRy}\).  A coarse implementation shall use SI joules
internally for material energies and tesla for fields; it shall *not* insert
`mub/mry` in a field kernel.

There is one pre-existing atomistic conformance issue that later tasks must
not mistake for a convention choice.  DMI input is copied into `dm_vect`
without a sign change, and the GPU `hamdev::dm_field`, GPU lattice-convolution
DMI tensor, and `applyhamiltonian` evaluate
\(\mathbf D_{ij}\times\mathbf M_j\), as specified above.  The current CPU
`HamiltonianActions:dzyaloshinskii_moriya_field` evaluates the opposite cross
product, \(\mathbf M_j\times\mathbf D_{ij}\).  The dimer fixture in section 4
is authoritative; CPU/GPU conformance must be fixed and tested before a
coarse DMI operator is accepted.  CG-01 itself does not alter production
kernels.

The field identity for every coarse term is consequently

\[
 \boxed{\quad\mathbf B_{b\alpha}=-\frac{1}{\mu_B\mu_{b\alpha}}
 \frac{\partial E}{\partial\mathbf m_{b\alpha}}\quad}                 \tag{1}
\]

where \(\mu_{b\alpha}\) is in units of \(\mu_B\), \(E\) is in joules,
and \(\mathbf m\) is a unit direction.  Tangent projection belongs to the
LLG update, not to the energy derivative.  Equation (1) gives the atomistic
formula above when a channel contains one atom.

`gama` is the positive magnitude \(\gamma_0=1.760859644\times10^{11}\)
s\(^{-1}\)T\(^{-1}\); the existing deterministic LLG code multiplies it by
the per-site `Landeg` and uses Gilbert damping `lambda1_array`.  The initial
coarse path shall use the same sign and Gilbert form,

\[
 \dot{\mathbf m}=
 -\frac{\gamma_0 g}{1+\lambda^2}
 [\mathbf m\times\mathbf B+\lambda\,\mathbf m\times(\mathbf m\times\mathbf B)]. \tag{2}
\]

## 2. State, channels, and dynamics

For spatial block \(b\), dynamical channel \(\alpha\), and ensemble \(e\),
restriction stores

\[
 \mathbf R_{b\alpha e}=\sum_{i\in(b,\alpha)}\mu_i\mathbf s_{ie},\qquad
 \mu_{b\alpha}=\sum_{i\in(b,\alpha)}\mu_i,\qquad
 p_{b\alpha}=|\mathbf R_{b\alpha}|/\mu_{b\alpha}.
\]

The active coarse direction is \(\mathbf m=\mathbf R/|\mathbf R|\) only
after the configured polarization criterion accepts the channel.  A zero or
insufficient resultant is never normalized.  \(\mu_{b\alpha}\) is the sum
of moment magnitudes, not \(|\mathbf R|\), and is held fixed in the initial
fixed-length model.

Initial enablement requires all atoms in a channel to have exactly the same
`Landeg` and `lambda1_array`; then `g` and \(\lambda\) are those common
values.  Heterogeneous Landé factors and damping are **rejected at setup**.
The natural angular-momentum result
\(g_\mathrm{eff}=\sum_i\mu_i/\sum_i(\mu_i/g_i)\) is recorded for future
work, but it does not authorize heterogeneous damping: its dissipative
reduction and fluctuation--dissipation partner need a separate derivation and
approval.  Thermal noise is not enabled.

## 3. Coordinates, indexing, and discrete energy

All spin components \(k,l\) and physical spatial components \(p,q\) are
Cartesian and range over \(x,y,z\).  Channel indices are \(\alpha,\beta\).
If the physical block vectors are columns of
\(L=[\mathbf B_1\ \mathbf B_2\ \mathbf B_3]\), then
\(\mathbf x=\mathbf x_0+L\boldsymbol\xi\),
\(V=|\det L|\), and

\[
 \partial_{x_p}= (L^{-T})_{pr}\partial_{\xi_r},\qquad
 G_p=(L^{-T})_{pr}\Delta_r.                                      \tag{3}
\]

`Delta_r` is the selected finite-volume/fixed-difference derivative in the
dimensionless block coordinate.  It and its exact transpose (including
boundary weights) are stored as an operator pair.  No kernel may treat a
skew block index as a Cartesian direction.  For uniform blocks, the frozen
all-coarse energy is

\[
\begin{aligned}
 E_\mathrm{xc}&=\sum_b V\sum_{\alpha\beta pq}
 A^{\alpha\beta}_{pq}(G_p\mathbf m_{b\alpha})\cdot(G_q\mathbf m_{b\beta}),\\
 E_\mathrm{D}&=\sum_b V\sum_{\alpha\beta kp}
 D^{\alpha\beta}_{kp}\,[\mathbf m_{b\alpha}\times
 G_p\mathbf m_{b\beta}]_k,\\
 E_\mathrm{ani}&=\sum_{b\alpha} V\,\varepsilon_{\mathrm{ani},b\alpha}
 (\mathbf m_{b\alpha}),\\
 E_\mathrm{Z}&=-\mu_B\sum_{b\alpha}\mu_{b\alpha}
 \mathbf B_\mathrm{ext,b\alpha}\cdot\mathbf m_{b\alpha},\\
 E_\mathrm{dip}&=E_\mathrm{FFT\ dipole}(\{\mu_{b\alpha}\mathbf m_{b\alpha}\}
 \mathbin{\cup}\{\mathbf M_i\}_\mathrm{fine}).
\end{aligned}                                                     \tag{4}
\]

`exchange_stiffness(p,q,alpha,beta)` means exactly
\(A^{\alpha\beta}_{pq}\): first/second Cartesian derivative index,
then target/source channel in the displayed bilinear form.  It must satisfy
\(A^{\alpha\beta}_{pq}=A^{\beta\alpha}_{qp}\).  `spiralization(k,p,alpha,beta)`
means exactly \(D^{\alpha\beta}_{kp}\): first Cartesian spin component of
the cross product, then Cartesian derivative direction, then left/right
channel.  Descriptor metadata must name this ordering and its SI units
(\(A\): J m\(^{-1}\), \(D\): J m\(^{-2}\)).

With \(G_p^T\) denoting the Euclidean transpose and (V\) retained explicitly,
the exact unconstrained derivatives of (4) are

\[
\begin{aligned}
 \frac{\partial E_\mathrm{xc}}{\partial\mathbf m_\gamma}
 &=2\sum_{\beta pq}G_p^T\!\left[V A^{\gamma\beta}_{pq}
 G_q\mathbf m_\beta\right],\\
 \frac{\partial E_\mathrm{D}}{\partial\mathbf m_\gamma}
 &=\sum_{\beta kp}\{-V D^{\gamma\beta}_{kp}
 (\mathbf e_k\times G_p\mathbf m_\beta)
 +G_p^T[V D^{\beta\gamma}_{kp}(\mathbf e_k\times\mathbf m_\beta)]\},\\
 \frac{\partial E_\mathrm{Z}}{\partial\mathbf m_{b\alpha}}
 &=-\mu_B\mu_{b\alpha}\mathbf B_\mathrm{ext,b\alpha}.
\end{aligned}                                                     \tag{5}
\]

Equation (1) applied to (5), plus the derivative supplied by each anisotropy
and FFT-dipole operator, is the only field contract.  For constant one-channel
coefficients and periodic central differences, \(G_p^T=-G_p\), so
\(\mathbf B_D=2D_{kp}\mathbf e_k\times G_p\mathbf m/(\mu_B\mu)\).
The factor two is physical: the DMI density in (4) contains no hidden 1/2.

Initially supported anisotropy is only a descriptor whose energy and
derivative are supplied together.  The first required descriptor is the
UppASD-compatible uniaxial form (one or two axes): for
\(c=\mathbf m\cdot\mathbf e\),
\[
 \varepsilon_\mathrm{ani}=K_1c^2+2K_2c^2-K_2c^4,
 \quad \partial\varepsilon/\partial\mathbf m=
 2c[K_1+2K_2(1-c^2)]\mathbf e.                                   \tag{6}
\]
The corresponding field is minus (1) times this derivative.  This matches
the current `uniaxial_anisotropy_field`; `K1` and `K2` must be converted from
the current field-times-\(\mu_B\) representation to J m\(^{-3}\) before
use.  Cubic, tensor exchange, pseudo-dipolar, biquadratic, ring, chiral,
LSF/induced-moment, and torque terms have no initial coarse descriptor and
are rejected whenever a coarse block would need them.

## 4. Analytic sign and prefactor fixtures

These fixtures are mandatory unit tests before an operator is accepted.

| fixture | energy / expected result |
| --- | --- |
| reciprocal Heisenberg dimer | \(E=-\mu_BJ\mathbf M_1\cdot\mathbf M_2\), \(\mathbf B_1=J\mathbf M_2\); confirms no extra pair factor |
| DMI dimer | \(E=\mu_BD\,\hat z\cdot(\mathbf M_1\times\mathbf M_2)\), \(\mathbf B_1=D\hat z\times\mathbf M_2\); fixes the active UppASD handedness |
| periodic DMI chain | \(\mathbf m=(\cos qx,\sin qx,0)\) gives \(E/V=Aq^2+D_{zx}q\), hence \(q_\min=-D_{zx}/(2A)\) |
| one coarse Zeeman channel | \(E=-\mu_B\mu\mathbf B\cdot\mathbf m\) gives exactly \(\mathbf B\) through (1) |
| uniform state | exchange is zero; DMI energy/torque is zero for periodic constant \(\mathbf m\) |
| skew affine spiral | evaluate (3) for a linear phase; its physical \(\nabla\phi\), not fractional slope, must enter (4) |

Each test also performs a central tangent finite difference of the energy
against \(-\mu_B\mu\,\mathbf B\cdot\mathbf t\).  The DMI dimer is the
authoritative handedness fixture; input-vector mapping may not silently flip
it.

## 5. Non-overlapping energy ownership

For a static resolution mask, form a buffer by dilating atomistic blocks by
the maximum short-range atomistic interaction radius, rounded outward to
blocks, plus the configured safety layer.  The total energy is

\[
 E=E_\mathrm{atom}[F\cup I;\,\text{real fine spins and coarse ghosts}]
   +E_\mathrm{coarse}[C\setminus I]
   +E_\mathrm{dip,FFT}[\text{all blocks}],                         \tag{7}
\]

where \(F\) is atomistic interior, \(I\) the buffer/interface, and \(C\)
the coarse set.  `E_atom` owns every short-range atomistic interaction with
at least one real fine/buffer atom; ghost spins only provide its boundary
arguments.  Its reaction field is restricted with the derivative-adjoint of
the prolongation.  `E_coarse` owns only stencils wholly in coarse interior;
all coarse terms touching \(I\) are suppressed.  FFT dipole is evaluated
once on its fixed regular grid, receives all block moments, and is never
added through either short-range owner.  This avoids both pair double counts
and a missing interface reaction torque.

The initial interface uses smooth trilinear prolongation in fractional block
coordinates followed by normalization.  Its field restriction must use the
derivative of that normalization; a linear-weight transpose alone is not an
energy derivative.  A sharp or piecewise-constant coarse exchange interface
is diagnostic-only.

## 6. Initial capability matrix

| category | accepted initially | rejected initially |
| --- | --- | --- |
| geometry | complete divisible `REGULAR_REPLICATED_CELL`, including skew physical cells | `NA=Natom` explicit-device layouts, partial blocks, arbitrary populations |
| channels | one ferromagnetic dynamical channel per block; channel dimension allocated everywhere | multi-channel dynamics, AFM/ferri coarse physics, zero-resultant normalization |
| short-range Hamiltonian | scalar Heisenberg, DMI, uniaxial descriptor after the fixtures pass | tensor exchange, SA/PD/BQ/BIQDM/ring/chiral, LSF, induced moments |
| long range / applied | existing uniform regular-grid FFT dipole; static external field | adaptive FFT, local dipole correction, unspecified time-dependent field coupling |
| dynamics | deterministic fixed-length Gilbert LLG, common global step, existing deterministic solver family | MC, GNEB, SLD, \(\mu\)ASD evolution, STT/SHE/SOT, subcycling |
| temperature | \(T=0\), no stochastic field | finite temperature, QHB, temperature-dependent material reduction |
| boundaries | periodic block derivatives; open boundaries only once the selected FFT mode and adjoint one-sided stencil are jointly validated | mixed/implicit boundary behavior, unvalidated open stencils |
| mask/restart | static mask; restart explicitly rejected until state serialization exists | adaptive transitions and restart claiming production support |

Feature-off behavior remains exactly atomistic.  An enabled input outside this
matrix must fail during setup with the violated capability named.

## 7. Multi-channel descriptor and storage contract

Topology and runtime must keep spatial blocks, basis sites, FFT source
channels, and dynamical channels distinct.  A conforming descriptor contains
the following logical fields; an implementation may use sparse stencils or
opaque handles instead of the dense notation, but the Fortran and C views
must have identical meaning.

```text
n_spatial_blocks, n_basis, n_fft_channels, n_dynamic_channels, n_ensembles
block_shape(3), block_vectors(3,3), block_volume(block), boundary_mode(3)
derivative_operator(physical_space,...), derivative_transpose(physical_space,...)

atom_to_block(atom), atom_to_basis(atom), atom_to_fft_channel(atom)
atom_to_dynamic_channel(atom)        ! -1 only for nonmagnetic atoms
block_channel_population(channel, block)
moment(3, channel, block, ensemble)  ! R before normalization
moment_sum(channel, block)           ! mu in muB
direction(3, channel, block, ensemble), field(3, channel, block, ensemble)
g_factor(channel, block), damping(channel, block), polarization(channel, block)
static_external_field(3, channel, block, ensemble)

fft_channel_weight(fft_channel, channel, block) ! muB mapped to each FFT source
fft_moment(3, fft_channel, block, ensemble)     ! deposited source in muB

exchange_stiffness(3,3,channel,channel)
spiralization(3,3,channel,channel)
local_exchange(channel,channel)

anisotropy_kind(channel,block)
anisotropy_axis_count(channel,block)            ! 0, 1, or 2
anisotropy_axis(3,2,channel,block)
anisotropy_K1(2,channel,block), anisotropy_K2(2,channel,block)

resolution_mask(block), energy_owner(block), channel_capability_mask(channel)
ordering_enums, units_enums, dmi_energy_sign, field_derivative_sign
convention_version
```

`block_vectors(:,r)` is the physical vector for one block step in fractional
direction \(r\), in metres.  `derivative_operator` is \(G_p\) from (3), and
`derivative_transpose` is its exact discrete transpose, including boundary
weights.  Periodic boundaries are the only initially enabled value.

`exchange_stiffness` has ordering `(derivative_p, derivative_q, left_channel,
right_channel)` and units J m\(^{-1}\).  `spiralization` has ordering
`(cross_product_spin_k, derivative_p, left_channel, right_channel)` and units
J m\(^{-2}\).  `local_exchange` is the coefficient \(\Lambda_{\alpha\beta}\)
of an energy density and has units J m\(^{-3}\).  The initial uniaxial
anisotropy kind stores one or two normalized Cartesian axes and one `(K1,K2)`
pair per axis, also in J m\(^{-3}\), with the energy and derivative fixed by
(6).  `static_external_field` and `field` are in tesla.  Moment sums, raw
moments, FFT weights, and FFT moments are numerical multiples of \(\mu_B\);
directions, polarization, \(g\), and damping are dimensionless.
`dmi_energy_sign` is the plus sign in (4), and `field_derivative_sign` is the
minus sign in (1); these are validated enum values rather than comments.

FFT sources are accumulated independently from dynamical channels while
fine spins exist.  In a coarse block they obey

\[
 \mathbf M^{\rm FFT}_{bf e} =
 \sum_\alpha w_{bf\alpha}\mathbf m_{b\alpha e},
\]

where `fft_channel_weight` is derived from `atom_to_fft_channel` and
`atom_to_dynamic_channel`.  This explicit weighted mapping permits several
dynamical channels to contribute to one FFT source and one dynamical channel
to contribute to several FFT sources.  It prevents an FFT channel from being
inferred from a basis or dynamical-channel index.

The implementability review against the approved capability matrix is:

| approved item | descriptor representation |
| --- | --- |
| scalar Heisenberg / continuum exchange | `exchange_stiffness`; channel-ready `local_exchange`; exact \(G_p,G_p^T,V\) |
| DMI | `spiralization` plus ordering and convention metadata |
| one- or two-axis uniaxial anisotropy | explicit kind, axis count, axes, `K1`, and `K2`; energy/derivative pair fixed by (6) |
| static external field | `static_external_field`, separate from material tensors |
| regular-grid FFT dipole | independent atom-to-FFT map, channel weights, and `fft_moment` |
| fixed-length deterministic LLG | `moment_sum`, `direction`, `field`, common `g_factor`, and `damping` |
| skew regular geometry / periodic derivatives | physical block vectors, volume, boundary enum, and adjoint operator pair |
| static fine/interface/coarse ownership | resolution mask and explicit energy owner |

Every initially approved term therefore has the state, coefficient, geometry,
units, ordering, and convention metadata needed by both CPU and GPU
implementations.  Deferred multi-channel physics is representable without
collapsing a compensated material to a net macrospin.  No code may infer a
dynamic channel from basis number, normalize a zero resultant, or enable a
term whose capability bit is absent.

## 8. Approved decisions and implementation gates

1. The DMI energy sign and \(q_\min=-D/(2A)\) handedness fixture are the public
   convention.  The no-sign-change input mapping and the source conformance
   issue are recorded in section 1.
2. Initial coarse support is limited to the matrix above, including deferral
   of cubic anisotropy and time-dependent fields.
3. Only periodic derivatives are enabled initially.  An open-boundary/FFT
   pairing remains a later, separately reviewed extension.
4. Polarization and resolution safety thresholds remain explicit
   configuration/acceptance parameters; no undocumented universal value is
   implied by this contract.
5. Heterogeneous \(g\), heterogeneous damping, and finite-temperature
   reduction require a separate angular-momentum/noise derivation and human
   approval.

**Required sign-off:** Human physics owner approval of sections 1--6; Sol
review that the section 7 descriptor can represent every accepted term.

**Sign-off:** Contract approved by human review.

**Implementability sign-off:** Descriptor review complete.  The table in
section 7 confirms representation of every initially approved term; the DMI
backend mismatch in section 1 is a required production-conformance test, not
a missing descriptor.
