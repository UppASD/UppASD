# CG-01: initial adaptive coarse-graining physics conventions

**Status:** proposed implementation contract; requires the two approvals at
the end before any physics implementation is enabled.  This document freezes
the conventions used by CG-02 onward, not a claim that the operators exist.

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

This is the convention implemented by `heisenberg_field` and
`dzyaloshinskii_moriya_field`.  `Energy:update_ene` reports
\(-f\,\mathbf M\cdot\mathbf B\): use \(f=1/2\) for a reciprocal pair
term, dipole, or DMI, and \(f=1\) for Zeeman.  It subsequently multiplies by
`mub/mry`, so an energy written in field-times-\(\mu_B\) units is converted
to mRy by \(\mu_B/\mathrm{mRy}\).  A coarse implementation shall use SI
joules internally for material energies and tesla for fields; it shall *not*
insert `mub/mry` in a field kernel.

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

## 7. Multi-channel storage contract (allocated before multi-channel physics)

Topology and runtime must keep spatial blocks, basis sites, FFT source
channels, and dynamical channels distinct.  The conceptual arrays are

```text
atom_to_block(atom), atom_to_basis(atom), atom_to_fft_channel(atom)
atom_to_dynamic_channel(atom)        ! -1 only for nonmagnetic atoms
block_channel_population(channel, block)
moment(3, channel, block, ensemble)  ! R before normalization
moment_sum(channel, block)           ! mu in muB
direction(3, channel, block, ensemble), field(3, channel, block, ensemble)
g_factor(channel, block), damping(channel, block), polarization(channel, block)
exchange_stiffness(3,3,channel,channel)
spiralization(3,3,channel,channel)
local_exchange(channel,channel), anisotropy_descriptor(channel,...)
```

Fortran/C descriptors must also carry `n_spatial_blocks`, `n_basis`,
`n_fft_channels`, `n_dynamic_channels`, ordering/units enums, and a
per-channel capability mask.  FFT sources are deposited from their own map;
the total FFT moment is the sum of mapped channel moments only where that map
declares it.  No code may infer a dynamic channel from basis number or use a
net block moment as a compensated-material channel.

## 8. Decisions awaiting human approval

1. Approve the DMI energy sign and the \(q_\min=-D/(2A)\) handedness fixture
   as the public convention, after comparing a real UppASD DMI input vector
   with the dimer fixture.
2. Approve that initial coarse support is limited to the matrix above,
   particularly deferral of cubic anisotropy and time-dependent fields.
3. Select the first approved open-boundary/FFT pairing; until then only
   periodic derivatives are enabled.
4. Approve the polarization threshold and resolution safety policy; these are
   numerical acceptance criteria, not derivable constants.
5. Approve any future heterogeneous \(g\), damping, or finite-temperature
   reduction only with a separate angular-momentum and noise derivation.

**Required sign-off:** Human physics owner approval of sections 1--6; Sol
review that the section 7 descriptor can represent every accepted term.  Both
remain pending at this documentation-only milestone.
