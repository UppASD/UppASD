# CG-03: typed stiffness and spiralization extraction

**Status:** implementation complete and accepted. Focused analytic tests pass,
and human review approved the signs, prefactors, units, conventions, and
runtime gate.

## Scope and compatibility

`source/SpinWaves/stiffness.f90` now exposes a side-effect-light material
service alongside the legacy `ferro_stiffness`, `DMI_stiffness`, and
`ferro_random_stiffness` entry points. The new calculation, fitting,
validation, and runtime-gate routines perform no printing, do not use fixed
file units, and do not store their result in the legacy module globals.

The existing `stiffness_wrapper` reporting path and its legacy numerical
reductions are unchanged. CG-03 therefore makes no separately approved change
to historical stiffness output. In particular, it does not change the
pre-existing CPU atomistic DMI cross-product issue recorded by CG-01.

The direct UppASD adapter currently accepts an ordered regular replicated
cell. Random-alloy output remains on its legacy reporting path; a random-alloy
runtime material requires a separately reviewed averaging/channel model.

## Directed-pair convention and SI lattice sums

The typed input stores Cartesian displacements in metres and directed pair
coefficients in joules:

\[
 E_J=-\frac12\sum_{ij}K_{ij}\,\mathbf m_i\cdot\mathbf m_j,\qquad
 E_D=+\frac12\sum_{ij}\mathbf L_{ij}\cdot
       (\mathbf m_i\times\mathbf m_j).
\]

For a central-cell outgoing pair table with reciprocal partners retained, the
raw continuum coefficients are

\[
 A^{\alpha\beta}_{pq}(\eta)=
 \frac{1}{4V}\sum_{i\in\alpha,j\in\beta}
 K_{ij}r_{ij,p}r_{ij,q}e^{-\eta |r_{ij}|},
\]

\[
 D^{\alpha\beta}_{kp}(\eta)=
 \frac{1}{2V}\sum_{i\in\alpha,j\in\beta}
 L_{ij,k}r_{ij,p}e^{-\eta |r_{ij}|}.
\]

`A` is stored as
`(space_p,space_q,left_channel,right_channel)` in J m\(^{-1}\).
`D` is stored as
`(cross_spin_k,space_p,left_channel,right_channel)` in J m\(^{-2}\).
The inter-channel local exchange is the coefficient of
\(\Lambda_{\alpha\beta}\mathbf m_\alpha\cdot\mathbf m_\beta\) for
\(\alpha<\beta\), in J m\(^{-3}\). It is symmetrized from both directed
channel sectors. No sign is inferred from an initial moment direction.

The UppASD adapter obtains the directed SI coefficients from the active field
Hamiltonian:

```text
K_ij     = mu_B * ncoup_ij * |mu_i| * |mu_j|
L_ij(:)  = mu_B * dm_vect_ij(:) * |mu_i| * |mu_j|
r_ij(:)  = alat * wrapped_coordinate_difference(:)
```

It rejects reduced-unit constants, invalid neighbour/basis indices, and
non-regular central-cell ordering.

## API and result contract

The public types and entry points are:

```text
coarse_lattice_input_type
coarse_lattice_sums_type
coarse_material_type

calculate_coarse_lattice_sums
fit_coarse_material
extract_coarse_material
extract_coarse_material_from_uppasd
validate_coarse_material_small_q
coarse_material_runtime_status
```

`coarse_material_type%raw` retains basis-resolved and channel-aggregated
local exchange, `A`, and `D` for every requested convergence parameter. No
eigenvalue collapse is performed. The fitted channel tensors, channel
moments, channel gamma, damping, units, tensor ordering, Hamiltonian signs,
pair convention, and convention version are stored in the same typed result.
Channel gamma is in s\(^{-1}\) T\(^{-1}\), damping is dimensionless, and a
basis-to-channel value of `-1` excludes a nonmagnetic basis site without
creating a dynamical channel.

Fitting is a separate quadratic least-squares extrapolation in \(\eta\), with
scaled independent variables. The result records the fit range, sample count,
coefficients, component-wise RMS residuals, the closest-sample-to-intercept
delta, and exchange reciprocity error. With fewer than three fit samples, an
explicit \(\eta=0\) sample is required and used directly.

## Independent validation

`validate_coarse_material_small_q` evaluates the pair Hamiltonian directly for
planar spin spirals and compares its energy-density change with the extracted
continuum tensors. It stores atomistic and continuum energies, absolute
errors, the maximum tested wave vector, and a tolerance-scaled error. Static
phase offsets permit controlled acoustic and optical energy-sector fixtures,
but are not reported as dynamical mode extraction.

`tests/coarse_graining/test_stiffness_material.f90` covers:

- the analytic simple-cubic nearest-neighbour ferromagnetic stiffness;
- direction-dependent exchange and its diagonal anisotropic tensor;
- Cartesian cross terms in a skew cell;
- DMI sign and handedness for `D(+x)=+z`;
- direct positive- and negative-\(q\) atomistic spiral energies;
- regularization parameters, fit coefficients, residuals, and raw samples;
- a two-basis/two-channel chain with separate `A12`, `A21`, local exchange,
  acoustic phase-sector energies, and optical phase-sector energies;
- explicit rejection of ferri/AFM runtime consumption.

## Runtime gate and production consumption

A one-channel ferromagnetic descriptor passes the extraction gate only after
direct small-\(q\) energy validation and validation of positive channel gamma
and nonnegative damping.

Every ferrimagnetic/antiferromagnetic runtime request is rejected unless the
descriptor records both acoustic and optical **mode extraction** as validated.
CG-03 intentionally provides no operation that sets those mode-extraction
flags. Static two-basis spiral-energy agreement is useful evidence, but it is
not a substitute for the dynamical-matrix/eigenvector work assigned to the
later multi-channel task.

The one-channel descriptor is now consumed by
`AdaptiveCGProduction`: exchange stiffness and DMI spiralization are
constructed during normal `do_adaptive_cg=Y` setup and passed into the
coarse tensor operator. Production requires an explicit positive SI `alat`
because the physical displacement and cell volume in these sums determine
the J m\(^{-1}\) and J m\(^{-2}\) scales.

The following remain deferred:

- random-alloy material averaging for runtime use;
- dynamical acoustic/optical branch extraction and normalization;
- ferrimagnetic and antiferromagnetic runtime enablement;
- the separately recorded atomistic CPU/GPU DMI conformance fix.

## Acceptance

**Human sign-off:** The CG-03 signs, prefactors, units, conventions, and
runtime gate are approved.
