# WP10.7a Sol report — finite uniform coarse projection

Date: 2026-07-28

## Result

The backend-neutral OPEN_FFT host builder now accepts uniform, exactly
divisible coarse blocks.  It constructs the exact finite projection of the
accepted block-one point-dipole operator.  No Ewald, surface, continuum-cell,
Newell, or point-macrospin self-demagnetizing term is present.

Production enablement remains WP10.7b.  The existing production gate still
rejects non-unit OPEN_FFT blocks.

## Derivation

Let `A=(g,a)` and `B=(h,b)` denote basis-resolved coarse channels, let
`s=(s1,s2,s3)` be the common block shape, and let
`n=s1*s2*s3` be the common population.  For block moment

```text
M_A = sum_(i in A) m_i,
```

the projected subspace assigns `m_i=M_A/n`.  If `R_iA=1/n` for `i in A` and
zero otherwise, restriction of the finite point Hamiltonian gives

```text
K_coarse = R^T K_point R

K_coarse(d,a,b) =
  1/n^2 sum_(u in block,v in block)
    T(sum_axis (d_axis*s_axis + u_axis-v_axis)*C_axis
      + Bas(a)-Bas(b)).
```

Here `d=g-h` and `T(r)=3 rr^T/|r|^5-I/|r|^3`.  The only omitted summand is
`d=0`, `a=b`, and `u=v`, which is the true point self `K_point(i,i)=0`.
Consequently, `K_coarse(0,a,a)` generally is not zero: it contains every
ordered distinct intra-block atom pair.  Its normalization is `1/n^2`
because both the source lift and target conjugate-field restriction contribute
`1/n`.

At `s=(1,1,1)`, `n=1`, `u=v=0`, so the expression reduces exactly to the
accepted basis-resolved block-one finite operator.

## Geometry contract

The builder now receives both the complete primitive-cell
`atomistic_grid=N` and coarse `active_grid=G`.  Before tensor allocation it
requires, for every axis,

```text
s_axis > 0
N_axis % s_axis == 0
G_axis == N_axis/s_axis.
```

Thus every implicit active basis channel has the identical complete Cartesian
offset set and population `n`; partial edges and inconsistent macro grids are
rejected deterministically.  Checked arithmetic covers atomistic, active,
FFT, block-population, pair, batch, element, and byte counts.

Diagnostics distinguish:

- `max_point_self_abs`, which remains exactly zero; and
- `max_projected_diagonal_abs`, which records the finite coarse diagonal.

For a single orthogonal two-point block along x, the tested projected
diagonal is exactly `diag(1,-1/2,-1/2)`.  This comes solely from the two
ordered distinct point pairs.

## Independent cross-checks

`test_host_open_dipole_kernel.cpp` independently forms `P^T K_point P` from
absolute fine-point positions for a representative target/source block at
every signed coarse displacement.  Expected tensors do not call or inspect
the host builder.

The current fp64 maxima are:

| Case | Shape and basis | Maximum tensor difference |
|---|---|---:|
| Non-cubic orthogonal | `N=4x2x1`, `s=2x1x1`, `NA=1` | `4.441e-16` |
| Skew | `N=2x4x2`, `s=1x2x2`, `NA=2` | `9.770e-15` |

For arbitrary uniform-per-block moments, independently evaluated fine fields
were averaged back to their conjugate macro fields.  The projected
convolution agreed to `4.441e-16`, and fine/coarse energies agreed to
`1.110e-16`.

The pre-existing complete block-one impulse matrix remains compared with
Luna's direct finite point sum; its maximum difference is `7.816e-14`.
An additional block-one projection check has zero projected same-basis
diagonal and agrees with the independent point-operator construction within
the same fp64 construction budget.

## Nonuniform-moment approximation and convergence

For moments that vary inside a block, the solver replaces them by the block
average `M_A/n`; it is therefore an approximation outside the projected
subspace.  A fixed skew `4x2x1`, `NA=1` fixture was evaluated with affine,
nonuniform vector moments.  Comparing the unrestricted fine operator with
coarse fields assigned to member atoms gave:

| x block width | Maximum field error | Energy error per atom |
|---:|---:|---:|
| 4 | `5.144e-1` | `4.058e-4` |
| 2 | `4.377e-1` | `2.833e-4` |
| 1 | `6.661e-16` | `1.388e-17` |

The fixture demonstrates decreasing error as the block width is refined and
exact recovery at width one.  The numerical values are fixture-dependent;
the general statement is that refinement enlarges the represented
piecewise-uniform subspace, with the unrestricted finite operator recovered
at block one.
