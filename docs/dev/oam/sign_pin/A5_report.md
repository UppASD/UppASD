# A5 sign pin: Fishman FM honeycomb

Date: 2026-09-28

This is the GA report for A5. The user instruction not to commit artifacts or
blobs is followed: the A5 commit contains only the text inputs and this report;
no source changes, binaries, run directories, or output files were added.

## Fishman model and conventions

The reference is Fishman, Berlijn, Villanova and Lindsay, *Phys. Rev. B* 107,
214434 (2023), arXiv:2304.07379:

<https://arxiv.org/abs/2304.07379>

The equations used for the comparison are:

- Eq. (8): `H_2 = sum_k' v_k^dagger L(k) v_k`.
- Eqs. (17)-(18): the magnon OAM uses the paraunitary eigenvector
  `X^(-1)` and `O_n(k) = -(i hbar/2) [k x <u_n|d/dk|u_n>] . z`.
- Eqs. (25) and (27):
  `F_n(k) = (1/(2 pi)) integral dphi O_n(k,phi)` and
  `O_n,av(k) = (2/k^2) integral_0^k dq q F_n(q)`.
- Eq. (28): the FM honeycomb quadratic Hamiltonian has prefactor
  `3 J S / 2`, with `G_k = d Theta_k`; Eqs. (29)-(30) define `Theta_k`
  and the NN structure factor `Gamma_k`.
- Figure 1(c) and the accompanying text define
  `d = -2 D / (3 J)`, with `J > 0` and `D < 0` for `d > 0`.
- Eqs. (31)-(32) give `hbar omega_n = 3 J S (1 + kappa +/- eta_k)`;
  the paper states that the easy-axis anisotropy shifts the energies but does
  not affect the OAM.
- The paper uses `a` as the NN distance in the honeycomb wave-vector
  convention and gives `k_max = (2 sqrt(3)/9)(2 pi/a)`, approximately
  `0.385 (2 pi/a)`. In the UppASD coordinate choice here the NN distance is
  `a = 1/sqrt(3)`, so `|K| = 4 pi/3 = 4.1887902047863905`.

### Mapping into UppASD

The model uses `J=1` and `S=1` in the input units, with the NNN DMI attached
to the three oriented A-A and B-B bonds and their reverse bonds. The paper's
published `d=0.1` therefore requires

```text
D_paper/J_paper = -(3/2) d = -0.15.
```

The factors can also be tracked in the Fourier reduction. The three NN
vectors give `sum_delta exp(i k.delta) = 3 Gamma_k`, producing the factor 3
in the honeycomb normalisation. Pairing each oriented NNN vector with its
reverse produces the sine term with the factor 2; matching that term to the
`3JS/2` matrix in Eq. (28) gives the displayed `-2D/(3J)` relation (the
remaining sign is the choice of NNN arrow orientation). The NN and NNN
quadratic terms carry `JS` and `DS`; hence `S` cancels from `d`, and changing
`S` while keeping `D/J` fixed does not change the eigenvectors or the OAM. The input `kfile` contains
`K_1=-0.05`, the UppASD easy-axis form corresponding to the paper's
`-K sum_i S_iz^2`; it only removes the Gamma zero mode.

The two inputs deliberately test both signs of the UppASD `dmfile` coefficient.
This also exposes the orientation issue that GA must settle: for the fixed
oriented bonds in `dmfile`, the positive UppASD sign produces positive band-1
OAM, while the algebraic paper value `D_paper/J=-0.15` has the opposite sign
under the same written bond orientation. The sign result is therefore reported
without silently relabelling the bond orientation.

## Runs

Binary: `build-oam/bin/uppasd` (current A3 worktree).

Exact command, run once in each case directory:

```text
env OMP_NUM_THREADS=1 /Users/andersb/Jobb/UppASD_6.1/build-oam/bin/uppasd > run.log 2>&1
```

The input requests `do_oam_lswt Y`, `f_oam_kmax 4.1887902047863905`,
`oam_nphi 128`, and `oam_nr 32`. The `f_oam_kmax` deprecation warning and
the warning that the rings extend beyond the inscribed radius are expected.

Relevant run-log tail lines for both signs:

```text
Maximum paraunitarity error:  7.86138E-15
Warning: Fishman OAM rings extend beyond the inscribed BZ radius:  4.18879E+00  3.62760E+00
Maximum paraunitarity error:  7.86138E-15
Fishman OAM polar mesh (Nphi,Nr,kmax,k_in):    128     32  4.18879E+00  3.62760E+00
Fishman OAM written to oam_lswt.honey.out
Chern calculation done.
```

The band-1 output row at the K endpoint (`k=4.18879020478639`) is:

| UppASD `D/J` | `F_1(K)/hbar` | `O_1,av(K)/hbar` | peak location | comparison with +0.236 |
|---:|---:|---:|---|---:|
| `+0.15` | `+0.758610060486273` | `+0.190903724924406` | `K` | `-19.1%` |
| `-0.15` | `-0.758610060486273` | `-0.190903724924406` | `K` | opposite sign, same magnitude |

The extrema occur at the K endpoint for both signs. Neither sign reproduces
the published magnitude `0.236 hbar` to about 5%; the magnitude discrepancy is
approximately 19.1%. The positive UppASD sign reproduces the published sign,
but not its magnitude. GA is therefore not a complete absolute-sign pin until
the maintainer decides how the UppASD oriented DM-file convention maps onto the
paper's drawn NNN orientation and resolves the independent magnitude mismatch.

The intended work-package commit is `[OAM-A5]` and contains only this report
and the ASCII input files under `sign_pin/`. No generated artifact or blob is
part of it.

## A5b rerun — corrected D mapping (2026-09-28)

R4 first reran the original `D/J = +/-0.15` inputs to check the K-point
normalisation from UppASD output. At `K = 4.188790204786391`, both signs gave
`E_1 = 205.6867800742812 meV` and `E_2 = 120.8500727235184 meV`, hence
`E_1 - E_2 = 84.8367073507628 meV = 1.5588 (4*ry_ev)`. This is
`6 sqrt(3) D S` for `D = 0.15` and `S = 1`; the paper's relation is twice
that value, so `D_UppASD = 2 D_paper` and the published `d = 0.1` requires
`|D_UppASD|/J = 0.30`.

Every D coefficient in `plus_D/dmfile` and `minus_D/dmfile` was then doubled
from `+/-0.15` to `+/-0.30`. The corrected runs used isolated copies of the
two input directories and:

```text
env OMP_NUM_THREADS=1 /Users/andersb/Jobb/UppASD_6.1/build-oam/bin/uppasd > run.log 2>&1
```

Both runs reported a maximum paraunitarity error of `7.86138E-15`, the
expected warning that `kmax` exceeds the inscribed BZ radius, and wrote the
Fishman OAM output. The corrected K-point gap was
`169.6734211892256 meV`, twice the baseline gap.

The band-1 K-point rows are:

| UppASD `D/J` | `F_1(K)/hbar` | `O_1,av(K)/hbar` | peak location | comparison with +0.236 |
|---:|---:|---:|---|---:|
| `+0.30` | `+0.761021860661387` | `+0.227959688312056` | `K` | `-3.4069%` |
| `-0.30` | `-0.761021860661387` | `-0.227959688312056` | `K` | opposite sign, same magnitude |

The positive result agrees with the independent oracle's approximately
`0.228` and is about 3.4% below the published `0.236 hbar`. For the sign
oscillation described in the paper, the corrected `F_1(k)` samples remain
nonnegative for `+D` and nonpositive for `-D` throughout `0 <= k <= K`, apart
from the zero at `k=0`; no sign oscillation is observed. This observation is
recorded without resolving the discrepancy. GA remains with the maintainer.
