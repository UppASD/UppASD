# OAM conventions contract

Binding for every change made under Blueprints A and B. Each clause is invisible to a "does it run" test, which is why it is written down before any code. Clauses marked **OPEN** need maintainer sign-off (gate G0) before the work that depends on them starts.

If the code or a physics question is not settled by this file, stop and write the question into `docs/dev/oam/OAM_QUESTIONS.md`. Do not choose a branch and continue.

---

## Part I — Trajectory OAM (`do_oam_traj`, pyswatter formulation)

**C1 — Wavefunction and sign.** `psi_i = m_x,i + i*m_y,i`, where m is the unit moment direction expressed in the frame fixed at `oam_init` (C8). Under Holstein–Primakoff about +z, `S+ = S_x + i*S_y = sqrt(2S)*a`, so `psi` is proportional to the magnon annihilation field and `lambda_L > 0` means magnon OAM along +z. This matches pyswatter. It is the complex conjugate of the `m_x - i*m_y` convention used in parts of the literature; say so in the module header so nobody "corrects" it.

**C2 — No deviation subtraction.** `psi` is built from the transverse components directly: `|psi_i|^2 = 1 - m_z,i^2`. It is never built from `m - S0`. A reference direction only defines the frame. Evidence for why: the current `oam_tri_improved` measures deviation from the configuration at setup, so a prescribed vortex texture loaded by restart file returns exactly 0.

**C3 — Operator and dual reference.** With linear-FEM gradients (C7):

```
ell_z(i)  = Im[ conj(psi_i) * ( (x_i - x0) * dpsi_dy(i) - (y_i - y0) * dpsi_dx(i) ) ]
lambda_L  = sum_i ell_z(i) w_i / sum_i |psi_i|^2 w_i
```

Report two values every sample:

- `lambda_L_origin`: `x0` = `oam_origin` (input key, 3 reals). Default: arithmetic centroid of the sampled sites, fixed at `oam_init`. Lever arm = raw coordinates minus `x0` (no minimum image).
- `lambda_L_centroid`: `x0 = R`, the instantaneous `|psi|^2 w` centroid. Along periodic directions `R` is the circular mean in reduced coordinates, mapped back into the cell, and lever arms are taken by minimum image to `R`. Along open directions it is the arithmetic mean.

Referencing to `R` removes the drift term `(R x P)_z`, **not** the extrinsic part: `lambda_L_centroid` still contains the envelope winding `l`. Say exactly this in the output header. A stationary packet (P = 0) has `lambda_L_origin = lambda_L_centroid` for any origin; only a moving packet separates them.

**C4 — Weighting. OPEN.** Two options: `oam_weight = site` (`w_i = 1`) or `area` (`w_i = A_i = (1/3) sum_{D in i} A_D`). Proposed default: `site`, matching the site-summed magnon count `N_m` and pyswatter's `spin-oam-balance`. pyswatter's standalone `lz` defaults to `area`, so the choice must be recorded in the output header. Maintainer confirms the default at G0.

**C5 — Guards.**
- Norm `sum |psi|^2 w` below `1e-14`: both lambda columns `NaN`, one warning per run, no division.
- Spread `sigma_psi` (the `|psi|^2 w` RMS radius about `R`) above `oam_sigma_max` × half the shorter cell side (default `0.6`): `lambda_L_centroid = NaN`, `lambda_L_origin` stays valid. For reference, a field spread uniformly over the cell has `sigma/half ≈ 0.82`, and a localised l=2 vortex of width 5a in a 41a cell has `sigma/half ≈ 0.4`.
- Sites with no gradient support (`site_wsum == 0`) are dropped from both numerator and denominator.

**C6 — Output contract.** File `oam_traj.<simid>.out`. The header is written at the first *sample*, not at `mstep <= 1`, so restarts get a header. Header lines start with `#` and state: the psi convention, the weighting, `g`, the origin, and the C3 sentence about the centroid. Then one row per sample with exactly these columns:

```
step  lambda_L_origin  lambda_L_centroid  N_m  Lz_tot_hbar  dSz_hbar  balance  R_x  R_y  sigma_psi
```

- `N_m = sum_i (mmom_i/g) (1 - m_z,i)` (site sum)
- `dSz_hbar = N_m`
- `Lz_tot_hbar = N_m * lambda_L_centroid`
- `balance = dSz_hbar + Lz_tot_hbar`
- Non-finite values are written as `NaN`.
- When more than one sublattice is present and `oam_sublattice` is not set, append one block of `lambda_L_origin lambda_L_centroid N_m` per sublattice.

Sampling fires when `mod(mstep - rstep - 1, oam_step) == 0`, i.e. relative to the start of the measurement phase. `oam_buff` rows are buffered before writing.

**C7 — Mesh.** One module (`mesh2d`) owns the triangulation for every consumer: OAM, skyrmion number and chirality.
- Vertices are placed by minimum image relative to vertex 1 along periodic directions; wrap cells are dropped along open directions.
- Triangles are oriented counter-clockwise, so the signed area is positive.
- `SystemData::coord` is never modified.
- At setup, print exactly one line, which `run_traj_checks.py` parses:

```
Mesh2D: ntri= <n> total_area= <x> cell_area= <y> degenerate= <k>
```

A periodic mesh must tile the cell exactly (`total_area = N1*N2*|C1 x C2|`). An open mesh covers `(N1-1)*(N2-1)*|C1 x C2|`.

**C8 — Frame.** Phase 1: a single global frame, with `e_z` along the normalised average moment at `oam_init` and `e_x, e_y` completing a right-handed set. If `min_i m_i . <m> < 0.9`, refuse with a clear message and produce no numbers. Smooth local frames for non-collinear states are out of scope, and must not be attempted with the existing `local_frame`: its `if (abs(ez(1))<0.9)` seed choice is discontinuous and injects phase jumps into the gradient.

## Part II — LSWT magnon OAM (`do_oam_lswt`, Fishman formulation)

**C9 — Bogoliubov matrix and formula.** `T` is the paraunitary Colpa matrix (`a = T b`), which is Fishman's `X^{-1}`. Its positive-energy columns are `1:NA` and it satisfies `T^dagger eta T = eta`. The pointwise OAM is

```
O_n(k)/hbar = -1/2 Im[ T_n^dagger eta dT_n/dphi ]
```

This is algebraically Fishman PRL 129, 167202, Eq. 11, where the operator acts on `X^{-1*}`; converting to `T^dagger eta dT` gives the −½. The sign does **not** come from UppASD's `exp(-i k.R)` convention: in 2D, `k -> -k` is a rotation by π and leaves ring averages unchanged. Fix the comment that says otherwise.

**C10 — Gauge-invariant outputs.**
- `F_n(k)` is the average of `O_n` over a circle of Cartesian radius `k` about Γ, in a gauge single-valued on the disk.
- `O_n,av(k) = (2/k^2) ∫_0^k q F_n(q) dq`.
- Identity used as a test: `F_n(k)` equals the Berry phase of the ring divided by 4π, i.e. the Berry flux through the disk over 4π (Stokes).
- The pointwise `O_n(k)` depends on the gauge. It may only be written to a diagnostics file, and only on request.

**C11 — Reciprocal coordinates.** Everything that reaches `setup_Jtens_q` / `setup_ektij` is a Cartesian q in units of 2π/alat, i.e. `q = k/(2π)`. Reduced (reciprocal-basis) coordinates are never passed to the Hamiltonian.

**C12 — Neighbour lists.** Exchange, DM, symmetric-anisotropic (SA) and pseudo-dipolar (PD) couplings each live on their own list (`nlist`, `dmlist`, `salist`, `pdlist`), which in general differ in length and order. A coupling vector must always be indexed with its own list.

**C13 — Absolute sign of F. OPEN.** Pin it against Fishman, Berlijn, Villanova, Lindsay, PRB 107, 214434 (2023): for the FM honeycomb with NNN DM, O_1,av peaks at 0.236ħ for d = 0.1, where a is the NN distance and the rings extend to K. The maintainer maps d to UppASD's D/J and signs off. Until then the sign is labelled "convention: see OAM_QUESTIONS.md" in the output header.

## Part III — Both

**C14 — Two different observables.** The trajectory path measures intrinsic plus extrinsic OAM at arbitrary amplitude. The LSWT path measures band- and k-resolved intrinsic OAM in the harmonic limit, which is identically zero for a single-sublattice ferromagnet. They meet only in the narrow-wavepacket limit, `lambda_L_centroid ≈ l_envelope + F_n(k0)/hbar`. **Never write a test asserting they are equal.** The only comparison allowed is the bridge test (B5.4).

**C15 — Naming.**

| Old key | New key | Old file | New file |
|---|---|---|---|
| `do_oam` | `do_oam_traj` | `oam.<simid>.out` | `oam_traj.<simid>.out` |
| `do_magnon_oam` | `do_oam_lswt` | `f_oam.<simid>.out` | `oam_lswt.<simid>.out` |
| — | — | `f_oam_diagnostics.<simid>.out` | `oam_lswt_diagnostics.<simid>.out` |
| — | — | `oam_k.<simid>.out` | only with `oam_lswt_pointwise Y`, appended to diagnostics |

Old keys stay as aliases: they set the new variable and print one line containing the word "deprecated". Other files that used the legacy trajectory flag (for example the cumulant JSON in `prn_averages.f90`) read the new variable.

**C16 — Escalation.** Any disagreement between the implementation, an oracle and this contract goes to `OAM_QUESTIONS.md` and stops that work package. Never adjust an oracle to match the Fortran.
