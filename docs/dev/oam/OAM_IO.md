# OAM input and output

UppASD has two OAM observables. `do_oam_traj Y` measures the real-space,
finite-amplitude trajectory field and writes `oam_traj.<simid>.out`. It reports
intrinsic and envelope (extrinsic) OAM, so it is not the same quantity as the
LSWT result. `do_oam_lswt Y` measures band- and wave-vector-resolved intrinsic
OAM in the harmonic limit and writes `oam_lswt.<simid>.out`; it requires
`do_chern Y`. The LSWT diagnostics file is
`oam_lswt_diagnostics.<simid>.out`. The gauge-dependent pointwise output is
appended only when `oam_lswt_pointwise Y` is requested.

Trajectory inputs are:

- `do_oam_traj`: enable trajectory OAM (`Y`/`N`). The old `do_oam` key is a
  deprecated alias and enables the same path while printing a deprecation
  warning.
- `oam_step`: sampling interval, default `100`.
- `oam_buff`: number of rows buffered before writing, default `10`.
- `oam_origin x y z`: fixed origin for `lambda_L_origin`; without it, the
  arithmetic site centroid at initialization is used.
- `oam_weight site|area`: use unit site weights or FEM site areas; default
  `site`.
- `oam_sigma_max`: centroid-spread guard in units of half the shorter cell
  side; default `0.6`.
- `oam_gfactor`: positive g factor for `N_m`; zero uses `Landeg_glob`.
- `oam_sublattice i j ...`: optional one-based unit-cell sublattice list. If
  omitted, all sites are included. For `NA > 1`, omission also appends one
  `lambda_L_origin_sN lambda_L_centroid_sN N_m_sN` block per sublattice.

The first ten trajectory columns are:

```text
step lambda_L_origin lambda_L_centroid N_m Lz_tot_hbar dSz_hbar balance R_x R_y sigma_psi
```

Additional sublattice blocks are appended after these columns. Non-finite
values are written as `NaN`. The cumulant JSON key
`orbital_angular_momentum` is the running mean of finite trajectory
`lambda_L_centroid` samples, or `null` when there are none.

LSWT inputs are `do_oam_lswt`, `oam_nphi` (angular points per ring), `oam_nr`
(radial points), `oam_kmax` (Cartesian radius, with the inscribed-BZ default),
and `oam_lswt_pointwise`. The old names `do_magnon_oam`, `f_oam_nphi`,
`f_oam_nr`, and `f_oam_kmax` remain deprecated aliases for their corresponding
new keys.

The trajectory observable and LSWT observable are deliberately different
(C14): trajectory OAM contains intrinsic and envelope winding at arbitrary
amplitude, while LSWT OAM is intrinsic band OAM in the harmonic limit. In the
particle-only narrow-wavepacket bridge convention used by B5.4,
`lambda_L_centroid = l_envelope - 2 F_n(k0)/hbar`; this is a derived bridge
relation, not an assertion that the two standalone observables are equal.
