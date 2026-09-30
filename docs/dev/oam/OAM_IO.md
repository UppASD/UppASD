# OAM input and output

UppASD has two OAM observables. Trajectory OAM supports full (non-dilute)
lattices only; dilute systems are refused. `do_oam_traj Y` measures the real-space,
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
- `oam_gradient auto|fem|spectral`: gradient method, default `auto`. `auto`
  resolves to `spectral` for periodic in-plane cells on `USE_FFTW` builds when
  the first-Brillouin-zone image search is inside its bounds; otherwise it
  resolves to `fem`. MKL-FFT builds therefore resolve `auto` to `fem`. An
  explicit `spectral` requires `BC1 = BC2 = P` and `USE_FFTW`, and refuses
  when those conditions or the finite image search fail. The spectral method
  uses the first-Brillouin-zone fold defined in C17.
  Spectral trajectory OAM is valid only when the field content lies strictly
  inside each sublattice's first Brillouin zone. K- and M-centred packets are
  out of scope for both spectral and FEM gradients.
- `oam_axis x y z`: optional frame axis, normalised at `oam_init`; by default
  the axis is the normalised average moment at `oam_init`. A zero axis or an
  initial state with minimum alignment below `0.9` is refused. For boosted,
  driven or restart-loaded states, set this to the ground-state axis so the
  frame does not follow the packet.

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
FEM λ is biased at short wavelength; `auto` selects spectral where possible.
On a square lattice the measured FEM bias is about `−8.15%` at `k·a = 0.5`
and `−30.60%` at `1.0`, relative to `k·a = 0`.

When explicit `spectral` mode refuses because the image search reaches its
finite boundary, the in-plane cell vectors are not reduced. Re-express `C1,
C2` as a reduced basis; for example replace `C2 = (1.5, 0.3)` by
`C2 − C1 = (0.5, 0.3)`.

With `Mensemble > 1`, λ columns use norm-weighted first moments over valid
ensembles. `N_m` and `dSz_hbar` are arithmetic means; `Lz_tot_hbar` and
`balance` are then recomputed from the aggregate, and centroid diagnostics
use norm weights over centroid-valid ensembles. `R_x`, `R_y` are combined in
reduced coordinates over the centroid-valid ensembles, weighted by `n_k`.
Along a periodic axis the combination is a circular mean:
`s = atan2(Σ_k n_k sin 2π s_k, Σ_k n_k cos 2π s_k)/2π mod 1`. Along an open
axis it is the linear mean. The result is converted back to Cartesian.
`sigma_psi` stays the `n_k`-weighted mean. A warning reports excluded
ensembles once. For `Mensemble = 1`, output is unchanged.

The observable is intended for coherent packets. In a finite-temperature run,
the origin columns can remain finite while the centroid spread guard rejects
most or all `lambda_L_centroid` values for delocalised thermal magnons; this
is expected, not a crash. An exactly saturated initial FM sample has zero
transverse norm and is reported as `NaN` by the same guard before thermal
fluctuations develop.

Release timing measurements use `OMP_NUM_THREADS=1`, FFTW, and a frozen
localized packet on an Apple MacBook Pro (M1, 8 CPU cores, macOS 26.7.1).
No-OAM integration cost is the difference between matched 101-step and
one-step runs divided by 100. OAM cost per sample is the corresponding
101-step versus one-step difference with OAM enabled, minus that no-OAM
baseline, divided by 100; this cancels startup, mesh, FFTW-plan, and first
sample costs. These are reference timings, not an acceptance gate.

| grid | NA | no-OAM integration step | FEM OAM/sample | spectral OAM/sample |
|---:|---:|---:|---:|---:|
| 256² | 1 | 6.15 ms | 3.28 ms | 7.30 ms |
| 256² | 2 | 11.07 ms | 14.65 ms | 21.41 ms |
| 512² | 1 | 22.55 ms | 21.48 ms | 35.16 ms |
| 512² | 2 | 50.08 ms | 58.57 ms | 95.78 ms |

Use `oam_step` ≥ 10 for large systems; at `oam_step 1` spectral OAM can cost
several integration steps per step.

For boosted, driven or restart-loaded states, use `oam_axis` for the
ground-state axis. In the C13 check the default axis tilts by about `0.2°`
and shifts `lambda_L_centroid` by about `6%` (`0.903` versus `0.957`).
