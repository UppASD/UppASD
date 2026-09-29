# B5 pyswatter cases

The trajectory harness creates these case names under `traj_checks/` when it
is run with a real UppASD binary:

| Case | Purpose | UppASD setting |
|---|---|---|
| `ell+1` | Counter-clockwise vortex sign pin (`m_x+i m_y` winding) | site-weighted, centred vortex |
| `shift_b` | Fixed-origin comparison after a rigid basis shift | `oam_origin 0 0 0` |
| `w_area` | Explicit FEM-area weighting comparison | `oam_weight area` |

Run the real harness from a scratch directory first so it writes the
`coord.oamtest.out` and `restart.oamtest.out` inputs expected by pyswatter:

```sh
mkdir b5-pyswatter-run
cd b5-pyswatter-run
UPPASD=/path/to/uppasd python3 /path/to/UppASD/tests/SpinWaves/oam_lswt/run_traj_checks.py
```

From that scratch directory, invoke
`/path/to/UppASD/docs/dev/oam/pyswatter/run_pyswatter_checks.sh`. Its `site`
integration flag is the exact B5.2 command from the blueprint; the `w_area`
case records the explicit UppASD area-weighting variant for the maintainer to
compare and report.

For a like-for-like `shift_b` check, rotate the moment file into UppASD's C8
global frame and compare the shared origin-referenced observables:

```sh
python3 docs/dev/oam/pyswatter/compare_shift_b_c8.py traj_checks/shift_b
```

This writes `restart.oamtest.c8.out` and `ref_c8.csv` in the case directory.
It compares `lambda_L`, `N_m`, and `delta_Sz_over_hbar`; total `Lz` and balance
are intentionally not compared because pyswatter uses the origin-referenced
lambda while the Fortran trajectory total uses the centroid-referenced value.

## B5.2 disposition

`ell+1` and `w_area` pass the `lambda_L` comparison within `1e-3`. `shift_b`
is **N/A for pyswatter cross-validation**: it requires periodic minimum-image
triangulation for the rigidly shifted periodic mesh, which the current
pyswatter CLI does not implement. Its C8-rotated `N_m` agrees with UppASD and
the independent oracle agrees with UppASD's `lambda_L`; the remaining
pyswatter mismatch is therefore a missing oracle feature, not an UppASD
failure. Do not use `shift_b` to adjust the Fortran implementation.

Do not check in `traj_checks/`, `coord.oamtest.out`, `restart.oamtest.out`,
`restart.oamtest.c8.out`, `ref.csv`, or `ref_c8.csv`; they are run outputs.
