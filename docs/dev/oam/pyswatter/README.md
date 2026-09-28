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

Do not check in `traj_checks/`, `coord.oamtest.out`, `restart.oamtest.out`, or
`ref.csv`; they are run outputs.
