# B4 — Input, aliases, sublattices, averages, docs  ·  Sonnet  ·  commit `[OAM-B4]`

Implement Blueprint B section **B4**. Paste the shared preamble first. Requires `[OAM-B3]`.

## Tasks

- Add the keys `oam_gfactor` and `oam_sublattice`, and the per-sublattice output block (C6).
- Make `do_oam` a deprecated alias that turns on `do_oam_traj` and prints a line containing "deprecated". **Do not** leave it as an orphaned variable fixed at 'N' — that was the mistake in the earlier attempt.
- Make `mesh2d_build` a single call site in `uppasd.f90` for all mesh consumers.
- Buffering: no rows lost when the run ends mid-buffer.
- `prn_averages.f90`: `"orbital_angular_momentum"` is the running mean of the finite `lambda_L_centroid` values, or `null`.
- Docs: a short section covering both OAM paths, every key, the output columns, and the C14 warning that they are different observables.

## Acceptance

- `run_traj_checks.py`: ALL PASS, including C10.
- A two-sublattice case (e.g. the honeycomb from `tests/SpinWaves/oam_lswt/mkhoney.py` with a prescribed texture) writes per-sublattice columns, and their `N_m` sum to the total.
- The `tests/` regression suite is unchanged.

**Commit `[OAM-B4]`, report.**
