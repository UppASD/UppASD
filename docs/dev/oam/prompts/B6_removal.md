# B6 — Remove the pseudo-OAM  ·  Sonnet  ·  commit `[OAM-B6]`

Blueprint B section **B6**. Paste the shared preamble first. Requires G3 sign-off.

Delete the listed routines, commented-out blocks and module variables from `topology.f90`. Keep `minimal_image_correction` only if `git grep` shows a remaining caller.

Comment rule: at most one line in `topology.f90`, e.g. `! Trajectory OAM: see oam.f90`. No history, no work-package references.

## Acceptance

- Clean build with `-Wall`: no new unused-variable warnings from `topology.f90`, `oam.f90` or `mesh2d.f90`.
- `run_traj_checks.py` and `run_lswt_checks.py`: ALL PASS.
- `git grep -n "calculate_oam\|oam_tri\|Lz_csum\|S0_arr"` returns nothing under `source/`.
- `git log --oneline v6.1.0rc..HEAD` shows the full series `[OAM-A1]` … `[OAM-B6]`.

**Commit `[OAM-B6]`, report.** The maintainer squashes or merges.
