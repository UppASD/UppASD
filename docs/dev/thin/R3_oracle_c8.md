# R3 — Oracle C8 frame · Sonnet · `[OAM-R3]`

The maintainer authorises this oracle change. It resolves the B3 G2 escalation: the Fortran is right; the oracle skipped C8.

In `tests/SpinWaves/oam_lswt/oracle_traj.py`, add `to_c8_frame(m)` and call it at the top of `evaluate()`:
- `e_z = sum(m)/|sum(m)|`;
- seed x̂, or ŷ if `|e_z·x̂| >= 0.9`;
- `e_x` = seed orthogonalised against `e_z` and normalised;
- `e_y = e_z × e_x`;
- return `(e_x·m, e_y·m, e_z·m)`.

This matches `oam_init` in `oam.f90`. Update the docstring to say the harness texture defines the frame. Touch nothing else.

**Acceptance:** `--selftest` ALL PASS; against the branch binary, C3 now passes and every check is ALL PASS.

**Don't** add an `oam_axis` key. That is a contract change: log it in `OAM_QUESTIONS.md` as a proposal (explicit frame axis, default = average moment at init, for restart-loaded or driven states).
