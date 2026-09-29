# R1 — Recover or redo A1/A2 · Sonnet · `[OAM-A1]`

A1/A2 never reached the branch. `diamag.f90` and `ams.f90` are unchanged from `b60b762`, and `!diamag_eps=-1.0_dblprec` is still commented out at `diamag.f90:110`. The work existed: the A5 numbers need the DM fix. So:

1. **Search first.** Run `git worktree list`, `git branch -a`, `git stash list` and `git reflog | grep -i "A1\|diamag"`. If a commit or stash with the A1/A2 changes exists, rebase it onto the branch tip and resolve conflicts.
2. **Otherwise redo it** from Blueprint A §A1–A2, using `docs/dev/oam/reference/lswt_reference_fixes.diff` for the DM loop and the default. It doesn't cover SA, PD, AMS, A2.2 or A2.3.

**Acceptance:** from a scratch directory, `python3 <repo>/tests/SpinWaves/oam_lswt/run_lswt_checks.py --binary <build>/bin/uppasd` gives 4/4 PASS. The branch currently gives 1/4. Report which route (recovered or redone) you took.
