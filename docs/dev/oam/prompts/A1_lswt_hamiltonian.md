# A1 + A2 — LSWT Hamiltonian and diagonaliser  ·  Sonnet  ·  commit `[OAM-A1]`

Implement Blueprint A sections **A1** and **A2** exactly as written. Paste the shared preamble first.

## Before coding

1. Build `b60b762` and run the harness to record the RED baseline:
   ```
   R=$(git rev-parse --show-toplevel); mkdir -p /tmp/oam-lswt && cd /tmp/oam-lswt
   UPPASD=$R/build-oam/bin/uppasd python3 $R/tests/SpinWaves/oam_lswt/run_lswt_checks.py
   ```
   Always run the harnesses from a scratch directory: they create run directories in the current directory.
   Expect [1] Goldstone, [2] invariance, [3] D=0 run and [4] oracle all to FAIL. Keep the output for the report.
2. Read `reference/lswt_reference_fixes.diff`. It shows the DM loop and the `diamag_eps` default. It doesn't cover SA, PD, AMS, A2.2 or A2.3.

## Tasks

- **A1.** Give DM, SA and PD their own neighbour-list loops in `setup_Jtens_q`, and fix the AMS distance index. Report what you found in `setup_Jtens2_q`.
- **A2.1.** Restore the `diamag_eps` default, and list the call paths that reach `setup_diamag`.
- **A2.2.** Add the instability warning / error stop, and the key `nc_allow_unstable`.
- **A2.3.** Run the fallback on the same K as the main path, plus one test that forces it.
- **A2.4.** Add the one-line note when `hfield` is set; add a docs sentence.

## Physics you must not change

- The sign conventions `J_n = -D_n`, `+sa2tens(-sa_vect)`, `+pd2tens(-pd_vect)`. They are right; only the bond indexing is wrong.
- The `2*dia_eps` diagonal shift. Restoring its default is the fix; don't redesign it.

## Acceptance

- Harness checks [1] and [3] GREEN. [2] and [4] stay RED until A3; say so explicitly.
- Kagome example: `bphase*.out` bitwise identical before and after, with the shortened settings in Blueprint A1.
- The new SA test passes, or is escalated if its expected physics is unclear.

**Commit `[OAM-A1]`, report, and continue to A3** — there is no gate here.
