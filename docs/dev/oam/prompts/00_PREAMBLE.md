# Shared preamble — paste above every prompt in this folder

You are working on the UppASD repository (Fortran, CMake) on a branch created from `v6.1.0rc` at `b60b762`. The planning documents are in `docs/dev/oam/`:
- `CONVENTIONS_OAM.md` — the binding physics contract;
- `BLUEPRINT_A_LSWT.md` and `BLUEPRINT_B_TRAJECTORY.md` — the work packages;
- `tests/SpinWaves/oam_lswt/` — independent reference implementations and acceptance harnesses;
- `OAM_QUESTIONS.md` — the escalation log. Create it if missing.

Read the contract and your work package in full before touching code.

## Standing rules

1. **Branch and commits.** Work on `feat/oam-refactor`; create it from `v6.1.0rc` if it does not exist. End every work package with **one commit** whose message starts `[OAM-<WP>]`. Never commit to `v6.1.0rc`, never merge, never push, never tag. The maintainer reviews and pushes. Before committing:
   - `git status --short --untracked-files=all` shows only intended files, and every new file is added;
   - the full tree configures and builds from scratch: `cmake -S . -B build-oam && cmake --build build-oam -j`.
2. **The oracles are not yours to fix.** If the Fortran disagrees with anything in `tests/SpinWaves/oam_lswt/`, the Fortran is wrong until proven otherwise. If you believe an oracle is wrong, write the argument into `OAM_QUESTIONS.md` and stop. Never change a tolerance.
3. **Escalate, don't resolve.** Any ambiguity in the contract, any physics choice not written down, and any test that can only pass by changing expected values: log it in `OAM_QUESTIONS.md` (date, WP, question, options, your recommendation) and stop that item.
4. **Scope.** Change only the files your work package names, plus build files and tests. If you find an unrelated bug, note it in `OAM_QUESTIONS.md`; don't fix it.
5. **Code style.** Follow the surrounding UppASD style: `dblprec`, `memocc` on every allocate/deallocate, `!>` Doxygen headers, OpenMP clauses spelled out. Comments describe what the code does and why. They never narrate the refactor ("moved in WP3", "see question log").
6. **Evidence, not assertion.** Your final report contains:
   - the commit hash and `git show --stat`;
   - the exact commands you ran and the tail of their output (harness GREEN/RED lines);
   - every entry you added to `OAM_QUESTIONS.md`;
   - anything you did not finish.
7. **Gates.** If your work package ends in a gate (G0–G3, GA), stop after the report. Do not start the next work package.
