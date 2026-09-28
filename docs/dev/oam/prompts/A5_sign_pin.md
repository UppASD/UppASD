# A5 — Absolute sign of F_n against Fishman (2023)  ·  Opus  ·  gate GA  ·  commit `[OAM-A5]` (inputs and report only)

Paste the shared preamble first. Requires `[OAM-A3]`. **No source changes.**

Goal: give the maintainer what they need to fix the sign label (contract C13).

1. Read Fishman, Berlijn, Villanova, Lindsay, PRB 107, 214434 (2023), arXiv 2304.07379. Extract the FM honeycomb Hamiltonian with NNN DM, the definition of `d`, and the lattice constant convention (a = NN distance). Quote the equations you rely on in the report.
2. Build the model in UppASD units (start from `tests/SpinWaves/oam_lswt/mkhoney.py`). State the mapping `d ↔ D/J` explicitly, including every factor of 2 and the S dependence, and justify each factor.
3. Run `do_oam_lswt` with `f_oam_kmax = |K|` (the warning is expected) for `D = +|D|` and `D = −|D|`, at the published `d = 0.1`.
4. Report `O_1,av(k)` for both signs next to the published peak of 0.236ħ, and whether the peak sits at K.

Do not edit the −½, the ring algorithm or any oracle. If neither sign reproduces the magnitude to about 5%, report that as the finding.

Add inputs under `docs/dev/oam/sign_pin/`. **Commit, report, stop.** The maintainer decides.
