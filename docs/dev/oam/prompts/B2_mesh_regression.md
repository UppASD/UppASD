# B2 — Mesh-fix regression  ·  Sonnet  ·  gate G1  ·  commit `[OAM-B2]` (data and log only)

Blueprint B section **B2**. Paste the shared preamble first. **No source changes.**

1. Build `b60b762` and `[OAM-B1]` side by side.
2. Run skyrmion number (`skyno T`) and chirality (`do_chiral Y`) on `tests/kagome`, `tests/HeisStripe`, `tests/triang2D` and one skyrmion example (pick one under `examples/` and name it). Use identical seeds and a short Nstep.
3. Add a table to `OAM_QUESTIONS.md`: case, observable, before, after, relative change.

The skyrmion number must be unchanged. If it moves, write "SECOND DEFECT SUSPECTED", with the case, and stop. Chirality may change by roughly the fraction of wrap triangles; compute that fraction and put it next to each chirality delta.

**Commit the log and inputs, report, stop for maintainer sign-off.**
