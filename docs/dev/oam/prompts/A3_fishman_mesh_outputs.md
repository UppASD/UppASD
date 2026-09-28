# A3 + A4 — Fishman polar mesh, outputs, tests  ·  Sonnet  ·  commit `[OAM-A3]`

Implement Blueprint A sections **A3** and **A4**. Paste the shared preamble first. Requires `[OAM-A1]`.

## Tasks

- **A3.1.** Build polar q as Cartesian `k/(2π)`. Delete the reduced-coordinate helpers. Fix the comments that claim UppASD takes reduced q.
- **A3.2.** Correct the sign comment (contract C9).
- **A3.3.** Rename outputs and the key (C15). Make the pointwise `oam_k` output opt-in via `oam_lswt_pointwise`. Put the C13 sign label in the header.
- **A3.4** (optional). Implement the cost savings only if they are straightforward. Report measured memory and time before and after.
- **A4.** Put the harness into `tests/SpinWaves/oam_lswt/` and wire it into CTest. Update `test_fishman_oam.py` as described.

## Trap to avoid

The test `test_cartesian_reduced_round_trip_for_oblique_cell` passes whether or not the bug exists: it tests the conversion, not the consumer. Don't use it as evidence. The evidence is harness check [2]: the same lattice described with two different cells must give bit-identical F.

## Acceptance

- `run_lswt_checks.py`: ALL PASS. Paste the four lines.
- `ctest -L oam-lswt` passes.
- `f_oam`/`oam_k` no longer appear in `source/` except in the deprecated-alias branch.

**Commit `[OAM-A3]`, report, stop.** Gate GA follows: the maintainer decides the sign question (A5).
