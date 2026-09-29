# UppASD OAM fix package

Regenerated 27 Sep 2026 against `v6.1.0rc` at `b60b762`. It replaces the 18 Sep blueprint, none of which reached the remote.

Copy this folder into the repository as `docs/dev/oam/` and commit it on `feat/oam-refactor` before dispatching the first prompt. The prompts refer to that path.

## Contents

| File | Purpose |
|---|---|
| `CONVENTIONS_OAM.md` | Binding physics contract, C1–C16. C4, C13 and C14 are confirmed; C14 uses the particle-field bridge `l - 2F_n`. |
| `BLUEPRINT_A_LSWT.md` | Fixes upstream of the Fishman magnon OAM: DM/SA/PD neighbour lists, `diamag_eps`, polar mesh, outputs, sign pin. |
| `BLUEPRINT_B_TRAJECTORY.md` | Replaces the real-space pseudo-OAM with the pyswatter formulation: mesh, kernel, I/O, validation, removal. |
| `prompts/00_PREAMBLE.md` | Standing rules; paste above every prompt. |
| `prompts/*.md` | One paste-ready prompt per work package. |
| `tests/SpinWaves/oam_lswt/oracle_traj.py` | Independent trajectory-OAM reference. `python3 oracle_traj.py` runs its self-test. |
| `tests/SpinWaves/oam_lswt/mkfixture_traj.py` | Writes UppASD run directories with a prescribed texture, loaded by restart file. |
| `tests/SpinWaves/oam_lswt/run_traj_checks.py` | Trajectory acceptance harness: mesh C0, C1–C10. `--selftest` checks the harness itself. |
| `tests/SpinWaves/oam_lswt/oracle_honey.py`, `mkhoney.py` | Independent LSWT + Wilson-loop oracle; FM honeycomb with Haldane DMI inputs. |
| `tests/SpinWaves/oam_lswt/oracle_bridge.py` | Independent all-site WLS and exact gauge-invariant k-space oracle for the B5.4 particle-field bridge; `--validate` checks `l - 2F_n`. |
| `tests/SpinWaves/oam_lswt/run_lswt_checks.py` | LSWT acceptance harness: Goldstone, cell invariance, D=0 run, oracle. |
| `tests/SpinWaves/oam_lswt/run_bridge_checks.py` | B5.4 diagnostic: generates LSWT-referenced honeycomb packets, measures frequency and trajectory OAM, and reports the production-mesh diagnostic separately from the settled bridge oracle. |
| `reference/lswt_reference_fixes.diff` | Verified worked example for the DM loop, `diamag_eps` and polar q. |

## Harness status, verified before delivery

| Harness | `b60b762` | Reference patch | `--selftest` |
|---|---|---|---|
| `run_lswt_checks.py` | 4/4 FAIL | 4/4 PASS | — |
| `run_traj_checks.py` | 18/18 FAIL (4 mesh + 14 OAM; no `oam_traj` output) | — | 13/13 PASS (C0 and C10 need a binary) |
| `run_bridge_checks.py` | — | trajectory frequency plus independent particle-field bridge oracle | — |

Run both from a scratch directory; they create run directories in the current directory.

## Dispatch order

```
G0  maintainer: settle C4 (weighting default) and read the contract
 |
 +-- Track A (source/SpinWaves) -------+   Track B (source/Measurement, Input, drivers)
 |   A1  Sonnet  [OAM-A1]              |   B1  Sonnet  [OAM-B1]
 |   A3  Sonnet  [OAM-A3]              |   B2  Sonnet  data only  -> G1 maintainer
 |   A5  Opus    sign pin -> GA        |   B3  Opus    [OAM-B3]  (G2 only if tolerances)
 |                                     |   B4  Sonnet  [OAM-B4]
 +-------------------------------------+
                     |
                    B5  Sonnet (CI, pyswatter prep) + Opus (bridge)  -> G3 maintainer
                    B6  Sonnet  [OAM-B6] removal
```

Tracks A and B touch disjoint files and can run in parallel. The bridge test (B5.4) needs both tracks done and the GA decision.

## Routing

- **Opus:** B3 (kernel and conventions restatement), B5.4 (bridge derivation), A5 (literature mapping).
- **Sonnet:** everything mechanical with an oracle attached — A1, A3, B1, B2, B4, B5.1–B5.3, B6.

## Gate status

- **G1 — accepted (2026-09-29).** `HeisStripe` with an empty 2D mesh is
  classified as zero/N/A chirality; the current binary returns finite zero.
- **G3 — accepted (2026-09-29).** B5.1–B5.4 are signed off. B5.2 accepts
  `ell+1` and `w_area`, with `shift_b` N/A for current pyswatter; B5.3 has
  the positive counter-clockwise sign; B5.4 passes the particle-field bridge
  oracle using `l - 2F_n`. The production per-sublattice trajectory value
  remains a separately labelled diagnostic.

## Change from the 18 Sep instructions

Agents now make one local commit per work package on `feat/oam-refactor`, and must leave no untracked files. They still never push, merge or touch `v6.1.0rc`. The earlier "never commit" rule is how the WP6 fragment ended up orphaned in an unstaged working tree. If you'd rather keep working-tree-only, change rule 1 in the preamble. Everything else is unaffected.
