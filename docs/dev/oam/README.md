# UppASD OAM fix package

Regenerated 27 Sep 2026 against `v6.1.0rc` at `b60b762`. It replaces the 18 Sep blueprint, none of which reached the remote.

Copy this folder into the repository as `docs/dev/oam/` and commit it on `feat/oam-refactor` before dispatching the first prompt. The prompts refer to that path.

## Contents

| File | Purpose |
|---|---|
| `CONVENTIONS_OAM.md` | Binding physics contract, C1–C16. Two clauses are OPEN for you (C4, C13). |
| `BLUEPRINT_A_LSWT.md` | Fixes upstream of the Fishman magnon OAM: DM/SA/PD neighbour lists, `diamag_eps`, polar mesh, outputs, sign pin. |
| `BLUEPRINT_B_TRAJECTORY.md` | Replaces the real-space pseudo-OAM with the pyswatter formulation: mesh, kernel, I/O, validation, removal. |
| `prompts/00_PREAMBLE.md` | Standing rules; paste above every prompt. |
| `prompts/*.md` | One paste-ready prompt per work package. |
| `tests/SpinWaves/oam_lswt/oracle_traj.py` | Independent trajectory-OAM reference. `python3 oracle_traj.py` runs its self-test. |
| `tests/SpinWaves/oam_lswt/mkfixture_traj.py` | Writes UppASD run directories with a prescribed texture, loaded by restart file. |
| `tests/SpinWaves/oam_lswt/run_traj_checks.py` | Trajectory acceptance harness: mesh C0, C1–C10. `--selftest` checks the harness itself. |
| `tests/SpinWaves/oam_lswt/oracle_honey.py`, `mkhoney.py` | Independent LSWT + Wilson-loop oracle; FM honeycomb with Haldane DMI inputs. |
| `tests/SpinWaves/oam_lswt/run_lswt_checks.py` | LSWT acceptance harness: Goldstone, cell invariance, D=0 run, oracle. |
| `reference/lswt_reference_fixes.diff` | Verified worked example for the DM loop, `diamag_eps` and polar q. |

## Harness status, verified before delivery

| Harness | `b60b762` | Reference patch | `--selftest` |
|---|---|---|---|
| `run_lswt_checks.py` | 4/4 FAIL | 4/4 PASS | — |
| `run_traj_checks.py` | 18/18 FAIL (4 mesh + 14 OAM; no `oam_traj` output) | — | 13/13 PASS (C0 and C10 need a binary) |

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

## Decisions waiting for you

1. **C4.** Default weighting `site` (proposed, matches the magnon count and `spin-oam-balance`) or `area` (pyswatter `lz` default).
2. **C13 / GA.** Absolute sign of F against Fishman 2023, after A5 reports.
3. **G1.** Skyrmion-number and chirality deltas from the mesh fix.
4. **G3.** pyswatter agreement and the bridge test.

## Change from the 18 Sep instructions

Agents now make one local commit per work package on `feat/oam-refactor`, and must leave no untracked files. They still never push, merge or touch `v6.1.0rc`. The earlier "never commit" rule is how the WP6 fragment ended up orphaned in an unstaged working tree. If you'd rather keep working-tree-only, change rule 1 in the preamble. Everything else is unaffected.
