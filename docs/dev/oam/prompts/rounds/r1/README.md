# OAM follow-up prompts (after audit of `59ea123`)

Each prompt is thin: it names the task, the evidence and the acceptance check, and relies on the package docs for everything else. Paste `00_PREAMBLE.md` (from the original package) above each one. R0 puts that preamble in the repo.

| # | Prompt | Model | Depends on | Commit |
|---|---|---|---|---|
| R0 | Commit the planning docs | Sonnet | — | `[OAM-R0]` |
| R1 | Recover or redo A1/A2 | Sonnet | R0 | `[OAM-A1]` |
| R2 | Multilayer mesh + empty-mesh guard | Sonnet | R0 | `[OAM-R2]` |
| R3 | Oracle C8 frame | Sonnet | R0 | `[OAM-R3]` |
| R4 | A5 rerun with corrected D mapping | Sonnet | R1 | `[OAM-A5b]` |
| R5 | Small cleanups | Sonnet | R2 | `[OAM-R5]` |
| R6 | Bridge derivation (no code) | Opus | R1, R4 | `[OAM-R6]` |
| R7 | Remove pseudo-OAM | Sonnet | G3 | `[OAM-B6]` |

R1, R2 and R3 touch different files and can run in parallel.

**Maintainer-only:** confirm C4 (`site` default); close GA after R4; run `docs/dev/oam/pyswatter/run_pyswatter_checks.sh`. The `ell+1` and `w_area` comparisons are the B5.2 checks. Treat `shift_b` as N/A/skip for current pyswatter because its periodic minimum-image mesh is not implemented; do not adjust Fortran to fit that mismatch.
