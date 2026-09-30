# OAM follow-up prompts, round 2 (after audit of `49f1b08`)

Paste `docs/dev/oam/prompts/00_PREAMBLE.md` above each prompt. The archived round is `docs/dev/oam/prompts/rounds/r2/`.

| # | Prompt | Model | Commit |
|---|---|---|---|
| T1 | Fix `zggev` workspace | Sonnet | `[OAM-T1]` |
| T2 | Correct BRIDGE §3 and add the validity statement | Sonnet | `[OAM-T2]` |
| T3 | Spectral gradients and `oam_axis` | Opus (a), Sonnet (b) | `[OAM-T3a]`, `[OAM-T3b]` |
| T4 | Revert the C3 fixture change and trim the LSWT header | Sonnet | `[OAM-T4]` |

- **Order:** T1, T2 and T4 are independent and small. T3 goes last: it touches `oam.f90`, the contract and the oracle, and its validity numbers feed back into T2's text.
- **Gate:** the maintainer signs off the T3a contract clause (C17) before T3a code is merged.
