# OAM follow-up prompts, round 3 (after the audit of `3dac8e0`)

Paste `docs/dev/oam/prompts/00_PREAMBLE.md` above each prompt. Commit these files to `docs/dev/thin3/` first.

| # | Prompt | Model | Commit |
|---|---|---|---|
| U1 | Fold spectral k to the Brillouin zone (contract, Fortran, oracle, check C12c) | Opus | `[OAM-U1]` |
| U2 | Finish the T2 docs; quote λ bias, not phase-gradient bias | Sonnet | `[OAM-U2]` |
| U3 | Rerun and report the bridge, with the correct targets | Sonnet | `[OAM-U3]` |
| U4 | Document the zone-boundary (K/M) limitation | Sonnet | `[OAM-U4]` |

- **Order:** U1 first. U3 needs U1. U2 and U4 can run in parallel after U1, because they quote the amended C17 text.
- **One commit per prompt, with its tag.** Don't combine prompts in one commit.
- **Bit-identity comparisons use `CMAKE_BUILD_TYPE=Debug` builds.** Release adds `-Ofast`, and any recompile then drifts by about 1e-12. The audit confirmed this drift is not a code change: 27/27 cases were byte-identical to `49f1b08` at `-O0`.
- **Gate:** none. The C17 amendment text in U1 is approved by the maintainer as written.
