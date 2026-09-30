# OAM finalization prompt, round 4 (after the audit of `877fcbb`)

| # | Prompt | Model | Commits |
|---|---|---|---|
| U5 | Finalize: B5.5 bridge check, `auto` default, ensemble aggregation, GPU copy bug, strictness, smoke tests, docs and cleanup | Opus | `[OAM-U5a]` … `[OAM-U5g]` |

- `reference/b55_band_packet_reference.py` is the audit script behind U5a. Run it from `tests/SpinWaves/oam_lswt` with a binary path; it takes about 30 s.
- **Maintainer decisions** (2026-09-29):
  - `oam_gradient` defaults to `auto` (spectral where possible, otherwise fem);
  - ensembles are aggregated norm-weighted (new clause C18).
- **Order:**
  - U5a first: it pins the physics before the default changes.
  - U5b and U5c touch the harness defaults, so run them one after the other.
  - U5d and U5e are independent.
  - U5f and U5g come last.
