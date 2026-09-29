# R0 — Commit the planning docs · Sonnet · `[OAM-R0]`

The contract, blueprints and prompts were never committed. Copy `CONVENTIONS_OAM.md`, `BLUEPRINT_A_LSWT.md`, `BLUEPRINT_B_TRAJECTORY.md`, `README.md`, `prompts/` and `reference/` from the package into `docs/dev/oam/` on `feat/oam-refactor`.

Leave the oracles where they are (`tests/SpinWaves/oam_lswt/`). Fix any path in the docs that still points to `docs/dev/oam/oracles/`.

No source changes. Commit, report `git show --stat`.
