# R7 — Remove pseudo-OAM · Sonnet · `[OAM-B6]`

After G3. Run `docs/dev/oam/prompts/B6_removal.md` unchanged.

Facts for this branch:
- `calculate_oam`, `oam_tri`, `oam_tri_improved`, `oam_tri_phase`, `make_orthonormal_basis`, `local_frame` and `project_transverse` are all still in `topology.f90` and no longer called;
- `S0_arr`, `mu_arr`, `Lz_*` and `step_counter` remain.

**Acceptance:** `git grep -n "calculate_oam\|oam_tri\|Lz_csum\|S0_arr" source/` is empty; both harnesses ALL PASS; regression suite 20/20.
