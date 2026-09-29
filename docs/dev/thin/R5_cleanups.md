# R5 — Small cleanups · Sonnet · `[OAM-R5]`

1. **`mesh2d_report`.** Print `cell_area` as the supercell area, `N1*N2*|C1×C2|`, so it compares directly with `total_area` (C7). It currently prints the unit-cell area.
2. **Warnings.** Remove the unused dummy argument `atype` from `oam_sample` and its caller, and unused `use` imports in `chirality_tri` and `calculate_oam`. Rebuild with `-Wall`: no warnings from `oam.f90` or `mesh2d.f90`.
3. **Alias.** The `do_oam` alias only ever switches trajectory OAM *on*. `do_oam N` must not override an earlier `do_oam_traj Y`.
4. **Mesh trigger.** In `uppasd.f90`, build the mesh for `skyno=='T'` but not for `skyno=='Y'`, which uses the stencil path.

**Acceptance:** both harnesses ALL PASS; regression suite 20/20; clean `-Wall` for the touched files.
