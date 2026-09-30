# R2 — Multilayer mesh and empty-mesh guard · Sonnet · `[OAM-R2]`

**Regression.** `mesh2d_build` indexes only the first z-layer (`i00 = NA*((y-1)*N1+(x-1))+it`), while the old `delaunay_tri_tri` looped over all layers. A two-layer skyrmion texture gives skyrmion number −2 on `b60b762` and −1 on the branch. Trajectory OAM then takes λ from layer 1 but `N_m` from all layers.

**Tasks.**
1. Triangulate every z-layer: add the layer offset `NA*N1*N2*(z-1)` to the site index. `max_tri` grows by a factor of N3. Triangles never cross layers.
2. Keep the old skyrmion-number semantics (sum over layers). If you think per-layer output is needed, log it in `OAM_QUESTIONS.md`; don't implement it.
3. `chirality_tri` divides by `nsimp`. When `nsimp == 0`, print one warning and return zeros. HeisStripe (xz-plane stripe, `ntri= 0`) currently gives NaN.
4. `mesh2d_build`: warn once when `N2 == 1 .and. N3 > 1` (the system isn't in the xy plane).
5. You may add one new check, C11, to `run_traj_checks.py`: a two-layer copy of the `ell+1` case must give the same λ as one layer, with `N_m` doubled. Don't change existing checks or tolerances.

**Acceptance:**
- trajectory harness ALL PASS, including C11;
- regression suite 20/20;
- the two-layer skyrmion test (build it with `mkfixture_traj.py`, `skyno T`, N3 = 2) gives −2.
