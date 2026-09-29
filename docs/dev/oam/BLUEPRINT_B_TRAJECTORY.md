# Blueprint B — Trajectory OAM (replace the real-space pseudo-OAM)

**Target:** `UppASD/UppASD`, branch `v6.1.0rc` at `b60b762`. Line numbers refer to that commit.
**Contract:** `CONVENTIONS_OAM.md`, Part I (C1–C8) and Part III.
**Harness:** `tests/SpinWaves/oam_lswt/run_traj_checks.py` (with `oracle_traj.py`, `mkfixture_traj.py`).
- `--selftest` is GREEN.
- Against `b60b762` all 14 functional checks are RED; the mesh checks C0 are also RED.
- The harness loads prescribed textures via restart file and compares against an independent oracle to 1e-6.

This blueprint replaces the 18 Sep version. Changes since then:
- The output contract, mesh diagnostic line and sampling rule are pinned (C6, C7).
- The `sigma_max` default is 0.6, not 0.35: 0.35 rejected a well-localised l = 2 vortex.
- Commit policy is one commit per work package on a feature branch.
- The lessons from the partial WP6 diff are written in (see "What went wrong last time").

## What is wrong now

The active routine `oam_tri_improved` (`topology.f90` 605–671) computes `Σ_t ρ_t Ω_t A_t`. `Ω_t` is already the solid angle through triangle t, so multiplying by `A_t` integrates area twice, and `ρ_t = <|m−S0|²>` makes it quartic in amplitude. Measured against `λ = l`:
- it tracks only the sign of l;
- the magnitudes run 1 : 1.53 : 2.35 for l = 1, 2, 3;
- it scales as amplitude⁴;
- it returns exactly 0 for a texture loaded at start-up, because it measures deviation from the initial state.

The phase-winding variants `oam_tri` and `oam_tri_phase` sum vorticity, which sits on a single core triangle.

The mesh (`delaunay_tri_tri`, 329–386) wraps vertex indices periodically but uses unwrapped coordinates. On 41×41 the total area is 6400 against a true 1600, and 4.8% of the triangles carry 75% of the area weight. This affects every area-weighted consumer.

Plumbing defects:
- `oam_step` and `oam_buff` are never parsed (only `do_oam` is, `inputhandler.f90:1457`).
- The header is written only if `mstep <= 1`, which a restart never reaches (`sd_driver.f90:504` sets `mstep = rstep+1`).
- `S0_arr` is captured at setup, before the initial phase (`uppasd.f90:1311`).
- `mu_arr`, `Lz_i`, `Lz_t` are dead.
- `deallocate` in `calculate_oam(...,2)` has no guard, and `flush_measurements` is called twice by `sld_driver.f90`.
- `prn_averages.f90:1125` divides by `n_Lz_cavg` unguarded.

## What went wrong last time

A partial diff of WP6 surfaced as unstaged changes. Its base blob (`9d9fb32`) was never pushed, and it didn't compile on its own. It removed symbols still used in `measurements.f90`, `uppasd.f90` and `prn_averages.f90`, and it referenced `oam.f90` and `OAM_QUESTIONS.md`, which weren't in the diff. So:

1. **Every work package ends in a commit on the feature branch** (see the preamble), and the whole tree must build.
2. **New files are committed**: `git status --short --untracked-files=all` must be empty at the end of a WP.
3. **Source comments describe the code, not the process.** No "removed in WP6", no references to question logs. One pointer comment at most.
4. **Legacy keys become deprecated aliases** (C15), not orphaned variables fixed at 'N'.

---

## B1 — `source/Measurement/mesh2d.f90` [Sonnet]

**Module state:**

```fortran
integer                    :: nsimp = 0
integer, allocatable       :: simp(:,:)          ! (3,nsimp), counter-clockwise
real(dblprec), allocatable :: tri_area(:)        ! (nsimp) > 0
real(dblprec), allocatable :: site_area(:)       ! (Natom) A_i = sum_{D in i} A_D / 3
real(dblprec), allocatable :: site_wsum(:)       ! (Natom) sum_{D in i} A_D
real(dblprec), allocatable :: grad_b(:,:), grad_c(:,:)   ! (3,nsimp) FEM d/dx, d/dy
integer, allocatable       :: site_tri_ptr(:), site_tri_idx(:)  ! CSR site -> triangles
```

**Tasks.**
1. **Move and fix.** Move `delaunay_tri_tri` here as `mesh2d_build(N1,N2,N3,NA,coord,C1,C2,C3,BC1,BC2,BC3)`.
   - Build each triangle from minimum-image vertex positions relative to vertex 1, using the supercell vectors `N1*C1`, `N2*C2`.
   - Along an open direction, drop the wrap cells.
   - The quad-diagonal choice must use the minimum-image positions.
   - `coord` is never written.
2. **Orientation.** Use the signed area `2A = (x2−x1)(y3−y1) − (y2−y1)(x3−x1)`; swap vertices 2 and 3 when it is negative. Count degenerate triangles (`A < 1e-12`).
3. **FEM coefficients at setup:** `b = (y2−y3, y3−y1, y1−y2)/(2A)`, `c = (x3−x2, x1−x3, x2−x1)/(2A)`, in minimum-image local coordinates.
4. **CSR adjacency** site → triangles, so the gradient can be a gather: OpenMP-safe without atomics and bitwise reproducible.
5. **`mesh2d_report()`** prints exactly the C7 line.
6. **Callers.** `chirality_tri`, `pontryagin_tri*`, the skyrmion-number path and `print_triangulation_mesh` use `mesh2d`. `topology.f90` keeps its physics routines.
7. **Hygiene.** `memocc` on every allocate and deallocate (`topology.f90` has none today). Guard re-allocation with `allocated()`.

**Acceptance.**
- Check C0 in `run_traj_checks.py` (square and hex, periodic and open) is GREEN. The OAM kernel isn't needed for it.
- Add `tests/Measurement/test_mesh2d` (Python driving the binary, or a Fortran unit driver) covering FEM exactness: for `psi = α + βx + γy`, the site gradient equals `(β, γ)` to 1e-12 on every site of a periodic mesh.
- **Stop at gate G1 (B2).**

## B2 — Regression on the mesh fix [Sonnet runs, maintainer decides]

Run the skyrmion-number and chirality measurements (`skyno T`, `do_chiral Y`) on `tests/kagome`, `tests/HeisStripe`, `tests/triang2D` and one skyrmion example, before and after B1, with identical seeds. Tabulate the deltas in `docs/dev/oam/OAM_QUESTIONS.md`.

- **Expected:** the skyrmion number is unchanged. It sums `Ω_t` with no area weight, and the wrap triangles' solid angles are unaffected by the coordinate fix.
- **If it moves, stop.** That is a second defect, not an improvement.
- Chirality may change by about the fraction of wrap triangles (`nsimp` changes for open directions).

**Gate G1:** the maintainer signs off the deltas before B3.

## B3 — `source/Measurement/oam.f90`, module `orbital_angular_momentum` [Opus]

**Public:** `oam_init`, `oam_sample`, `oam_flush`, input variables.

**`oam_init`** (called after the initial phase):
- build the C8 frame and refuse non-collinear states;
- resolve the default `oam_origin`;
- open no file yet.

**`oam_sample(mstep, rstep, emom, mmom, atype)`**, per ensemble k:
1. `psi_i = m_i·e_x + i m_i·e_y` (C1, C2).
2. Triangle gradients `Σ_v b_v ψ_v`, `Σ_v c_v ψ_v`; site gradients by CSR gather, weighted by `tri_area/site_wsum`.
3. Weights `w_i` (C4); drop sites with `site_wsum == 0` (C5).
4. `R` by the C3 rule; `sigma_psi`.
5. `λ` about the origin and about `R` (C3); `N_m` and the derived columns (C6); guards (C5).
6. Average over ensembles *after* computing the per-ensemble λ.
7. Buffer the row; flush every `oam_buff` rows.

**`oam_flush`:** write the remaining rows. Deallocation is guarded with `allocated()` and `stat=`, and a second call is harmless.

**`g`:** from `oam_gfactor` if set, else `Landeg_glob`.

**Sublattices:** `oam_sublattice` (integer list, default all) restricts every sum. With no filter and more than one sublattice, append the per-sublattice block (C6).

**Header text (C3, C6).** It includes the sentence: *"lambda_L_centroid is referenced to the |psi|^2 centroid R; this removes the drift term (R x P)_z but not the envelope winding l."*

**Minimal plumbing, included in B3** so the harness can run:
- parse `do_oam_traj`, `oam_step`, `oam_buff`, `oam_origin`, `oam_weight`, `oam_sigma_max` in `inputhandler.f90`;
- call `oam_init` after the initial phase (`uppasd.f90`);
- call `oam_sample` and `oam_flush` in `measurements.f90` (lines 200 and 306), passing `rstep` or storing it at init.

The old `calculate_oam` calls are removed from those call sites here. The routines themselves are deleted in B6.

**Acceptance:** `run_traj_checks.py` C1–C9 GREEN; C0 still GREEN. **Stop at gate G2 if any check needs a tolerance change.** Tolerances are part of the contract.

## B4 — Input, call sites, output, naming [Sonnet]

1. **`inputhandler.f90`.** Add `oam_gfactor` and `oam_sublattice`, including the per-sublattice output block (C6). Make `do_oam` a deprecated alias (C15).
2. **`uppasd.f90` 1302–1312.** The mesh build becomes a single `mesh2d_build` call whenever any mesh consumer is on (skyrmion number, chirality, OAM). There's no second build inside consumers.
3. **Buffering.** Check that `oam_buff` buffering matches the `prn_trajectories` pattern, and that a run ending mid-buffer loses no rows.
4. **`prn_averages.f90` 1124–1126.** The cumulant JSON key `"orbital_angular_momentum"` becomes the running mean of `lambda_L_centroid` over finite samples; write `null` when there are none.
5. **`prn_topology.f90` 248–251.** Defaults move to the new module; `oam_step` defaults to 100.
6. **Docs.** Document all keys, both OAM paths and the C14 distinction in `docs/` (short section) and in the `inpsd.dat` keyword reference.

**Acceptance:** harness C10 GREEN (legacy alias with a "deprecated" message); full harness GREEN; `tests/` regression suite unchanged.

## B5 — Validation beyond the harness

**B5.1 — CI** [Sonnet]. Register `run_traj_checks.py` as a CTest test labelled `oam-traj`.

**B5.2 — pyswatter cross-validation** [maintainer runs, Sonnet prepares]. Take any harness case directory, e.g. `traj_checks/ell+1`, and run:

```
pyswatter-animate spin-oam-balance coord.oamtest.out restart.oamtest.out \
    --lz-integration site --output ref.csv
```

The column layouts match pyswatter's extended position format and its moment format. Strip the `#` header block first if pyswatter's reader rejects it.

Agreement is required to 1e-3 on `lambda_L`, once the origin is matched (`--shift`, or `oam_origin` set to pyswatter's origin). This check is two-way: a disagreement may be a pyswatter default. Report it; don't adjust the Fortran.

**B5.3 — Sign cross-check** [Sonnet]. Take a single vortex written with `psi ∝ e^{+iφ}` in the C1 frame, i.e. `m_x + i m_y` winding counter-clockwise. Its `lambda_L` must be positive in both UppASD and pyswatter.

**B5.4 — Bridge test** [Opus; requires Blueprint A merged]. Two-sublattice honeycomb with Haldane DMI (`tests/SpinWaves/oam_lswt/mkhoney.py`), where `F_n ≠ 0`:
1. Run `do_oam_lswt`; record `F_n(k0)` for a chosen band n and small k0.
2. Build the trajectory initial state from the full particle band packet, including the angular variation of `T_n(k)` around the ring (Holstein–Primakoff, small amplitude), times a broad Gaussian envelope with l = 0. A fixed `T_n(k0)` spinor is only a control and cannot carry the Berry term. Write it with `mkfixture_traj`-style restart files. Run at zero damping and temperature, sampling for several precession periods.
3. In the independent reciprocal-supercell oracle, use the analytic spatial derivative to assert
   `lambda_L_centroid → -2 F_n(k0)/ħ` for l = 0, within the finite-packet tolerance, and verify that the packet frequency is `E_n(k0)/ħ`.
4. Repeat with an l = 1 envelope and assert a shift of +1. The current per-sublattice production mesh is reported as a separate discretisation diagnostic until an all-site bridge observable is implemented; it must not be relabelled as `F_n`.

This is the only place the two observables are compared (C14). The accepted particle-field convention is `l - 2F_n`; a sign or factor mismatch goes to `OAM_QUESTIONS.md` and is not hidden by changing the oracle.

**Gate G3:** the maintainer signs off B5.2–B5.4.

## B6 — Removal [Sonnet]

Only after G3. From `topology.f90` delete:
- `calculate_oam`, `oam_tri`, `oam_tri_improved`, `oam_tri_phase`;
- `make_orthonormal_basis`, `local_frame`, `project_transverse`;
- the ~110 commented-out lines at 760–867;
- the module variables `S0_arr`, `mu_arr`, `Lz_*`, `n_Lz_cavg`, `step_counter`.

Keep `minimal_image_correction` only if something still uses it.

If a fluctuation-weighted topological-charge diagnostic is wanted, it becomes a separate, differently named observable without the `A_t` factor. Ask; don't keep the old code as a hedge.

**Acceptance:**
- clean build with `-Wall` producing no new unused-variable warnings from these files;
- harness GREEN;
- `git grep -n "calculate_oam\|oam_tri\|Lz_csum"` empty outside docs and history.

## Out of scope

- Smooth local frames for non-collinear states (C8, phase 2).
- A GPU implementation.
- The Zeeman term in LSWT (Blueprint A2.4).
