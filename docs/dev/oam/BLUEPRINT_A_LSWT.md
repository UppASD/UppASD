# Blueprint A — LSWT Hamiltonian and Fishman magnon OAM fixes

**Target:** `UppASD/UppASD`, branch `v6.1.0rc` at `b60b762` (4 Sep 2026). Line numbers below refer to that commit.
**Contract:** `CONVENTIONS_OAM.md`, Part II (C9–C13) and C16.
**Harness:** `tests/SpinWaves/oam_lswt/run_lswt_checks.py` (with `mkhoney.py`, `oracle_honey.py`). On `b60b762` all four checks fail; with the reference patch applied, all four pass.
**Reference patch:** `reference/lswt_reference_fixes.diff` fixes A1 (DM only), A2.1 and A3.1, and was verified with the harness. Treat it as a worked example, not a patch to apply blindly: it does not cover SA, PD or AMS.

## Why this goes first

The Fishman algebra in `chern_number.f90` is correct: it matches an independent Wilson-loop oracle to 0.3% on the k-points it actually samples, and bands 1 and 2 sum to zero to 2e-16. The numbers are still wrong for every system with nonzero magnon OAM, because of three defects upstream of the formula. The same defects corrupt dispersions, S(q,ω), Chern numbers and κxy, so this blueprint is a release fix in its own right. Blueprint B's bridge test (B5.4) depends on it.

It touches only `source/SpinWaves/` and can run in parallel with Blueprint B1–B4, which touch `source/Measurement/`, `source/Input/` and the drivers.

## Evidence

The test model is an FM honeycomb with NN J = 1 mRy and Haldane NNN DMI D = 0.1 mRy along z. The oracle is a separate LSWT with a Wilson-loop Berry phase.

| Defect | Observation on `b60b762` | After fix |
|---|---|---|
| DM on the wrong bonds | Goldstone gap 0.091 meV (D=0.1), 0.81 meV (D=0.3); imaginary NN term at Γ; Chern numbers sign-flipped | ω(Γ) = 1.1e-4 meV (regularisation only) |
| Reduced coordinates in the polar mesh | F₁(k=2.34) = 8.40e-3 with C₂=(0.5, 0.866) but 0.163 with C₂=(1.5, 0.866), for the same lattice | Both 7.180e-3; oracle 7.205e-3 |
| `diamag_eps` never initialised | Isotropic FM aborts: `ERROR STOP diamag: invalid paraunitary eigenvectors`, paraunitarity error 6.8e-3 at Γ | Runs; error 5e-15 |

---

## A1 — Coupling vectors indexed by their own neighbour list (C12)

**Files:** `source/SpinWaves/diamag.f90` (`setup_Jtens_q`, 1695–1841), `source/SpinWaves/ams.f90`.

**Defect.** The loop at 1740 runs over the exchange list (`ham%nlist`, `ham%nlistsize`), but reads:
- `ham%dm_vect(:,j,ih)` at 1756;
- `ham%sa_vect(:,j,ih)` at 1761;
- `ham%pd_vect(:,j,ih)` at 1766.

Each of those arrays belongs to its own list (`dmlist/dmlistsize`, `salist/salistsize`, `pdlist/pdlistsize`, allocated in `hamiltoniandata.f90` 267–272, 299–303, 398–402). When the lists differ, couplings land on the wrong bonds, some are dropped, and `j` can run past `max_no_dmneigh`. The Kagome Chern example is unaffected only because its DM and exchange bonds coincide; its output is bit-identical with and without the fix.

**Fix.**
1. Keep the exchange loop for the isotropic `ncoup` only.
2. Add three loops after it, one per enabled coupling (`do_dm`, `do_sa`, `do_pd`). Each takes `ja` from its own list, recomputes `jat`, `dist` (`f_wrap_coord_diff`), `R_n` (`find_R`) and `FTfac`, and accumulates into `Jtens_q(:,:,ia,jat,iq)` with the sign already used:
   - DM: `J_n = -D_n` with `D_n = dm2tens(dm_vect)`;
   - SA: `+sa2tens(-sa_vect)`;
   - PD: `+pd2tens(-pd_vect)`.
3. In `ams.f90`, the DM loop at 322 correctly iterates `dmlist` but takes the distance from `ham%nlist(j,mutemp)` at 325. Use `ham%dmlist(j,mutemp)`. Check the other DM loops in the file for the same pattern.
4. Check `setup_Jtens2_q` (`do_jtensor=1`) and state in the report whether it has the same issue.

**Acceptance.**
- `run_lswt_checks.py` check [1] passes: ω(Γ) < 1e-2 meV for the isotropic FM with DMI.
- The Kagome Chern example output is unchanged: `bphase*.out` bitwise. Run it shortened with `ip_mode N`, `Nstep 1`, `ncell 12 12 1`, `kgrid 30 30 1`, `do_sc N`, which was verified identical before/after.
- Add one SA case to `tests/SpinWaves`: an FM chain with SA on a second-neighbour list different from the J list. The Goldstone mode must stay gapless when the SA tensor is traceless and symmetric about z. If it's unclear what the SA term should do physically, escalate rather than guess.

## A2 — Diagonaliser robustness

**File:** `source/SpinWaves/diamag.f90`.

**A2.1 — Restore the regularisation default.** At line 110 `!diamag_eps=-1.0_dblprec` is commented out, so the module variable is uninitialised. gfortran gives 0, so the intended `dia_eps = 1e-6` branch at 502–506 is never taken. At an exact Goldstone point the Cholesky factor is then singular to machine precision and T loses paraunitarity. Uncomment the line. In the report, list every `setup_diamag` call path so the default is guaranteed on every route into `diagonalize_quad_hamiltonian`, including `do_chern` without `do_diamag`.

**A2.2 — Don't hide instabilities.** Line 541 shifts the spectrum by `-minval(eig_val)` whenever K has a negative eigenvalue. This hides genuine instabilities: flipping the sign of the anisotropy constant gave the same Γ gap and no warning.
- Warn when `-λmin > 1e-8 * max|λ(K)|`, giving iq and the magnitude.
- When `require_paraunitary` is true (Chern/OAM callers), make it `error stop` unless the new key `nc_allow_unstable Y` is set.
- The tiny negative eigenvalues from roundoff at Goldstone points (about −1e-15 relative) must stay silent.

**A2.3 — Fallback units.** `fallback_bosonic_diag` (658–758) solves on `h_in`, while the main path uses `K = -h_in*fcinv`. Its eigenvalues therefore come out with the wrong sign and scale (factor −1/fcinv ≈ −470), and the "positive" branch it picks is the hole branch. Pass it the same shifted `K`. The path is effectively unreachable after A2.1; one unit test that forces it, by a flag or a direct call, is enough.

**A2.4 — Zeeman term (document only).** `hfield` does not enter the LSWT Hamiltonian: 10 T left every energy unchanged. Don't implement it here. Print one line when `hfield ≠ 0` together with `do_diamag`/`do_chern`, and add a sentence to the docs. Record implementation as an item in `OAM_QUESTIONS.md`.

## A3 — Fishman polar mesh and outputs

**File:** `source/SpinWaves/chern_number.f90`.

**A3.1 — Cartesian q (C11).** Line 931 passes `polar_cartesian_to_reduced(...)`, i.e. `q_i = k·C_i/2π`, but `setup_Jtens_q`/`setup_ektij` treat q as Cartesian (`q·2π·dist`). For any non-orthonormal cell, each ring becomes a sheared ellipse. Replace it with `polar_q(:,iq) = (/ kr*cos(phi), kr*sin(phi), 0 /) / (2π)`. Then:
- delete `polar_cartesian_to_reduced` (1150) and `fishman_cartesian_to_reduced` (1139);
- delete `fishman_reduced_to_cartesian` (1161), which is dead code;
- keep `fishman_reciprocal_basis` and `fishman_inscribed_radius`, which are correct;
- fix the comments at 1110–1111 and 1136–1138, which assert that UppASD passes reduced q.

The shipped `test_cartesian_reduced_round_trip_for_oblique_cell` passed throughout, because it tests the conversion, not what consumes it. Replace it with the harness's cell-description invariance check.

**A3.2 — Sign comment (C9).** Rewrite 1099–1102 and the matching comment in `fishman_ring_oam`: the −½ is Fishman Eq. 11 with T = X⁻¹; the Fourier sign plays no role for ring averages in 2D.

**A3.3 — Outputs and naming (C10, C15).**
- `f_oam.<simid>.out` → `oam_lswt.<simid>.out`; `f_oam_diagnostics` → `oam_lswt_diagnostics`.
- The pointwise `oam_k` output (450–464) is gauge-dependent: write it only when `oam_lswt_pointwise Y`, appended to the diagnostics file, with its header kept.
- `do_magnon_oam` → `do_oam_lswt` in `read_parameters_chern_number` (1309), with the deprecated alias.
- Put the C13 sign label in the `oam_lswt` header.

**A3.4 — Cost (optional, low priority).** For `flag /= 0`, `setup_tensor_hamiltonian` diagonalises 3·nq points and allocates `S_prime(2NA,2NA,3,3,6nq)`. With default polar settings that is about 2 GB at NA = 12.
- Skip `iq > nq` when `norm2(diamag_qvect) == 0` (the ±Q copies are then identical).
- Don't transform or allocate `S_prime` when `flag /= 0`: make the argument optional in `diagonalize_quad_hamiltonian`.
- Measure the saving and report it.

## A4 — Tests

1. Use the harnesses in `tests/SpinWaves/oam_lswt/`. Register `run_lswt_checks.py` as a CTest test (binary path from `${TestBinary}`) under `RUN_REG_TESTS`, labelled `oam-lswt`.
2. In `test_fishman_oam.py`:
   - keep the algebra tests;
   - replace the string-grep test `test_fortran_oam_is_not_the_berry_flux_proxy` with an assertion that the polar mesh builds q as `k/(2π)`;
   - add the Berry-flux identity (C10) as a pure-Python test using `oracle_honey.py`.
3. The harness must be RED on `b60b762` and GREEN after A1–A3. Report both outputs.

## A5 — Absolute sign of F (gate GA, maintainer)

Opus prepares and the maintainer decides:
1. Set up Fishman's FM honeycomb in UppASD units, with rings extended to K (`f_oam_kmax` = |K|; the warning about exceeding the inscribed radius is expected).
2. Report O₁,av(k) for both signs of D next to the published 0.236ħ peak for d = 0.1.
3. The maintainer maps d to D/J and fixes the sign label in the header.

Do not change the −½ to make a number match. If neither sign reproduces the magnitude, that is a finding to escalate.

## Release note (draft)

> LSWT: DM, symmetric-anisotropic and pseudo-dipolar couplings were attached to exchange-list bonds when their neighbour lists differed. Dispersions, S(q,ω), Chern numbers and κxy for such inputs change. Inputs whose DM bonds coincide with the exchange bonds (e.g. the Kagome example) are unaffected. Fishman F_OAM rings were ellipses for non-orthonormal cells. `do_magnon_oam` is now `do_oam_lswt`.
