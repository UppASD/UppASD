# B3 — Trajectory OAM kernel  ·  Opus  ·  commit `[OAM-B3]`  ·  gate G2 if tolerances are touched

Implement Blueprint B section **B3**, including its minimal plumbing. Paste the shared preamble first. Requires G1 sign-off.

## Before coding

Write into `OAM_QUESTIONS.md`, as your first entry, a restatement **in your own words** of:
- C1 (why `m_x + i m_y` and what a positive λ means);
- C2 (why there is no `m − S0`);
- C3 (what `lambda_L_centroid` removes and what it does not).

Derive the origin dependence yourself: `λ(x0') − λ(x0) = ((x0 − x0') × P)_z / N`, where `P = Σ Im[ψ* ∇ψ] w`. Say why a stationary vortex shows no origin dependence.

If your restatement disagrees with the contract anywhere, stop there.

## Implementation notes

- **Per-sample cost** is O(nsimp) complex multiply-adds using the precomputed `grad_b`/`grad_c`. There are no linear solves inside the time loop.
- **Gradient to sites:** gather via CSR, weighted by `tri_area / site_wsum`. Use `!$omp parallel do` over sites; no atomics.
- **Centroid:** circular mean along periodic axes, in reduced coordinates of the supercell `(N1*C1, N2*C2)`. Minimum-image lever arms to R; raw lever arms to the fixed origin.
- **Ensembles:** compute λ per ensemble, then average.
- **Guards:** exactly C5. `NaN` via `ieee_value(0._dblprec, ieee_quiet_nan)`.
- **Header:** exactly C6, including the centroid sentence.
- **Sampling:** `mod(mstep - rstep - 1, oam_step) == 0`.

## Acceptance

- `UPPASD=... python3 run_traj_checks.py`: C0–C9 GREEN. C10 (legacy alias) belongs to B4 and may FAIL.
- Run `--selftest` too and include it in the report.
- A single-thread and a 4-thread run give identical `oam_traj` rows.

**Commit `[OAM-B3]`, report, stop.** If any check could only pass by changing a tolerance or an oracle, that is gate G2: stop and log it.
