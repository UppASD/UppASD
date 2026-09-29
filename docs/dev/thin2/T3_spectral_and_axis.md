# T3 — Spectral gradients and `oam_axis`

Two commits. Do (a) first.

## (a) Spectral per-sublattice gradient · Opus · `[OAM-T3a]` · gate on C17

**Why.** The audit showed that an FFT gradient on each sublattice's Bravais grid reproduces the exact λ to 1e-5 at k₀ = 1 and 2 (1.0006 and 1.0408 against exact 1.0006 and 1.0408), where linear FEM gives 0.882 and 0.602.

**First**, draft contract clause **C17** in `CONVENTIONS_OAM.md` and stop for sign-off:
- **Key:** `oam_gradient fem|spectral`, default `fem`.
- **Scope:** spectral requires `BC1 = BC2 = 'P'` and a build with `USE_FFTW`; otherwise `oam_init` refuses with a clear message and produces no numbers. There is no silent fallback.
- **Grid:** each (sublattice `it`, layer `z`) is an `N1×N2` grid in reduced coordinates. The Cartesian k is `(m1·b1 + m2·b2)/N` with `m` folded to `[-N/2, N/2)`, and b = 2π·reciprocal basis of `C1, C2`. The Nyquist component's derivative is set to zero for even N.
- **Semantics:** the gradient is exact for fields band-limited to that fold. The weights, centroid and every output column are unchanged from C3–C6. The `oam_traj` header records the gradient method.

**After sign-off:**
1. Implement in `oam.f90`, following the `#ifdef USE_FFTW` pattern of `fftdipole_fftw.f90`. Use 2D complex plans created once in `oam_init` and destroyed in `oam_flush`. Plans are per grid shape, not per sublattice.
2. **Oracle (authorised addition).** Add `spectral_gradient()` to `oracle_traj.py`, written independently in numpy (`np.fft`), and a `gradient=` argument to `evaluate()`. Don't change the FEM path.
3. **Harness (authorised addition).** Add C12 to `run_traj_checks.py`, skipped with a notice if the binary lacks FFTW:
   - C12a: spectral output matches the spectral oracle to 1e-6 on the `ell+1` and `hex` cases;
   - C12b: for a vortex boosted to k·a = 1.5 on the square lattice, spectral λ matches the analytic λ of that plane-wave sum to 1e-4, and FEM is at least 0.1 away (demonstrating the bias).
4. **Bridge.** Rerun `run_bridge_checks.py` with `oam_gradient spectral` for the production packet. Report whether the production l = 1 minus l = 0 shift and the l = 0 value now meet the l − 2F targets. **Report, don't gate.**

**Acceptance:** C0–C12 ALL PASS (C12 run with an FFTW build); the FEM default is bit-identical to `49f1b08` outputs on all harness cases.

## (b) `oam_axis` · Sonnet · `[OAM-T3b]`

The proposal logged in `OAM_QUESTIONS.md` (R3) is approved. Add clause C8b to the contract:
- **Key:** `oam_axis x y z`.
- **Default:** the normalised average moment at init (current behaviour).
- When set, the axis is normalised and used as `e_z`. The collinearity check uses this axis.
- The header records the axis and whether it was defaulted.

**Tasks.**
1. Implement in `oam_init` and `inputhandler.f90`.
2. Add an `axis=` argument to the oracle's `to_c8_frame` / `evaluate` (authorised).
3. Add harness check C13: the original C3 boosted packet (no mean subtraction) with `oam_axis 0 0 1` matches the oracle evaluated in the lab frame to 1e-6. Without the key, it matches the default-frame oracle.

**Acceptance:** ALL PASS; regression suite 20/20.
