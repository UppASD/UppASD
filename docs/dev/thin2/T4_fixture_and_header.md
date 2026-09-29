# T4 — Revert the C3 fixture change and trim the LSWT header · Sonnet · `[OAM-T4]`

1. **C3 fixture.** In `run_traj_checks.py`, restore C3's `boosted()` to the original: no `psi - psi.mean()`, `m_z` taken from the vortex, and no `"timestep": "1e-30"` key. The change was unauthorised, and it removed the frame-tilt coverage. With the R3 oracle the original fixture passes: the audit measured Fortran against oracle at 1e-14.
   - Keep the 24.16E restart precision in `mkfixture_traj.py`.
2. **LSWT header.** In `chern_number.f90`, replace the two C13 header lines with one generic line:

   ```
   # Sign: O_n = -1/2 Im[T^dagger eta dT/dphi], T = X^-1 (Fishman); validated against PRB 107, 214434, see docs/dev/oam/sign_pin.
   ```

   The honeycomb-specific 0.228/0.236 numbers belong in `sign_pin/A5_report.md` only.

**Acceptance:** trajectory harness ALL PASS (C3 on the original fixture); `run_lswt_checks.py` ALL PASS; the `oam_lswt` header contains no model-specific numbers.
