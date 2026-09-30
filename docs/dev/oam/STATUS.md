# OAM refactor status

Updated 2026-09-30 for U6.

## Validated

- The trajectory harness passes C0–C16 and C18, including periodic/open and
  non-FFTW `oam_gradient auto` resolution, periodic circular C14c ensemble
  centroids, norm-weighted ensemble λ columns, arithmetic ensemble
  `N_m`/`dSz_hbar`, aggregate `Lz_tot_hbar`/`balance`, centroid diagnostics,
  the one-time exclusion warning, and the Debug-bounds C16 dilute refusal.
- The independent LSWT bridge passes its frequency and particle-field checks.
  B5.5 is registered as `oam-bridge-band`; its spectral values agree with the
  exact-gradient oracle to `4.46e-6`, the l-shift to `4.52e-6`, and LSWT to
  `1.35e-3`. The reported FEM values are the expected short-wavelength-biased
  control.
- The CPU builds pass the OAM CTest label on both FFTW and non-FFTW builds;
  non-FFTW B5.5 and explicit spectral probes skip/refuse with their documented
  notices. The harness self-test passes, and the Resaro regression suite
  reports 20 tests performed with 0 failures.
- The report-only undamped dynamic band check is available as
  `python3 run_bridge_checks.py --band-dynamics`. It is not registered in
  CTest and takes about 30 seconds. It uses the B5.5 honeycomb band packet at
  `k0 = 1.0`, `l = 0, 1`, `mkhoney.write(J=1, D=1)`, `N = 90`, zero damping,
  3001 steps at `1e-16 s`, `oam_step 300`, and `oam_axis 0 0 1`. The expected
  spectral `lambda_L_centroid(t)` remains within `3e-3` of its initial value
  (1.0006 for `l = 1`); the `l = 1` minus `l = 0` shift stays within `1e-3`
  of 1.000 at every sample, and `N_m` stays within `1e-4` relative. FEM has a
  time-independent short-wavelength bias, about 0.880 for `l = 1`. In a
  single band with isotropic `|psi_k|^2`, `⟨∂_φ ω⟩ = 0`, so `lambda` is
  conserved. The Release FFTW run measured spectral `l = 1` from 1.0006016
  to 0.9984446, a maximum l-shift error of `2.246e-5`, and maximum relative
  `N_m` drift of `3.642e-5`; FEM `l = 1` averaged 0.8806816.
- The CUDA measurement path was verified by inspection: its existing
  `fortran_do_measurements` copy gate now receives the public
  `oam_sample_due(mstep)` trigger. No CUDA toolchain was available for a
  numerical GPU comparison.

## Scope

The reference is a collinear ferromagnet. Spectral gradients require periodic
in-plane cells and field content strictly inside the first Brillouin zone;
coherent, localised packets are the intended observable. MKL-FFT/non-FFTW
builds resolve `oam_gradient auto` to FEM. Explicit `spectral` remains strict
and refuses unsupported cells or builds. Trajectory OAM supports full
(non-dilute) lattices only; dilute systems are refused.

## U5f physical smoke runs

The damped clean l=1 packet used `damping = 0.01`, `T = 0`, 2000 steps,
`oam_step = 10`, and `auto` on the FFTW Debug build. Values are sampled every
200 steps; `N_m` is the fourth output column.

| step | `N_m` | `lambda_L_centroid` |
|---:|---:|---:|
| 1 | 9.06924522e-03 | 0.99999997 |
| 201 | 8.13428985e-03 | 0.99211033 |
| 401 | 7.31654978e-03 | 0.98015579 |
| 601 | 6.60003900e-03 | 0.97129693 |
| 801 | 5.97101877e-03 | 0.97020114 |
| 1001 | 5.41767746e-03 | 0.97267985 |
| 1201 | 4.92985695e-03 | 0.97908606 |
| 1401 | 4.49881861e-03 | 0.98191751 |
| 1601 | 4.11704269e-03 | 0.96910203 |
| 1801 | 3.77805568e-03 | 0.94344747 |
| 2001 | 3.47628217e-03 | 0.91881295 |

The thermal run used a square FM restart, `T = 10 K`, damping `0.13`,
2000 steps, `Mensemble = 4`, `oam_step = 10`, and `auto`. It completed without
crashing and printed the expected one-time exclusion warning. The exact
saturated step-1 row is `NaN` under the zero-norm guard; after thermalisation,
origin λ is finite and centroid λ is rejected by the spread guard throughout.

| step | `lambda_L_origin` | `lambda_L_centroid` | `N_m` |
|---:|---:|---:|---:|
| 201 | 0.07969388 | NaN | 11.56634605 |
| 401 | -0.00751175 | NaN | 12.79003903 |
| 601 | -0.20414182 | NaN | 14.15533832 |
| 801 | 0.10494615 | NaN | 14.48687120 |
| 1001 | 0.00701361 | NaN | 15.18929574 |
| 1201 | -0.02988981 | NaN | 15.13122365 |
| 1401 | 0.34766273 | NaN | 16.25817899 |
| 1601 | 0.00314794 | NaN | 16.31254626 |
| 1801 | 0.05523970 | NaN | 16.17549598 |
| 2001 | 0.10182326 | NaN | 16.76702606 |

## Open items

- No numerical GPU comparison yet: run the C12b fixture with `oam_step 10` on
  a CUDA machine and compare with the CPU run to `1e-10`.
- The LSWT magnitude gap remains `O₁,av(K) = 0.228` versus `0.236`
  published (−3.4%); the sign is confirmed.
- The carrier-referenced gradient proposal remains open.
