# OAM refactor status

Updated 2026-09-29 for U5.

## Validated

- The trajectory harness passes C0–C15 and C18, including periodic/open and
  non-FFTW `oam_gradient auto` resolution, norm-weighted ensemble λ columns,
  arithmetic ensemble `N_m`/`dSz_hbar`, aggregate `Lz_tot_hbar`/`balance`,
  centroid diagnostics, and the one-time exclusion warning.
- The independent LSWT bridge passes its frequency and particle-field checks.
  B5.5 is registered as `oam-bridge-band`; its spectral values agree with the
  exact-gradient oracle to `4.46e-6`, the l-shift to `4.52e-6`, and LSWT to
  `1.35e-3`. The reported FEM values are the expected short-wavelength-biased
  control.
- The CPU builds pass the OAM CTest label on both FFTW and non-FFTW builds;
  non-FFTW B5.5 and explicit spectral probes skip/refuse with their documented
  notices. The harness self-test passes, and the Resaro regression suite
  reports 20 tests performed with 0 failures.
- The CUDA measurement path was verified by inspection: its existing
  `fortran_do_measurements` copy gate now receives the public
  `oam_sample_due(mstep)` trigger. No CUDA toolchain was available for a
  numerical GPU comparison.

## Scope

The reference is a collinear ferromagnet. Spectral gradients require periodic
in-plane cells and field content strictly inside the first Brillouin zone;
coherent, localised packets are the intended observable. MKL-FFT/non-FFTW
builds resolve `oam_gradient auto` to FEM. Explicit `spectral` remains strict
and refuses unsupported cells or builds.

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

- The LSWT magnitude gap remains `O₁,av(K) = 0.228` versus `0.236`
  published (−3.4%); the sign is confirmed.
- The carrier-referenced gradient proposal remains open.
