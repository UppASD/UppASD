# CPU-HAM-05 — CPU pair-backend crossover campaign

**Date:** 2026-09-02
**Source revision:** `f6b28c09717c15eb53d456b8f3e9d423`
**Status:** complete for the measured workload and size ranges; no automatic backend policy added.

## Decision summary

- **bcc Fe, medium range (`z = 96`, `N_A = 2`):** retain **DIRECT** as the measured default. Portable **SPARSE** did not produce a robust win; its 20³ result is an isolated high-residual fit and is not treated as a crossover.
- **dhcp Nd, long range (`z = 1338`, `N_A = 4`):** **CONVOLUTION** is the clear steady-state winner at every measured size. With FFTW it is about 8–11× faster than DIRECT on the 16³–25³ ladder and amortizes its fixed cost within approximately 0–14 production steps against the measured alternatives.
- **Short-range scalar-J control (`z = 6`, `N_A = 1`):** retain **DIRECT** for the measured 16³–64³ range. CONVOLUTION is slower in steady state at every size; REDUCED-DIRECT is close to DIRECT only at the largest size and remains experimental.
- **REDUCED-DIRECT:** benchmarked where eligible, retained as an opt-in experimental backend. The campaign does not justify a global promotion.
- **J+D skyrmions:** not included, as required; they remain DIRECT pending CPU-HAM-06.

The machine-readable measurements are in [`CPU_HAM_05_BACKEND_CROSSOVER.json`](CPU_HAM_05_BACKEND_CROSSOVER.json). Values below are fitted from complete production-process wall time, using `T = setup + steady_step * nstep`.

## Host, build, and method

| Item | Value |
|---|---|
| Machine | Apple M1, macOS Darwin 25.5, arm64, 8 physical/logical processors, 16 GiB |
| Compiler | `/opt/homebrew/bin/gfortran` 16.2.0 |
| Build | Release, double precision, CPU-only |
| FFT backend | FFTW enabled; provider FFTW, one FFT provider thread |
| MKL | Disabled and unavailable on this machine |
| OpenMP sweep | 1, 2, 4, 8 threads; `OMP_DYNAMIC=FALSE`, `OMP_PLACES=cores`, `OMP_PROC_BIND=close` |
| Samples | One complete run per `nstep`/backend/thread cell |
| Fit ladders | Main: 5/10/20 steps; short-range: 20/40/80; Nd convolution: 100/200/400 |
| Timed unit | Complete production process, including initialization and measurement phase |

No external MKL run was needed: sparse is portable and convolution was exercised through FFTW. A provider-specific MKL comparison would require another host and is outside this campaign.

The production executable does not export a pair-only timer or convolution stage timers. Therefore `pair_field_time_seconds` and exact convolution setup/transform/spectral-multiply timings are recorded as unavailable in the JSON. The authoritative steady-state quantity here is the complete production timestep. The driver also records an explicitly labelled `full_effective_field_time_estimate_seconds = steady_step / 2`, because the tested production integrators issue two effective-field calls per timestep; this is an estimate, not a replacement for an exported kernel timer.

## Eligibility and coverage

| Workload | DIRECT | REDUCED-DIRECT | SPARSE | CONVOLUTION | Measured sizes |
|---|---:|---:|---:|---:|---|
| bcc Fe | yes | no — non-scalar-J control is ineligible | yes | no — non-scalar-J control is ineligible | 13³, 20³, 32³ |
| dhcp Nd | yes | yes | yes | yes | 16³, 20³, 25³ |
| short-range scalar-J | yes | yes | yes | yes | 16³, 32³, 64³ |

`Ndirected` is `Natom × z`. The bcc Fe ladder has 4,394/16,000/65,536 atoms and 421,824/1,536,000/6,291,456 directed interactions. The Nd ladder has 16,384/32,000/62,500 atoms and 21,921,792/42,816,000/83,625,000 directed interactions. The short-range ladder has 4,096/32,768/262,144 atoms and 24,576/196,608/1,572,864 directed interactions. The JSON is the source of truth for all integer counts.

## Backend-BEST results

Setup is the fitted fixed cost in seconds. `Tsteady` is the fitted complete production timestep in seconds. Throughput is directed pair interactions per second, including the campaign ensemble factor where applicable.

### bcc Fe, medium range

| Size | Natom | z | Backend (best threads) | Setup s | Tsteady s | Spin-steps/s | M directed interactions/s |
|---|---:|---:|---|---:|---:|---:|---:|
| 13³ | 4,394 | 96 | DIRECT (8) | 0.1417 | 0.005722 | 174.8 | 1,474.4 |
| 13³ | 4,394 | 96 | SPARSE (4) | 0.2177 | 0.009787 | 102.2 | 862.0 |
| 20³ | 16,000 | 96 | DIRECT (4) | 0.4644 | 0.026471 | 37.8 | 1,160.5 |
| 20³ | 16,000 | 96 | SPARSE (4) | 1.0222 | 0.021954 | 45.6 | 1,399.3 |
| 32³ | 65,536 | 96 | DIRECT (8) | 1.9908 | 0.087556 | 11.4 | 1,437.1 |
| 32³ | 65,536 | 96 | SPARSE (8) | 2.1792 | 0.174921 | 5.7 | 719.3 |

The 20³ SPARSE slope is numerically lower than DIRECT, but its fit RMSE is 0.1032 s, larger than the fitted slope itself. With one sample per cell and no repeatable trend across the ladder, this is not sufficient evidence for a crossover.

### dhcp Nd, long range

| Size | Natom | z | Backend (best threads) | Setup s | Tsteady s | Spin-steps/s | M directed interactions/s |
|---|---:|---:|---|---:|---:|---:|---:|
| 16³ | 16,384 | 1,338 | CONVOLUTION (1) | 0.7973 | 0.001652 | 605.2 | 26,533.5 |
| 16³ | 16,384 | 1,338 | DIRECT (8) | 0.6514 | 0.013875 | 72.1 | 3,159.9 |
| 16³ | 16,384 | 1,338 | REDUCED-DIRECT (8) | 0.5943 | 0.016760 | 59.7 | 2,616.0 |
| 16³ | 16,384 | 1,338 | SPARSE (8) | 0.6696 | 0.039862 | 25.1 | 1,099.9 |
| 20³ | 32,000 | 1,338 | CONVOLUTION (4) | 0.8156 | 0.002664 | 375.4 | 32,146.8 |
| 20³ | 32,000 | 1,338 | DIRECT (4) | 0.9472 | 0.030388 | 32.9 | 2,818.0 |
| 20³ | 32,000 | 1,338 | REDUCED-DIRECT (8) | 0.7379 | 0.042283 | 23.7 | 2,025.2 |
| 20³ | 32,000 | 1,338 | SPARSE (8) | 0.9646 | 0.069134 | 14.5 | 1,238.6 |
| 25³ | 62,500 | 1,338 | CONVOLUTION (4) | 1.5778 | 0.005972 | 167.5 | 28,006.1 |
| 25³ | 62,500 | 1,338 | DIRECT (8) | 1.3678 | 0.051022 | 19.6 | 3,278.0 |
| 25³ | 62,500 | 1,338 | REDUCED-DIRECT (8) | 1.4108 | 0.065581 | 15.2 | 2,550.3 |
| 25³ | 62,500 | 1,338 | SPARSE (8) | 1.7722 | 0.133348 | 7.5 | 1,254.2 |

The 25³ non-convolution matrix is complete over 1/2/4/8 threads. Several lower-thread fits have high residuals; the backend-BEST rows use the measured eight-thread minima for DIRECT, REDUCED-DIRECT, and SPARSE. The convolution 25³ matrix is also complete over 1/2/4/8 threads.

For convolution, the FFT grid is the corresponding replicated periodic grid: 16³, 20³, or 25³, with the four-atom basis retained in the spectral multiply. FFT provider metadata is present in the JSON; internal transform-stage timings are not exported by the production executable.

### Short-range scalar-J control

| Size | Natom | z | Backend (best threads) | Setup s | Tsteady s | Spin-steps/s | M directed interactions/s |
|---|---:|---:|---|---:|---:|---:|---:|
| 16³ | 4,096 | 6 | DIRECT (1) | 0.0243 | 0.000222 | 4,513.2 | 221.8 |
| 16³ | 4,096 | 6 | REDUCED-DIRECT (1) | 0.0326 | 0.000317 | 3,157.1 | 155.2 |
| 16³ | 4,096 | 6 | SPARSE (1) | 0.0179 | 0.000401 | 2,494.1 | 122.6 |
| 16³ | 4,096 | 6 | CONVOLUTION (1) | 0.0164 | 0.000420 | 2,381.9 | 117.1 |
| 32³ | 32,768 | 6 | DIRECT (2) | 0.0832 | 0.001852 | 540.0 | 212.3 |
| 32³ | 32,768 | 6 | REDUCED-DIRECT (8) | 0.0619 | 0.002035 | 491.4 | 193.2 |
| 32³ | 32,768 | 6 | SPARSE (8) | 0.0628 | 0.002362 | 423.4 | 166.5 |
| 32³ | 32,768 | 6 | CONVOLUTION (2) | 0.0687 | 0.002984 | 335.1 | 131.8 |
| 64³ | 262,144 | 6 | REDUCED-DIRECT (8) | 0.4294 | 0.016069 | 62.2 | 195.8 |
| 64³ | 262,144 | 6 | DIRECT (8) | 0.4012 | 0.016808 | 59.5 | 187.2 |
| 64³ | 262,144 | 6 | SPARSE (4) | 0.4322 | 0.018747 | 53.3 | 167.8 |
| 64³ | 262,144 | 6 | CONVOLUTION (8) | 0.4449 | 0.025265 | 39.6 | 124.5 |

The 64³ REDUCED-DIRECT versus DIRECT difference is only about 4.4% in a one-sample run, so it is a tie/workload-dependent result, not a production-policy recommendation.

## Setup economics and measured crossovers

The break-even estimate is computed from the fitted models as `n = (setup_fast - setup_slow) / (Tsteady_slow - Tsteady_fast)` when the candidate has both a setup penalty and a steady-state advantage. A value of zero means the candidate is already no worse in fitted setup. Values are only reported for measured workload classes and are not extrapolated into an AUTO rule.

| Workload/size | Candidate | Compared with | Candidate best threads | Approx. n_break_even |
|---|---|---|---:|---:|
| Nd 16³ | CONVOLUTION | DIRECT | 1 vs 8 | 12 |
| Nd 16³ | CONVOLUTION | REDUCED-DIRECT | 1 vs 8 | 14 |
| Nd 16³ | CONVOLUTION | SPARSE | 1 vs 8 | 4 |
| Nd 20³ | CONVOLUTION | REDUCED-DIRECT | 4 vs 8 | 2 |
| Nd 20³ | CONVOLUTION | DIRECT/SPARSE | 4 vs 4/8 | 0 |
| Nd 25³ | CONVOLUTION | DIRECT | 4 vs 8 | 5 |
| Nd 25³ | CONVOLUTION | REDUCED-DIRECT | 4 vs 8 | 3 |
| Nd 25³ | CONVOLUTION | SPARSE | 4 vs 8 | 0 |

There is no robust measured SPARSE crossover for bcc Fe. In the low-z control, DIRECT (or the near-tied reduced variant at 64³) has the lower steady-state cost; convolution has no steady-state crossover in the measured range.

The workload boundary is therefore dominated by interaction density: long-range Nd has `z = 1338` and strongly favours convolution, while the `z = 6` control favours direct traversal. The bcc Fe `z = 96` results do not establish a sparse crossover on this machine.

The three-size ladders are sufficient for measured backend-BEST comparisons but not for a statistically defensible multi-parameter cross-size coefficient fit once each backend is allowed its own thread optimum. No cross-size model was forced; the report uses the supported per-size `nstep` fits and only measured-range comparisons.

## Locality/SFC conclusion

This campaign used the natural production ordering. CPU-HAM-02A remains the authoritative natural-versus-retained-locality comparison: the retained SFC/Morton experiment was workload-dependent, with no promotion justified for the Nd and J+D cases. The backend sweep does not overturn that result. Natural ordering remains the default; SFC remains opt-in and workload-dependent. No AUTO policy is implemented.

## Reproducibility and limitations

The campaign driver is [`cpu_ham05_driver.py`](../benchmarks/harness/cpu_ham05_driver.py). It copies the admitted templates to a temporary run directory, changes only the backend dispatch keywords plus deterministic zero-temperature run controls, verifies the backend activation marker, and fits complete production-process timings. Raw runs and summaries were kept outside the repository because benchmark outputs are ignored; the compact merged snapshot is committed as the JSON artifact above.

Important limitations are explicit rather than hidden in the conclusions:

- one sample per matrix cell;
- exact pair-only and convolution internal-stage timers are not exported by `sd.gfortran`;
- several lower-thread Nd 25³ non-convolution fits have high residuals despite the complete 1/2/4/8 thread matrix;
- hardware performance counters were unavailable on this host;
- no MKL provider comparison was performed or required;
- cross-size coefficient fits were not forced from three noisy size points;
- fits describe the measured ladders only and are not an AUTO decision model.
