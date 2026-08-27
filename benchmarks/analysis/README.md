# Benchmark analysis

Consumes records, produces understanding. Analysis code never produces
records and never re-runs a measurement.

Present (WP-05):

* `omp_scaling_report.py` — OpenMP scaling reports (blueprint section H):
  `build_scaling_summary` turns one or more problem sizes' OpenMP-scaling
  `aggregate` families (via `harness.omp_scaling`) into a machine-readable
  summary (`write_scaling_summary`, plain JSON) naming each size's
  speedup/efficiency table, CPU-1T/CPU-BEST and `p_best` -- plus one
  `p_best_vs_natom` table across every size, sorted by `natom` regardless of
  input order. `plot_omp_scaling` renders step-time/speedup/efficiency-vs-
  threads (one figure per size) and `p_best` versus `N_atom` as PNGs.

Present (WP-07):

* `crossover.py` -- `S_GPU_BESTCPU`/`S_GPU_1T` (`compute_gpu_speedup`, the
  ratio of a matched CPU/GPU `steady_step_seconds` aggregate pair's medians)
  and crossover-threshold detection (`find_crossover`/`find_all_crossovers`
  over a `build_speedup_curve`-sorted curve, default thresholds 1.0x/1.25x/
  2.0x/5.0x). Interpolation happens only between two measured neighbouring
  points that bracket the threshold, in log(natom) vs log(speedup) space (a
  documented modelling choice, not a measured fact); a threshold already
  exceeded at the smallest tested size, or never reached at the largest, is
  classified `below_tested_range`/`above_tested_range` with no fabricated
  `crossover_natom` -- never extrapolated and presented as measured. An
  interpolated `crossover_natom` is rounded to 2 significant figures
  (`round_significant`) so it is never displayed with spurious precision.
* `throughput.py` -- `spin_steps_per_second` (every workload),
  `directed_interaction_visits_per_second` (neighbour-list workloads only)
  and `fft_grid_points_per_second` (FFT/dipole workloads only); each raises
  rather than applying itself to the wrong `workload_class`.
  `compute_throughput` dispatches on `workload_class` and returns exactly the
  metrics that apply.
* `campaign_summary.py` -- section G's machine-readable summary:
  `build_case_summary` combines one case/variant's measured (size,
  precision) GPU cells (caller-supplied CPU-BEST/CPU-1T/GPU
  `steady_step_seconds` aggregates -- this module selects none of them
  itself) into per-precision speedup/throughput rows, crossover results,
  measured `natom` range and quality status, plus `R_GPU_32_64` precision
  penalties (`harness.precision_metrics`) at every size measured in both
  SINGLE and DOUBLE. `write_summary_json`/`write_summary_csv`/
  `write_crossover_csv`/`write_precision_penalty_csv` render it to disk.
* `scaling_plots.py` -- plots 1, 2 and 5-9 of section E (step time vs
  `N_atom`/`N_directed`, GPU/BESTCPU speedup vs `N_atom`/`N_directed`, SINGLE
  vs DOUBLE GPU performance, crossover summary, FFT/dipole time vs grid
  size), each returning `None` when the summary has nothing relevant to plot
  rather than an empty or misleading figure. Plots 3-4 (OpenMP scaling and
  efficiency) are `omp_scaling_report.plot_omp_scaling`'s job (WP-05) and are
  not reproduced here. Medians are drawn with MAD error bars wherever the
  summary carries one.
* `markdown_report.py` -- `generate_markdown_report` renders one case/
  variant's summary as Markdown with three explicitly separated sections:
  **FACT** (measured medians and ratios only), **INTERPOLATION** (crossover
  thresholds located between measured neighbours, or their
  below/above-tested-range classification), and **HYPOTHESIS** (hedged
  interpretive sentences, e.g. a large `R_GPU_32_64` flagged as *consistent
  with* a bandwidth-limited DOUBLE path) -- interpretations are never
  restated as measured or interpolated claims.

Generated output goes to `analysis/out/`, which is gitignored.
