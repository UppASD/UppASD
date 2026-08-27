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

Planned (WP-07):

* aggregation of raw samples into per-cell statistics (`sample_count`,
  `median`, `mad`, `minimum`, `maximum`) and `T(n) = T_fixed + n·t_step` fits
  -- already implemented in `harness/aggregate.py` and `harness/steady_state.py`
  (WP-04/WP-03); WP-07 is expected to build the corresponding reports here;
* GPU crossover against CPU-BEST, with crossings of 1.0×, 1.25×, 2.0× and
  5.0× classified as `below_tested_range`, `within_tested_range` or
  `above_tested_range`. Interpolation only between valid measured neighbouring
  points; never extrapolate and call the result a measured crossover;
* throughput metrics: spin-steps/s for all workloads,
  `directed_interaction_visits_per_second` for neighbour-list workloads, and a
  grid-normalized metric for FFT/dipole workloads. Metrics are not forced
  across workload classes;
* further plots and reports beyond OpenMP scaling.

Generated output goes to `analysis/out/`, which is gitignored.
