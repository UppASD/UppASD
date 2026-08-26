# Benchmark analysis

Consumes records, produces understanding. Analysis code never produces
records and never re-runs a measurement.

Planned (WP-07):

* aggregation of raw samples into per-cell statistics (`sample_count`,
  `median`, `mad`, `minimum`, `maximum`) and `T(n) = T_fixed + n·t_step` fits;
* OpenMP scaling: `S_OMP(p) = T(N,1)/T(N,p)`, `E_OMP(p) = S_OMP(p)/p`,
  CPU-1T and CPU-BEST;
* GPU crossover against CPU-BEST, with crossings of 1.0×, 1.25×, 2.0× and
  5.0× classified as `below_tested_range`, `within_tested_range` or
  `above_tested_range`. Interpolation only between valid measured neighbouring
  points; never extrapolate and call the result a measured crossover;
* throughput metrics: spin-steps/s for all workloads,
  `directed_interaction_visits_per_second` for neighbour-list workloads, and a
  grid-normalized metric for FFT/dipole workloads. Metrics are not forced
  across workload classes;
* plots and reports.

Generated output goes to `analysis/out/`, which is gitignored.
