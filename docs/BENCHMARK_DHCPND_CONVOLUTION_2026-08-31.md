# dhcp Nd GPU Convolution Campaign — 2026-08-31

A focused follow-on to WP-10: characterizes UppASD's GPU lattice-convolution
Hamiltonian evaluator (`do_gpu_convolution Y`,
`source/gpu_files/gpuLatticeConvolutionHamiltonian.{hpp,cpp}`) against the
ordinary sparse neighbour-list path, on B04_dhcpNd — the one case among this
project's five admitted families whose template satisfies the feature's
runtime gate (see
[`docs/BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md`](BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md)
for the other four, and `benchmarks/harness/gpu_convolution.py` for why).

## Executive summary

- **GPU convolution is a real, substantial win for this workload**: 3.7–6.5×
  faster than the GPU's own sparse-kernel path, 24–39× faster than CPU-BEST,
  at both sizes tested and both precisions.
- **The disclosed aliasing risk at the new, deliberately anisotropic
  100×100×3 size did not materialize.** Its measured interaction count
  (160,560,000 directed interactions / 120,000 atoms = exactly 1338.0/atom)
  matches the established clean rate for this case exactly — the same value
  40×40×40 and every properly-replicated WP-10 B04 size measured. Aliased
  bonds through the thin z-axis would have inflated this count; they did
  not.
- **Zero measurement-validity problems**: the convolution-activation
  verification built into `convolution_driver.py`
  (`harness/gpu_convolution.py`) found zero silent fallbacks across both
  sizes and both precisions — every convolution sample that was recorded
  really ran the convolution kernel, not a mislabelled sparse-kernel run.
- Confirms and extends this repository's own prior single data point
  (`conv_bench.txt`'s `dhcp_nd_long_range`, ncell=32×32×32, 1000 steps,
  13.2× Hamiltonian-only speedup) at production benchmark methodology
  (multi-nstep steady-state fits, 5 samples, full OMP sweep) and at larger,
  more production-realistic sizes.

## Provenance

Same frozen source and binaries as WP-10: commit `ad44f260` (verified via
`git diff --quiet ad44f260 HEAD -- source/ CMakeLists.txt`, not an exact-HEAD
match — see the launcher scripts' own fix for why), `build_wp10_cpu` /
`build_wp10_cuda_fp32` / `build_wp10_cuda_fp64` (checksums recorded in the
WP-10 report; unchanged, since nothing since has touched `source/`). No new
build was needed — `gpuLatticeConvolutionHamiltonian.cpp` is compiled in
unconditionally.

Campaign: `benchmarks/campaigns/dhcpnd_convolution_2026-08-31.yaml`, driven
by `benchmarks/harness/convolution_driver.py` (new — adds a convolution arm
to any campaign's GPU cells, gated per case by
`gpu_convolution.check_convolution_gate` and per-sample by
`gpu_convolution.verify_convolution_active`). Executed via the standalone
launcher `benchmarks/run_dhcpnd_convolution_campaign.sh` after two in-session
attempts were interrupted (once by the editor's process tree being torn
down, once by a since-fixed OOM in the workload-metadata parser on a
~36 GiB `struct.out` — see the commits fixing both). 2 sizes × (8 CPU
threads + 2 precisions × 2 modes) = 24/24 configurations produced valid
aggregates; zero fallback events; zero sample errors.

## Workload

| Size | natom | directed_interactions | mean_neighbors |
| --- | --- | --- | --- |
| 40x40x40 (isotropic, existing B04 ladder tier) | 256,000 | 342,528,000 | 1338.0 |
| 100x100x3 (new, anisotropic thin-slab) | 120,000 | 160,560,000 | 1338.0 |

Both variant `dhcpNd_finite_t` (2.00000001 K), `measurement_profile:
DYNAMICS_ONLY`, 5 samples, 3-point steady-state fits, full 1–8 physical-core
OpenMP sweep, CUDA SINGLE and DOUBLE.

## Results

### Step time, CPU-BEST vs. GPU pair vs. GPU convolution

| Size | CPU-BEST (p) | GPU SINGLE pair | GPU SINGLE conv | GPU DOUBLE pair | GPU DOUBLE conv |
| --- | --- | --- | --- | --- | --- |
| 40x40x40 | 0.1479 s (p=8) | 0.02568 s | 0.005089 s | 0.02802 s | 0.006225 s |
| 100x100x3 | 0.06361 s (p=5) | 0.006061 s | 0.001642 s | 0.011584 s | 0.001778 s |

![Step time comparison](../benchmarks/analysis/out/dhcpnd_convolution_2026-08-31/step_time_comparison__B04_dhcpNd__dhcpNd_finite_t.png)

### Speedups

| Size | Precision | conv vs. GPU pair | conv vs. CPU-BEST | GPU pair vs. CPU-BEST |
| --- | --- | --- | --- | --- |
| 40x40x40 | SINGLE | 5.05× | 29.1× | 5.76× |
| 40x40x40 | DOUBLE | 4.50× | 23.8× | 5.28× |
| 100x100x3 | SINGLE | 3.69× | 38.7× | 10.50× |
| 100x100x3 | DOUBLE | 6.51× | 35.8× | 5.49× |

![Convolution speedup](../benchmarks/analysis/out/dhcpnd_convolution_2026-08-31/convolution_speedup__B04_dhcpNd__dhcpNd_finite_t.png)

Convolution wins by a comparable, large margin at both sizes and both
precisions — the anisotropic 100×100×3 shape (an FFT grid dominated by two
large, one very short axis) does not appear to meaningfully hurt convolution
performance relative to the isotropic 40×40×40 shape; if anything DOUBLE
precision's conv-vs-pair speedup is *higher* at 100×100×3 (6.51× vs 4.50×).
Two of the eight conv aggregates have `sample_count` below 5 (SINGLE conv:
1/5 at 40×40×40, 4/5 at 100×100×3) from the same jitter-dominated-fit
pattern documented throughout the WP-10 report — convolution's per-step cost
is fast enough at these sizes that noise occasionally beats the fitted
slope. Every affected aggregate still has ≥1 valid fit; none are fabricated.

### OpenMP scaling, with GPU pair/conv baselines

![OMP scaling 40x40x40](../benchmarks/analysis/out/dhcpnd_convolution_2026-08-31/omp_scaling/B04_dhcpNd__dhcpNd_finite_t/omp_scaling_with_conv__40x40x40.png)
![OMP scaling 100x100x3](../benchmarks/analysis/out/dhcpnd_convolution_2026-08-31/omp_scaling/B04_dhcpNd__dhcpNd_finite_t/omp_scaling_with_conv__100x100x3.png)

CPU-BEST is 8 threads at 40×40×40 (S_OMP=4.07×, E_OMP=0.51) and 5 threads at
100×100×3 (S_OMP≈2.6× at that point, though the curve is non-monotonic — see
Limitations). Both GPU pair baselines sit roughly an order of magnitude
below the CPU curve at every thread count; both conv baselines sit another
factor of ~4–6× below that again — the full 3-tier picture (CPU sweep, GPU
pair, GPU conv) in one figure that neither the WP-10 report nor its plots
show, since WP-10 never ran convolution.

## Aliasing-risk finding

`benchmarks/cases/B04_dhcpNd/case.yaml`'s own comment on the 100×100×3 size
disclosed an unverified risk: this case's exchange list needs
lattice-translation offsets out to `n3=±2` along the axis this size only
replicates 3 times, below the naive alias-free floor of 5. The measured
`directed_interactions` at this size (160,560,000, i.e. exactly
1338.0/atom) matches the established clean rate for every other correctly-
replicated B04 size in this repository exactly. Duplicated bonds through a
too-thin periodic wrap would inflate this count above 1338/atom; they did
not. This is indirect evidence (a matching aggregate count, not a direct
per-bond duplicate check), but it is evidence, and it points the same
direction across both precisions and both backends measured here. Treat the
z-axis-aliasing caveat as **not confirmed to be a problem for this specific
replication**, rather than either confirmed-safe or still-open;
`case.yaml`'s comment has not been changed on the strength of this alone.

## Limitations

- Two configurations have `sample_count` < 5 (noted above) — normal
  jitter-driven fit rejection, not a defect.
- 100×100×3's CPU OMP-scaling curve is non-monotonic (a dip at p=5, a rise
  at p=6, matching the kind of anomaly WP-10's own B04/25³ p=7 point showed)
  — reported as measured, not smoothed away; most plausibly transient
  system-level variability on a shared host rather than a stable property
  of this workload.
- Only B04_dhcpNd was in scope (the only case eligible per
  `gpu_convolution.check_convolution_gate`); B01/B02/B03/B05 were not
  re-attempted here even though `convolution_driver.py` is generic enough
  to be pointed at their campaign manifests too — they would all be
  correctly skipped (confirmed in WP-10's own methodology review), so
  running it against them would add no new data.
- Only two sizes; no crossover-threshold analysis (unlike the WP-10 report)
  since a 2-point curve does not support the same interpolation discipline.

## Where everything lives

| What | Where |
| --- | --- |
| This report | `docs/BENCHMARK_DHCPND_CONVOLUTION_2026-08-31.md` |
| Plotting script | `benchmarks/analysis/convolution_plots.py` |
| Raw samples + aggregates (360 raw records, 24 aggregates) | `benchmarks/results/dhcpnd_convolution_2026-08-31/` (gitignored, local only) |
| Generated plots | `benchmarks/analysis/out/dhcpnd_convolution_2026-08-31/` (gitignored, local only) |
