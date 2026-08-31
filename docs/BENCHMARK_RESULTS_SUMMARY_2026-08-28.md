# WP-10 Benchmark Results — Summary

A quick-look companion to the full report,
[`BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md`](BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md)
(all 14 required sections, full provenance, full methodology review — read
that for anything not covered here). This page is the "what happened and
what does it look like" digest.

## What was run

5 production workload families (bcc Fe, 2D skyrmion, 3D skyrmion, dhcp Nd,
dipole/FFT), 3 sizes each, on one frozen commit
(`ad44f260aa918aca2304d2bf22bd544e11e199bc`): a full 1–8 OpenMP thread
sweep, CUDA SINGLE and DOUBLE, 5 samples/cell, 3-point steady-state fits.
180/180 planned configurations produced valid data (two bugs cost the
first pass 100% of B02/B03's GPU data and 100% of B05's data; both were
found, fixed, and the affected cells re-measured — see the full report's
§15 for details). Total: 2,701 production-executable runs, ~11.4 hours of
actual measurement.

## Regenerating the plots

```sh
python3 benchmarks/analysis/wp10_plots.py \
  --results-dir benchmarks/results/wp10_2026-08-28 \
  --out-dir benchmarks/analysis/out/wp10_2026-08-28
```

Requires the raw/aggregate records under `benchmarks/results/`, which are
**not** committed (`benchmarks/.gitignore`: raw and aggregated result
records are always generated, gitignored artifacts — see
`benchmarks/README.md`'s "Result-data policy"). The PNGs below are
likewise generated, not tracked; this page links to them by relative path
for local viewing after running the command above, and the images shown
here were captured from that run on this host.

48 PNGs across two directories:

- `benchmarks/analysis/out/wp10_2026-08-28/omp_scaling/<case>__<variant>/` —
  one 3-panel (step time / speedup / efficiency vs. thread count) figure
  per size, plus one CPU-BEST-thread-count-vs-N_atom figure, per
  case/variant (24 PNGs).
- `benchmarks/analysis/out/wp10_2026-08-28/crossover/` — step time,
  GPU/CPU-BEST speedup, SINGLE-vs-DOUBLE, crossover-threshold and (B05
  only) FFT-grid-size figures, one set per case/variant (24 PNGs plus 11
  case/variant combinations that legitimately have nothing to plot for a
  given figure — e.g. `step_time_vs_ndirected` on B05, which has no
  neighbour list — are skipped, not blank).

## Headline result: GPU already wins everywhere tested

Every family, every precision, beats CPU-BEST already at the *smallest*
size tested (15,625–16,384 atoms). None of these curves ever cross back
below 1×, so there is nothing resembling a classic "small problems favor
CPU" crossover within this campaign's size range — see the full report's
§7 for why the true 1× point is out of reach of this campaign's ladders,
not absent.

![B01 bcc Fe speedup vs N_atom](../benchmarks/analysis/out/wp10_2026-08-28/crossover/speedup_vs_natom__B01_bccFe__bcc_fe_t0.png)

*B01_bccFe (T=0): GPU/CPU-BEST speedup vs. atom count, both precisions.
Dotted/dashed lines mark the 1×/1.25×/2×/5× thresholds — every measured
point already clears 5× at the smallest tested size.*

## The biggest number in the whole campaign: B05's algorithmic crossover

B05's CPU arm uses the only validated CPU dipole mode, a direct `O(N²)`
sum; the GPU arm uses an FFT convolution. This is not primarily a
"faster chip" story — it's two different algorithmic complexity classes
racing each other, and it shows:

![B05 dipole/FFT speedup vs N_atom](../benchmarks/analysis/out/wp10_2026-08-28/crossover/speedup_vs_natom__B05_dipoleFFT__dipole_on.png)

*B05_dipoleFFT: GPU/CPU-BEST speedup reaches ~1,039× (DOUBLE) at the
largest tested size (16,384 atoms) — by a wide margin the largest
crossover measured anywhere in this campaign. See the full report §10 for
why (CPU cost grows ~quadratically with atom count; GPU cost stays flat,
dominated by the FFT rather than the atom count directly).*

## OpenMP scaling: efficiency drops off past 3–4 threads (mostly)

![B01 bcc Fe OMP scaling at 20x20x20](../benchmarks/analysis/out/wp10_2026-08-28/omp_scaling/B01_bccFe__bcc_fe_t0/omp_scaling__20x20x20.png)

*B01_bccFe (T=0, 16,000 atoms): step time, speedup and efficiency vs.
thread count. CPU-BEST here is 8 threads at only 3.0× speedup (E≈0.38) —
typical of most cells in this campaign. B04_dhcpNd is the one family that
scales meaningfully better (E_OMP 0.6–0.9 through 5–8 threads on most
sizes) — see the full report §5/§9 for why its explicit, pre-flattened
long-range bond list appears to help both single-thread cost *and*
OpenMP scaling.*

## Precision: SINGLE is a real, modest, matched-vs-unmatched win

- CUDA DOUBLE vs. CPU DOUBLE is a genuinely precision-matched comparison
  (both operate at the same numerical width).
- CUDA SINGLE vs. CPU DOUBLE is a real, supported production
  configuration but **not** precision-matched (CPU never runs SINGLE —
  `UPPASD_PRECISION` doesn't gate `dblprec`).
- Median R_GPU_32_64 (DOUBLE time / SINGLE time) across 17 of 18 cells is
  **2.2×** — see the full report §6/§8 for the one flagged low-sample-count
  outlier (B03/50³, `sample_count=3`, not treated as a real precision
  inversion).

![B01 bcc Fe SINGLE vs DOUBLE](../benchmarks/analysis/out/wp10_2026-08-28/crossover/single_vs_double__B01_bccFe__bcc_fe_t0.png)

## Interaction count vs. atom count: partial collapse, not universal

![B01 bcc Fe step time vs N_directed](../benchmarks/analysis/out/wp10_2026-08-28/crossover/step_time_vs_ndirected__B01_bccFe__bcc_fe_t0.png)

Normalizing by directed-interaction count collapses B01/B02/B03 (all
symmetry-expanded neighbour lists) onto a consistent ~9–13 ns/interaction,
a real ~20× tighter spread than raw atom count gives across those three
families — but it does *not* pull B04 onto that same curve: B04's
pre-flattened explicit bond list is ~27× more interaction-efficient at the
same interaction count. List *representation*, not just size, matters —
full evidence and numbers in the full report's §9.

## Known display quirk

`crossover_summary__*.png` plots (the threshold-marker figures) can look
visually crowded on families where most/all thresholds land at the same
"below tested range" edge — several annotation labels stack on the same
point. This is existing WP-07 plotting code (`analysis/scaling_plots.py`,
unmodified here), informationally correct but cosmetically dense; not
something this task changed.

## Where everything lives

| What | Where |
| --- | --- |
| Full 14-section report | `docs/BENCHMARK_PRODUCTION_CHARACTERIZATION_2026-08-28.md` |
| This summary | `docs/BENCHMARK_RESULTS_SUMMARY_2026-08-28.md` |
| Plotting script | `benchmarks/analysis/wp10_plots.py` |
| Machine-readable summary (OMP tables, crossover, precision penalty, throughput — everything the plots are drawn from) | `benchmarks/results/wp10_2026-08-28/analysis/summary.json` (regenerate with `benchmarks/analysis/wp10_summary.py`) |
| Raw samples + aggregates (2,701 + 360 records) | `benchmarks/results/wp10_2026-08-28/` (gitignored, local only) |
| Generated plots | `benchmarks/analysis/out/wp10_2026-08-28/` (gitignored, local only) |
