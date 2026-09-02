# CPU-HAM-05 — Benchmark CPU pair-Hamiltonian backend crossovers

**Model:** Luna

## Dependencies

`CPU-HAM-01B`, `CPU-HAM-02C`, `CPU-HAM-04B`, plus `CPU-HAM-03B` if sparse GO.

## Purpose

Determine which CPU pair backend is appropriate for which workload.

No new production optimization in this task.

## Backends

At minimum:
- DIRECT;
- REDUCED-DIRECT where eligible;
- CONVOLUTION where eligible.

Add SPARSE if available.

## A. Workloads

Required:
1. bcc Fe — medium-range;
2. dhcp Nd — long-range;
3. scalar-J short-range periodic control if needed for low-z regime.

`J+D` skyrmions remain DIRECT until CPU-HAM-06 unless a rigorous mixed comparison is possible.

## B. Sizes

Use general benchmark size ladders.

## C. Threads

Sweep physical OpenMP/FFT thread counts.

Each backend may have a different optimum.

Compare backend-BEST, not same arbitrary `p`.

## D. Metrics

Report:
- `Natom`;
- `Ndirected`;
- mean neighbours;
- backend;
- thread count;
- setup cost;
- pair-field time;
- full effective-field time;
- spin-steps/s;
- interaction throughput where meaningful.

For convolution also report FFT grid, FFT setup, transforms, and spectral multiply.

## E. Models

Fit only if supported:

DIRECT:
\[
T_D\sim a_DNz+b_DN+c_D
\]

SPARSE:
\[
T_S\sim a_SNz+b_SN+c_S
\]

CONVOLUTION:
\[
T_C\sim a_FN\log N+a_MN\cdot\text{basis factor}+c_F
\]

Do not force fits to create AUTO.

## F. Crossover

Determine measured/interpolated crossover regions.

Report `N`, `z/Ndirected`, `N_A`, thread count, and machine.

## G. Locality result

Include natural DIRECT versus retained locality variants.

Determine whether SFC remains useful after reduced-stencil work.

## H. Setup economics

Estimate `n_break_even` required to amortize SPARSE/CONVOLUTION setup.

## I. Output

Create `docs/CPU_HAM_05_BACKEND_CROSSOVER.md` and machine-readable data.

## J. Decision

Classify preferred backend for each workload class empirically.

Do not create AUTO yet.

## Checklist

**Result:** complete for the measured ranges. The detailed report is [`docs/CPU_HAM_05_BACKEND_CROSSOVER.md`](../../CPU_HAM_05_BACKEND_CROSSOVER.md), with machine-readable measurements in [`docs/CPU_HAM_05_BACKEND_CROSSOVER.json`](../../CPU_HAM_05_BACKEND_CROSSOVER.json). CPU-HAM-05 added no production backend or dispatch policy.

- [x] Implementations frozen.
- [x] DIRECT benchmarked.
- [x] REDUCED-DIRECT benchmarked where eligible; bcc Fe is explicitly ineligible.
- [x] SPARSE benchmarked if available.
- [x] CONVOLUTION benchmarked where eligible.
- [x] Nd matrix complete.
- [x] bcc Fe matrix complete.
- [x] Short-range low-z control complete.
- [x] Backend-optimal threads found.
- [x] Setup cost included.
- [x] Steady-state cost included.
- [x] Setup amortization estimated.
- [x] Crossovers identified only within measured range.
- [x] SFC/locality conclusion updated; natural ordering remains default and SFC remains workload-dependent/opt-in.
- [x] No AUTO policy implemented.

## Commit

`CPU-HAM-05: characterize CPU pair backend crossovers`
