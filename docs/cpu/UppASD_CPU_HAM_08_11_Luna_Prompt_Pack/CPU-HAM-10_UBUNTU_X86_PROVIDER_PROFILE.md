# CPU-HAM-10 — Ubuntu/x86 Provider and Hardware-Counter Campaign

**Model:** Luna

## Dependency

CPU-HAM-09 complete.

## Target environment

Ubuntu Linux on x86-64.

Preferred available tooling:
- Intel VTune Profiler;
- `perf`;
- LIKWID where available;
- oneMKL;
- FFTW;
- gfortran and/or Intel oneAPI compiler.

Do not require every tool simultaneously.

## Purpose

Answer the questions that could not be resolved on the Apple M1 campaign:

1. Is long-range DIRECT actually memory/cache/gather limited?
2. How does DIRECT scale with OpenMP on Linux/x86?
3. Does a true oneMKL inspector/executor sparse SpMM backend outperform the portable CSR loop?
4. How do FFTW and oneMKL FFT providers compare for CPU convolution?
5. What are the real backend crossovers on this Ubuntu machine?
6. How should thread ownership be configured?

This is an evidence campaign plus narrowly scoped provider implementation where necessary.

## A. Freeze environment and builds

Record:
- exact git commit;
- clean/dirty state;
- CPU model;
- sockets;
- physical cores;
- SMT;
- NUMA topology;
- RAM;
- compiler versions;
- oneMKL version;
- FFTW version;
- VTune version;
- relevant CMake options;
- optimization flags.

Create at least:
- GNU/FFTW build;
- GNU or Intel + oneMKL build where supported.

Do not compare binaries with materially different optimization levels without labelling them.

## B. CPU affinity and NUMA

Use explicit:
- `OMP_NUM_THREADS`;
- `OMP_PLACES=cores`;
- `OMP_PROC_BIND=close/spread` as tested;
- process affinity.

For multi-socket machines, separately characterize:
- one-socket runs;
- cross-socket/full-node runs.

Do not mix them silently.

## C. DIRECT VTune/perf characterization

Required workloads:
- dhcp Nd long-range;
- bcc Fe medium-range;
- short-range scalar-J or J+D control.

Thread sweep:
`1, 2, 4, 8, ...` physical cores.

Collect where available:
- memory bandwidth;
- L1/L2/LLC miss behavior;
- DRAM bandwidth;
- IPC/CPI;
- backend-bound/front-end-bound breakdown;
- memory-bound fraction;
- cache-bound fraction;
- branch behavior if relevant;
- vectorization/SIMD utilization;
- TLB metrics;
- thread imbalance;
- CPU frequency.

For VTune, useful analyses may include:
- Hotspots;
- Memory Access;
- HPC Performance Characterization.

Do not blindly run every profiler mode if one already answers the question.

## D. Bottleneck verdict

For each workload classify DIRECT as:
- DRAM bandwidth limited;
- cache/gather latency limited;
- compute limited;
- load-imbalance limited;
- synchronization/fork-join limited;
- mixed.

Support the verdict with counters.

Do not use the old working hypothesis as evidence.

## E. oneMKL sparse provider

The portable sparse backend remains the correctness reference for sparse representation.

Implement or revive an **optional** oneMKL sparse provider if current build infrastructure permits it cleanly.

Requirements:
- optional build capability;
- no mandatory oneMKL dependency;
- persistent sparse handle;
- persistent inspector/optimization state;
- scalar J initially;
- dense three-component RHS;
- prefer SpMM or equivalent inspector/executor apply;
- no structure rebuild per timestep.

If oneMKL API integration requires broad build-system surgery:
STOP and report.

Do not resurrect obsolete legacy MKL sparse APIs blindly.

## F. oneMKL sparse correctness

Compare:
- oneMKL sparse vs portable sparse;
- oneMKL sparse vs canonical DIRECT.

Use:
- Nd;
- Fe;
- random states;
- multi-basis if HAM-08 validated it;
- multiple ensembles.

Energy must remain field-derived.

## G. Sparse threading ownership

Benchmark:
- UppASD outer OpenMP threading;
- oneMKL internal threading.

Avoid nested oversubscription.

Preferred model is usually:
- one backend call;
- oneMKL owns threads.

But measure rather than assume.

## H. FFTW versus oneMKL convolution

If the HAM-09 convolution abstraction permits both providers, implement/enable optional oneMKL DFTI provider.

Requirements:
- same reduced-stencil spectral kernel;
- same normalization;
- same physics;
- persistent descriptors;
- same field/energy oracle.

Compare providers on:
- Nd;
- eligible Fe;
- scalar J;
- J+D if provider path naturally supports it.

## I. FFT provider threading

For each provider measure:
- 1;
- 2;
- 4;
- 8;
- etc. physical cores.

Record best provider thread count by size.

Do not assume FFTW and oneMKL have identical best threading.

## J. Production crossover matrix

Using repeated measurements and the general benchmark harness, compare:

- DIRECT;
- REDUCED-DIRECT only as an explicit experimental control;
- portable SPARSE;
- oneMKL SPARSE if implemented;
- FFTW CONVOLUTION;
- oneMKL CONVOLUTION if implemented.

Required workloads:
- Nd long-range;
- Fe medium-range;
- short-range control.

Use multiple system sizes.

Use enough repetitions to distinguish ~5-10% effects.

Report medians and dispersion.

## K. Setup amortization

For oneMKL sparse and both convolution providers report:
- setup time;
- steady-state time;
- break-even step count versus DIRECT.

## L. Final decisions

At task completion answer explicitly:

1. Is DIRECT bandwidth/gather limited for Nd?
2. Does optimized sparse library execution add value?
3. Which sparse implementation, if any, should remain production-supported?
4. FFTW or oneMKL: which provider is better on this machine?
5. Is one provider robust enough to recommend, or should both remain?
6. What are the backend crossover regions?
7. Does any evidence now justify AUTO?
   Do not implement AUTO in this task.

## M. Deliverable

Create:

`docs/CPU_HAM_10_UBUNTU_X86_PROFILE.md`

Include:
- machine/build provenance;
- VTune/perf findings;
- sparse provider findings;
- FFT provider findings;
- crossover tables;
- thread ownership recommendations;
- limitations.

## Checklist

- [x] Ubuntu/x86 environment frozen; provenance is in `docs/CPU_HAM_10_UBUNTU_X86_PROFILE.md`.
- [x] Affinity controlled.
- [x] NUMA topology recorded.
- [ ] DIRECT VTune/perf profile completed — `perf_event_paranoid=4` blocks counters; VTune/LIKWID are unavailable.
- [ ] Nd bottleneck proven — timing suggests memory/gather sensitivity, but no hardware-counter proof is possible on this host.
- [ ] Fe bottleneck proven — same counter limitation.
- [x] Short-range control profiled with the production harness and HAM-09 stage benchmark.
- [x] oneMKL sparse provider evaluated or explicit STOP rationale recorded — modern APIs were inventoried, but integration is not cleanly available through the current provider/build boundary.
- [ ] oneMKL sparse parity passes if implemented — not implemented; no parity claim made.
- [ ] Sparse provider thread ownership measured — oneMKL provider not implemented.
- [x] FFTW convolution measured.
- [x] oneMKL convolution measured if cleanly implementable — not cleanly implementable in the current CPUFFTProvider/CMake boundary; explicit STOP rationale is recorded.
- [ ] FFT provider parity passes — only FFTW is wired, so cross-provider parity is not applicable yet.
- [x] FFT provider threading sweep complete for the available FFTW provider (1/2/4/8 threads).
- [ ] Repeated crossover matrix complete — representative two-size production checks are complete; REDUCED-DIRECT remains an experimental control and the full provider matrix awaits oneMKL implementations.
- [x] Setup amortization reported.
- [x] Medians/dispersion reported.
- [x] No unsupported AUTO policy introduced.
- [x] Final provider/backend recommendations documented.

## Commit

`CPU-HAM-10: characterize CPU Hamiltonian on Ubuntu x86`
