# CPU-HAM-00 — Establish the production HamiltonianActions CPU baseline

**Model:** Luna
**Production optimization in this task:** none.

## Purpose

Before optimizing anything, determine exactly where CPU time goes and why OpenMP scaling degrades.

## Mandatory scope

Production CPU Hamiltonian implementation: `HamiltonianActions`.

Do **not** profile obsolete `heisge()` / `ApplyHamiltonian` as though they were production execution.

## Goal

Produce a call-path and hardware-counter baseline for at least:

1. dhcp Nd — very long-range `J`;
2. bcc Fe — medium-range `J`;
3. one short-range `J+D` model, preferably one benchmark skyrmion case.

Use representative production sizes, including the benchmark size that demonstrated poor Nd OpenMP scaling.

## A. Trace production call path

Starting from the actual LLG timestep, establish:
- which `effective_field` entry point is used;
- which `HamiltonianActions` routines dominate;
- number of effective-field calls per timestep;
- energy requested/not requested;
- relevant onsite versus pair contributions;
- any full-array clears;
- OpenMP region boundaries;
- fork/join/barrier count.

Produce a concise call graph.

## B. Profile term contributions

For each case determine approximate time split between:
- isotropic exchange;
- DMI;
- onsite terms;
- field finalization;
- energy accumulation/reduction;
- non-Hamiltonian integrator work.

Do not add intrusive synchronization merely to generate timings if a profiler can obtain them externally.

## C. Hardware counters

Where tools are available, collect as a function of OpenMP threads:

`1, 2, 4, 8, ...` up to available physical cores.

At minimum seek:
- memory bandwidth;
- IPC;
- LLC misses;
- cycles/stalls attributable to memory where available;
- vectorization status;
- CPU frequency;
- NUMA/affinity.

Acceptable tools include `perf`, LIKWID, VTune, or equivalent.

Do not add a hard dependency on any profiler.

If a metric is unavailable, record it as unavailable.

## D. OpenMP placement

Record:
- `OMP_NUM_THREADS`;
- `OMP_PLACES`;
- `OMP_PROC_BIND`;
- CPU affinity;
- NUMA topology.

Compare the current binding against one sensible alternative if supported.

## E. Test the memory-bound hypothesis

For each workload answer:
- Does memory bandwidth saturate as thread count increases?
- Does the point of saturation correlate with OpenMP speedup flattening?

Particularly compare Nd vs bcc Fe vs short-range `J+D`.

## F. Energy execution audit

Establish exactly:
- whether `effective_field` always computes energy;
- whether ordinary LLG predictor/corrector callers consume it;
- where that energy is currently used;
- whether energy reduction introduces measurable OpenMP overhead.

Do not change semantics yet.

## G. Full-array work audit

Identify all unconditional `O(N)` operations around `HamiltonianActions`, including array zeroing/copying, and classify whether each is semantically required.

## H. Deliverable

Create `docs/CPU_HAM_00_BASELINE_PROFILE.md`.

Include:
- call graph;
- workload table;
- OMP scaling;
- counter results;
- energy-use findings;
- full-array-pass findings;
- bottleneck classification.

Classify each case as primarily:
- `MEMORY_BANDWIDTH`
- `CACHE_LATENCY/GATHER`
- `COMPUTE`
- `OPENMP_OVERHEAD`
- `LOAD_IMBALANCE`
- `MIXED/UNRESOLVED`

Do not force a classification if evidence is insufficient.

## Checklist

- [x] Production call path proven.
- [x] Obsolete heisge excluded.
- [x] Nd profiled.
- [x] bcc Fe profiled.
- [x] Short-range `J+D` profiled.
- [x] Thread scaling measured.
- [x] Memory bandwidth investigated.
- [x] IPC investigated.
- [x] LLC/cache behaviour investigated.
- [x] Vectorization status investigated.
- [x] OpenMP placement recorded.
- [x] NUMA considerations recorded.
- [x] Energy execution semantics mapped.
- [x] Full-array passes inventoried.
- [x] Bottleneck classification documented.
- [x] No optimization committed.

## Commit

`CPU-HAM-00: profile production CPU Hamiltonian`
