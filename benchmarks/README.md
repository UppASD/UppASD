# UppASD production benchmark subsystem

This directory holds the benchmark framework described in
[docs/UppASD Production Benchmark Framework.md](../docs/UppASD%20Production%20Benchmark%20Framework.md),
the master blueprint that governs it. Where this README and the blueprint
disagree, the blueprint wins.

## Scope

The framework exists to produce reproducible, defensible answers about
production UppASD performance:

* CPU/OpenMP scaling and the best thread configuration per workload;
* CPU ↔ CUDA crossover as a function of problem size;
* SINGLE versus DOUBLE precision behaviour;
* interaction-range and connectivity dependence;
* neighbour-list versus FFT/dipole workloads;
* setup cost versus steady-state simulation cost;
* measurement overhead;
* machine-to-machine differences.

In scope as execution modes: **CPU/OpenMP** and **CUDA**. **HIP** is
representable in the data contract, but no HIP performance claim may be made
without execution on real, suitable hardware.

Explicitly out of scope: MPI (not implemented for this scope, not a benchmark
dimension, no `mpi_ranks` field exists), mixed-precision implementation,
production kernel optimization, and hard performance gating on shared CI
machines.

## Production timing principle

Headline numbers come from the **real UppASD production executable running a
complete production workflow**. The numerator and denominator of every
reported speedup must represent equivalent physical and numerical work.

For ordinary ASD the timed unit is the complete production timestep —
`H(m_n)` → Depondt predictor → `H(m*)` → Depondt corrector — including
finite-temperature Langevin/thermal-field work whenever the selected
configuration uses it. A speedup must never be derived by comparing a
complete timestep on one backend against a subset of the timestep on another.

Three timing regimes are distinguished:

| Quantity | Meaning |
| --- | --- |
| `process_wall_seconds` | Complete process wall time: startup, input parsing, allocation, Hamiltonian and neighbour setup, GPU init, FFT planning, simulation, requested measurements, finalization. |
| `setup_seconds` | Fitted fixed-cost intercept of `T(n) = T_fixed + n·t_step`. It is an *intercept*, not a measured pure setup cost. |
| `steady_step_seconds` | Fitted steady-state production timestep — the primary compute-performance quantity. |

The fitted quantities are derived, so they are always `null` in a raw sample
record and are carried by aggregate records instead.

## Microbenchmark versus production benchmark

| | Production benchmark | Microbenchmark / profile |
| --- | --- | --- |
| What runs | The real `sd` executable, complete workflow | A kernel, phase, or instrumented region |
| What it may claim | Speedup, crossover, precision effect, throughput | Diagnostics, attribution, hypotheses |
| Where it lives | This subsystem, under the result contract | Ad-hoc tooling, `docs/` evidence notes |

Kernel benchmarks, profiler traces and phase timings are **secondary
diagnostic evidence**. They explain a production result; they never replace
one. In particular, steady-state performance is never constructed by
subtracting unrelated profiler phase timings.

## Directory responsibilities

| Directory | Responsibility | Status |
| --- | --- | --- |
| `schema/` | Versioned machine-readable result contract (JSON Schema) and its changelog. | WP-01 (this) |
| `cases/` | Immutable case manifests and input templates: physics families, variants, legal size ladders, workload metadata, input override allow-lists. | WP-02 |
| `campaigns/` | Declarative campaign manifests: which cases, sizes, builds and run configurations a LEAN or FULL campaign executes. | WP-09 |
| `harness/` | Execution and validation code: record validation (now), run isolation, provenance capture, environment-quality gating, timing (later). | WP-01 → WP-06 |
| `analysis/` | Aggregation, fitting, crossover detection, throughput metrics, plots and reports. Consumes records; never produces them. | WP-07 |
| `tests/` | Infrastructure tests: schema validation, manifest validation, harness logic, aggregation logic. No performance thresholds. | WP-01 → |

## Result-data policy

1. **Templates are immutable inputs.** Benchmark execution must never rewrite
   a tracked production template, and generated benchmark output must never be
   written into a tracked physics-template directory (`cases/`, `examples/`,
   `tests/`).
2. **Every executable run uses its own generated work directory.** Work
   directories live under gitignored paths (`benchmarks/work/`,
   `benchmarks/runs/`) or an external campaign result directory.
3. **Generated results are not tracked.** Raw records, aggregates and rendered
   analysis products are written under gitignored paths
   (`benchmarks/results/`, `benchmarks/analysis/out/`); see
   [`.gitignore`](.gitignore). What is tracked is the *contract*: schema,
   manifests, harness, analysis code and fixtures.
4. **Raw samples stay individual.** Aggregation never replaces or overwrites a
   raw record, and every aggregate names the `run_id`s it was computed from.
   Never report only the fastest sample.
5. **Provenance travels with the number.** A record either carries its full
   source/build/hardware provenance or references an immutable build record
   that does.

## Further reading

* [`TERMINOLOGY.md`](TERMINOLOGY.md) — case, variant, size, build, run
  configuration, sample, campaign.
* [`VALIDITY.md`](VALIDITY.md) — numerical versus environment validity,
  quality flags, precision semantics.
* [`schema/README.md`](schema/README.md) — the result contract and its
  versioning policy.

## Running the infrastructure tests

From the repository root:

```sh
python3 -m pytest benchmarks/tests -q
```

Requires `jsonschema` and `pytest`. A single record file can also be checked
directly:

```sh
python3 benchmarks/harness/schema_validate.py path/to/record.json
```
