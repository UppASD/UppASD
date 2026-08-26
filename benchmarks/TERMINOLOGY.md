# Benchmark terminology

These seven terms are the axes of the framework. Every one of them maps onto
fields of the result contract, so that any record can be traced back to
exactly what was run. Two records may only be compared directly when they
differ in the dimensions under study and agree in the rest.

## CASE

A **physical input/model family** — a Hamiltonian class and lattice, not a
particular input file.

Examples: `B01_bccFe` (medium-range exchange), `B02_skyrmion2D` (short-range
J+D), `B04_dhcpNd` (very long-range), `B05_dipoleFFT`.

A case owns immutable input templates. Its physics is not simplified to make
benchmarking convenient: interaction shells are not truncated, cutoffs are not
reduced and weak interactions are not discarded.

Record field: `case_id`.

## VARIANT

An **alternative physics or runtime configuration of the same case**: the same
model family exercised in a materially different way.

Examples: `T0` versus `finiteT` for bcc Fe; `dipole_only` versus
`exchange_plus_dipole` for the FFT case.

A variant may change temperature, which Hamiltonian terms are active, or the
integrator regime. It does not change the case's identity, and it is not a
mechanism for silently reducing work.

Record field: `variant_id`.

## SIZE

A **legal geometric replication of the case** — a point on the case's declared
size ladder.

Examples: `16x16x16`, `32x32x32`, `64x64x16`.

Legality is a property of the case: only replications that preserve the
intended physics and periodicity are admitted. Size is not a scalar. Every
neighbour-list record therefore also carries `natom`,
`directed_interactions`, `mean_neighbors` and `max_neighbors`; every
FFT/dipole record also carries `fft_grid` and `fft_grid_points`.

Record fields: `size_id`, plus the workload-metadata block.

## BUILD

A **compiled binary configuration**: one executable, identified immutably.

A build is fixed by source state (`git_commit`, `git_dirty`,
`git_dirty_files`), toolchain (`compiler`, `compiler_version`,
`compile_flags`, `cmake_options`), backend and requested precision, and is
identified by `build_id` and `binary_checksum`.

Two runs from different builds are never pooled. Changing an optimization flag
produces a new build, hence a new `build_id`.

Record fields: `build_id`, `binary_checksum`, and the source/build block.

## RUN CONFIGURATION

The **runtime choices** applied to a given build on a given machine: which
backend is used, the OpenMP environment (`omp_threads`, `omp_places`,
`omp_proc_bind`, `omp_dynamic`, `process_affinity`), the measurement profile
(`DYNAMICS_ONLY` or `PRODUCTION_LIGHT`) and the simulation length (`nstep`).

Two canonical CPU run configurations are maintained per case and size:

* **CPU-1T** — one OpenMP thread; used for parallel efficiency and historical
  comparison.
* **CPU-BEST** — the fastest *valid measured* supported OpenMP configuration.
  This is the principal CPU baseline for GPU crossover. A GPU speedup against
  CPU-1T may be reported, but never in place of GPU versus CPU-BEST.

Record fields: `backend`, `gpu_backend`, the OpenMP block,
`measurement_profile`, `nstep`, `machine_id`.

## SAMPLE

**One executable measurement** — one complete process run, producing exactly
one raw record.

Samples are individual and stay individual. Repeated samples of the same cell
differ only in `sample_index` and `run_id`. Statistics over them
(`sample_count`, `median`, `mad`, `minimum`, `maximum`) live in a separate
aggregate record that names its contributing `run_id`s.

Record fields: `run_id`, `sample_index`, `timestamp_utc`, the timing block.

## CAMPAIGN

A **declarative set of run configurations**: which cases, variants, sizes,
builds and run configurations to execute, at how many samples, under which
environment-quality policy.

Two tiers exist:

* **LEAN** — harness development, developer sanity, gross regression
  detection. Limited cases, sizes, thread counts and sample counts. Every lean
  report must state prominently:
  **LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS.**
* **FULL / AUTHORITATIVE** — all core workload families, practical size
  ladders, broad OpenMP scaling, CUDA SINGLE and DOUBLE where supported,
  sufficient independent samples, strict provenance and strict
  environment-quality gating.

Record field: `campaign_id`; tier is a property of the campaign manifest
(WP-09).

## Relationship

```text
campaign
  └── case
        └── variant
              └── size
                    └── build  ×  run configuration
                          └── sample  ×  N   →  raw records
                                                  └── aggregate record
```

A **cell** is one (campaign, case, variant, size, build, run configuration)
combination. Aggregation happens strictly within a cell.
