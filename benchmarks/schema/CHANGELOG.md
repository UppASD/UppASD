# Benchmark result contract changelog

Versioning policy: see [README.md](README.md).

## campaign_manifest 1.0.0 — WP-09

Independent contract with its own `schema_version`; not a result record kind.

* `campaign_manifest.v1.schema.json` — declarative campaign manifests:
  `tier` (`LEAN`/`FULL`), `backend_scope` (`CPU_ONLY`/`CUDA_ONLY`/
  `CPU_AND_CUDA`), case/variant/size selection (with an `"ALL"` wildcard),
  a symbolic CPU `thread_policy` (never a concrete thread-count list —
  resolved against real host topology only once bound to a machine, via
  `benchmarks/harness/campaigns.py::resolve_cpu_sweep`, the same discipline
  `omp_sweep.py` established), GPU precisions, `sample_count`, `measurement_profile` and
  `environment_quality_mode`. Cross-field rules tie the tiers to the
  blueprint: `LEAN` requires the exact non-authoritative `report_banner`
  text and `DEVELOPER` environment mode; `FULL` forbids the banner, requires
  `sample_count >= 5` and `STRICT` environment mode; `CPU_AND_CUDA` requires
  an explicit `single_campaign_execution: true` acknowledgement that its CPU
  and CUDA cells run together under one `campaign_id`/build. Loaded by
  `benchmarks/harness/campaigns.py`, not by `schema_validate.py`.

## omp_sweep 1.0.0 — WP-05

Independent contract with its own `schema_version`; not a result record kind.

* `omp_sweep.v1.schema.json` — declarative OpenMP thread-sweep manifests:
  `thread_counts` (must include 1, so CPU-1T can be established), one
  `proc_bind` (`close` or `spread`) for the whole sweep rather than a blind
  sweep of both, `places` fixed to `cores` and `dynamic` fixed to `false`
  (blueprint section 9/C), and `smt_mode` (`physical_only` versus the
  separately-enabled `smt_extension`). Loaded by
  `benchmarks/harness/omp_sweep.py`, not by `schema_validate.py`.

## 1.1.0 — WP-04

* `benchmark_aggregate.v1.schema.json` — the `median`/`minimum`/`maximum`
  positivity constraint (`exclusiveMinimum: 0`) now applies only to metrics
  other than `setup_seconds`. The schema's own `setup_seconds` conditional
  already documented that a fitted fixed-cost intercept "may legitimately be
  negative or zero" (blueprint section 6.2), but the constraint had
  incorrectly stayed unconditional, silently rejecting exactly the records
  that comment described. Backward compatible: every previously valid record
  necessarily had positive statistics already and keeps validating; only
  newly permits a `setup_seconds` aggregate with a non-positive statistic.

## case_manifest 1.0.0 — WP-02

Independent contract with its own `schema_version`; not a result record kind.

* `case_manifest.v1.schema.json` — declarative case manifests: variants,
  legal size ladder (2D/3D/case-specific), `expected_physics` documentation,
  the allow-listed input-override vocabulary, and the workload-metadata hook
  name. Loaded by `benchmarks/harness/cases.py`, not by `schema_validate.py`.

## 1.0.0 — WP-01

Initial result contract.

* `benchmark_record.v1.schema.json` — raw sample records, covering identity,
  workload (neighbour-list and FFT/dipole), execution, precision,
  source/build provenance, hardware, simulation, timing and quality.
* `benchmark_aggregate.v1.schema.json` — aggregate records with
  `sample_count`, `median`, `mad`, `minimum`, `maximum` and the `run_ids` they
  were computed from.
* Requested and effective precision separated, with `MIXED` representable as an
  unsupported production mode.
* `numerical_valid` separated from `environment_valid`, with a closed set of
  environment quality flags.
* CUDA and HIP both representable; HIP without any measurement is representable
  as a skipped run carrying no timing.
* No MPI fields: MPI is not a benchmark dimension for this project, and unknown
  fields are rejected.
