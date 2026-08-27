# Benchmark result contract changelog

Versioning policy: see [README.md](README.md).

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
