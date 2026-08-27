# Benchmark result contract

The machine-readable contract every benchmark record must satisfy. It is the
interface between the harness (which writes records) and the analysis code
(which reads them), and it is what makes a number from six months ago still
interpretable.

## Files

| File | Record kind | Purpose |
| --- | --- | --- |
| `benchmark_record.v1.schema.json` | `raw_sample` | One executable measurement. Always individual. |
| `benchmark_aggregate.v1.schema.json` | `aggregate` | Statistics over the raw samples of one cell. |
| `case_manifest.v1.schema.json` | — (WP-02) | A case's variants, legal size ladder, workload metadata hook and allow-listed input overrides. Loaded by `benchmarks/harness/cases.py`, not by `schema_validate.py` — it describes a case, not a result record. |
| `omp_sweep.v1.schema.json` | — (WP-05) | A declarative OpenMP thread-sweep: thread counts, one binding policy for the whole sweep, and `physical_only`/`smt_extension` mode. Loaded by `benchmarks/harness/omp_sweep.py` — it describes intent, not a result record, and does not itself know what hardware it will run on. |
| `CHANGELOG.md` | — | Version history of the contract. |

Format is [JSON Schema draft 2020-12](https://json-schema.org/draft/2020-12/schema).
UppASD inputs are Fortran namelist-ish and its test suite uses YAML, but
neither is a validating format; JSON Schema gives real validation with a
dependency (`jsonschema`) that is already available in the project's Python
environment. Records are stored as JSON, one object per file, and remain
trivially convertible to a dataframe because the record is flat.

## Versioning policy

`schema_version` is semantic:

* **Patch** — documentation, description and `$comment` changes only.
* **Minor** — backward-compatible additions: new optional fields, widened
  enumerations. Existing records keep validating; the v1 schema file is
  updated in place.
* **Major** — any change that can invalidate an existing record: a new required
  field, a removed field, a narrowed enumeration, a tightened rule. A major
  bump creates a new schema file (`*.v2.schema.json`) and the old one is kept
  so historical records remain checkable.

The v1 schemas accept any `1.x` record. A record whose major version does not
match is rejected rather than silently coerced.

Both schemas set `additionalProperties: false`. Unknown fields are an error,
not a convenience: an unrecognized field means the writer and the contract
disagree, which is exactly the condition that silently corrupts a campaign.

## Field groups (raw sample)

| Group | Fields |
| --- | --- |
| Identity | `schema_version`, `record_kind`, `run_id`, `campaign_id`, `case_id`, `variant_id`, `size_id`, `sample_index`, `machine_id`, `timestamp_utc`, `workflow` |
| Workload | `workload_class`, `natom`, `directed_interactions`, `mean_neighbors`, `median_neighbors`, `max_neighbors`, `fft_grid`, `fft_grid_padded`, `fft_grid_points`, `dipole_backend` |
| Execution | `backend`, `gpu_backend`, `omp_threads`, `omp_places`, `omp_proc_bind`, `omp_dynamic`, `process_affinity` |
| Precision | `requested_precision`, `effective_cpu_precision`, `effective_gpu_precision`, `precision_support_state`, `comparison_precision_class` |
| Source/build | `git_commit`, `git_dirty`, `git_dirty_files`, `compiler`, `compiler_version`, `compile_flags`, `cmake_options`, `binary_checksum`, `build_id` |
| Hardware | `cpu_model`, `cpu_sockets`, `cpu_physical_cores`, `cpu_threads`, `numa_nodes`, `gpu_model`, `gpu_id`, `gpu_driver`, `gpu_runtime` |
| Simulation | `temperature`, `timestep`, `nstep`, `measurement_profile` |
| Timing | `run_status`, `timing_method`, `process_wall_seconds`, `setup_seconds`, `steady_step_seconds` |
| Quality | `numerical_valid`, `environment_valid`, `quality_flags`, `notes` |

All of these are required. Fields that do not apply to a record are present
and `null` — the contract distinguishes "not applicable" from "forgotten",
and several conditional rules turn on that distinction.

`workflow` is `ORDINARY_ASD` or the reserved `ADAPTIVE_CG`; both obey exactly
the same headline timing contract.

There is deliberately **no MPI field**. MPI is not a benchmark dimension for
this project, and `additionalProperties: false` means a stray `mpi_ranks` is
rejected rather than quietly accepted.

## Conditional rules

Encoded in the schema, so that a violating record cannot be written:

* CPU runs carry no GPU identity and no GPU precision; GPU runs must name
  `gpu_backend` (`CUDA` or `HIP`), and if they completed, the device, driver,
  runtime and effective GPU precision.
* Neighbour-list workloads must carry neighbour-workload metadata; pure
  neighbour-list workloads must not claim an FFT grid.
* FFT/dipole workloads must carry `fft_grid`, `fft_grid_points` and
  `dipole_backend`.
* `requested_precision: MIXED` ⇒ `precision_support_state: unsupported`; an
  unsupported mode is not numerically valid and has no comparison class.
* A GPU record may claim `PRECISION_MATCHED` only when the audited CPU and GPU
  precisions agree.
* `environment_valid: true` ⇒ `quality_flags` empty; `git_dirty: true` ⇒
  `dirty_tree` flag raised.
* `run_status: COMPLETED` ⇒ a real `process_wall_seconds` and a real timing
  method; any other status ⇒ neither validity flag may be true.
* `setup_seconds` and `steady_step_seconds` are reserved derived fields and are
  always `null` in a raw sample record.

## Cross-field checks

A few invariants cannot be expressed in JSON Schema and are enforced by
`benchmarks/harness/schema_validate.py` alongside the schema:

* `fft_grid_points` equals the product of the grid actually transformed;
* `mean_neighbors ≤ max_neighbors`;
* `directed_interactions == natom × mean_neighbors`;
* `sample_count` equals the number of `run_ids` in an aggregate;
* `minimum ≤ median ≤ maximum`;
* `timestamp_utc` parses as ISO 8601.

Validation therefore always goes through the harness helper, not through a
bare `jsonschema` call:

```python
from harness import schema_validate
schema_validate.validate_file("results/lean01/run-0001.json")
```

## Examples

Worked records for every supported combination live in
`benchmarks/tests/fixtures/valid/`; records that must be rejected, with the
reason each is rejected, live in `benchmarks/tests/fixtures/invalid/` and
`invalid/expected_failures.json`.
