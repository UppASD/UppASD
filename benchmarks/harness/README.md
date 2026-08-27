# Benchmark harness

Code that validates, and later executes and instruments, benchmark runs.

Present (WP-01):

* `schema_validate.py` — validates raw sample and aggregate records against
  the versioned contract in [`../schema/`](../schema/), applying both JSON
  Schema and the cross-field invariants that JSON Schema cannot express. Usable
  as a library (`validate_record`, `validate_file`) or as a CLI:

  ```sh
  python3 benchmarks/harness/schema_validate.py results/lean01/*.json
  ```

Present (WP-02):

* `cases.py` — loads and validates case manifests against
  [`../schema/case_manifest.v1.schema.json`](../schema/case_manifest.v1.schema.json),
  resolves variants/sizes, enforces the allow-listed input-override vocabulary
  (`GLOBALLY_SAFE_OVERRIDE_KEYS`), copies an immutable template into a
  generated work directory with overrides applied
  (`generate_run_directory`), and computes a deterministic SHA-256 template
  fingerprint (`compute_template_fingerprint`). Also provides the campaign-tier
  guard (`filter_cases_for_campaign`, `require_not_infrastructure_only`) that
  keeps `infrastructure_test_only` cases out of authoritative reports.
* `workload_metadata.py` — the case-specific hook mechanism of section E:
  `neighbor_list_from_struct_output` parses UppASD's own
  `struct.<simid>.out` diagnostic dump for real neighbour-workload metadata;
  `fft_grid_from_replication` derives the dipole-FFT grid from the case's own
  supercell replication using UppASD's own padding formula. Neither estimates
  a quantity real production input/output already supplies.

Planned:

| Work package | Adds |
| --- | --- |
| WP-03 | Isolated run directories, complete executable timing, multi-`nstep` steady-state estimation, measurement profiles |
| WP-04 | Source/build provenance, CPU/GPU metadata, contention and throttling detection, sample-quality classification |
| WP-05 | OpenMP thread sweeps, affinity and binding control |
| WP-06 | CUDA SINGLE/DOUBLE campaigns, precision audit, GPU memory-limit handling |

The harness writes records; it never edits tracked templates, and it never
writes generated data into a tracked directory.
