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

Present (WP-03):

* `measurement_profiles.py` — turns `DYNAMICS_ONLY`/`PRODUCTION_LIGHT` into
  ordinary, allow-listed `inpsd.dat` overrides (section E). `DYNAMICS_ONLY`
  pushes every measurement-cadence key past the run's own `nstep` so no
  optional measurement/output work triggers within the measured run;
  `PRODUCTION_LIGHT` applies one fixed, documented cadence
  (`PRODUCTION_LIGHT_CADENCE`). Neither edits production source.
* `steady_state.py` — the `T(n) = T_fixed + n*t_step` model (section B):
  `fit_multi_nstep` (least squares over ≥3 distinct `nstep`, the
  authoritative-campaign default), `two_point_estimate` (the lean developer
  subtraction, labelled `timing_method = TWO_POINT_SUBTRACTION`, section C),
  and `calibrate_step_span` (section D: widens the sampled step-count span
  when `T(n2) - T(n1)` is too small relative to measured run-to-run jitter,
  within configured bounds).
* `runner.py` — per-run execution (section A): generates an isolated work
  directory from a case template, invokes the real executable, captures full
  wall-clock time and stdout/stderr separately, applies the section F
  validity rules (nonzero exit, fatal error, missing completion evidence,
  NaN/Inf, atom-count mismatch, missing essential output), and writes one
  schema-valid raw sample record (`run_sample`). Workload metadata is
  resolved once per `(case, variant, size)` cell via
  `resolve_workload_metadata` and reused for every sample of that cell,
  including a later sample that fails -- it describes the workload that was
  configured, not something a possibly-crashed process reported.
  Provenance/hardware/precision auditing is WP-04's job; callers supply it as
  `context` (`developer_context` is a non-authoritative local/self-test
  stand-in that always raises `metadata_incomplete`).

Planned:

| Work package | Adds |
| --- | --- |
| WP-04 | Source/build provenance, CPU/GPU metadata, contention and throttling detection, sample-quality classification |
| WP-05 | OpenMP thread sweeps, affinity and binding control |
| WP-06 | CUDA SINGLE/DOUBLE campaigns, precision audit, GPU memory-limit handling |

The harness writes records; it never edits tracked templates, and it never
writes generated data into a tracked directory.
