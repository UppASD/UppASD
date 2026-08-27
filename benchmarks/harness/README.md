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
  Provenance/hardware/precision auditing is supplied as `context`:
  `provenance.build_static_context` (WP-04) is the authoritative source;
  `developer_context` remains a non-authoritative local/self-test stand-in
  that always raises `metadata_incomplete`. `run_sample`'s
  `require_clean_environment` argument is WP-04's strict mode.

Present (WP-04):

* `source_provenance.py` — section A: git commit/dirty state/dirty files
  (`capture_git_provenance`), the executable's SHA-256
  (`compute_binary_checksum`), and everything read from the build's own
  `CMakeCache.txt` (`parse_cmake_cache`/`build_source_context`) --
  `requested_precision` and `gpu_backend` come from `UPPASD_PRECISION`/
  `UPPASD_GPU_BACKEND` there, never guessed from the binary's filename.
* `cpu_provenance.py` — section B: CPU topology from `/proc/cpuinfo` and
  `/sys/devices/system/node` (`capture_cpu_topology`), the OpenMP
  environment (`capture_omp_environment`), process affinity
  (`capture_process_affinity`), per-core governor/frequency
  (`capture_cpu_frequency_state`) and load average
  (`capture_background_load`). Every filesystem root is a parameter, so
  tests run against a synthetic tree rather than the real machine.
* `gpu_provenance.py` — section C: CUDA device identity/state via
  `nvidia-smi` and toolkit version via `nvcc --version`
  (`query_cuda_devices`, `query_cuda_runtime_version`,
  `query_cuda_compute_processes`), plus contamination/throttle/clock-drift
  classification (`classify_cuda_contamination`,
  `detect_cuda_clock_instability`). HIP exposes the structurally analogous
  functions against `rocm-smi --showallinfo --json`
  (`query_hip_devices`, `classify_hip_contamination`) -- representable and
  tested via mock parsers without requiring AMD hardware, per the
  blueprint's "no HIP performance claim without real hardware" rule. Every
  external call goes through an injectable `run` callable.
* `provenance.py` — orchestrates the above into one `context` dict shaped
  exactly like `runner.run_sample`'s (`build_static_context`, resolved once
  per build+machine, the authoritative counterpart to
  `runner.developer_context`); section D's environment-quality flags
  (`classify_background_load`, `classify_cpu_frequency_stability`,
  `classify_gpu_environment`, bracketing a sample with
  `capture_gpu_snapshot` before and after); and section E's strict mode
  (`gate_sample`, raising `EnvironmentGateError` -- surfaced through
  `runner.run_sample`'s `require_clean_environment` argument, `False` by
  default so ordinary developer runs, including against a dirty tree, are
  never blocked).
* `aggregate.py` — section F: `compute_sample_statistics` (count/median/
  MAD/min/max, never just the fastest sample) and the two aggregate-record
  builders, `build_process_wall_aggregate` (directly over completed raw
  samples) and `build_fit_aggregate` (over several independent
  `steady_state.FitResult`s of the same cell, for `steady_step_seconds`/
  `setup_seconds`). Aggregation is strictly within one cell -- a
  heterogeneous mix of `run_id`s raises rather than silently pooling.

Planned:

| Work package | Adds |
| --- | --- |
| WP-05 | OpenMP thread sweeps, affinity and binding control |
| WP-06 | CUDA SINGLE/DOUBLE campaigns, precision audit, GPU memory-limit handling |

The harness writes records; it never edits tracked templates, and it never
writes generated data into a tracked directory.
