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

Planned:

| Work package | Adds |
| --- | --- |
| WP-02 | Case/template loading, size expansion, input override allow-lists |
| WP-03 | Isolated run directories, complete executable timing, multi-`nstep` steady-state estimation, measurement profiles |
| WP-04 | Source/build provenance, CPU/GPU metadata, contention and throttling detection, sample-quality classification |
| WP-05 | OpenMP thread sweeps, affinity and binding control |
| WP-06 | CUDA SINGLE/DOUBLE campaigns, precision audit, GPU memory-limit handling |

The harness writes records; it never edits tracked templates, and it never
writes generated data into a tracked directory.
