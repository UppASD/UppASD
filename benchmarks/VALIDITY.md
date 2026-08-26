# Validity, environment quality and precision semantics

## Two independent validity questions

The contract keeps two questions apart, because they have different answers
and different consequences.

| Field | Question | Consequence if false |
| --- | --- | --- |
| `numerical_valid` | Did the run complete and produce physically/numerically acceptable output? | The record is not evidence of anything. |
| `environment_valid` | Was the measurement environment suitable for authoritative timing? | The record is valid physics, but not authoritative performance evidence. |

The canonical example: a simulation that completes correctly while another
process occupies the GPU is **numerically valid** and **not environment
valid**. Its result is a real result; it is simply not admissible as
performance evidence.

Both are recorded per sample. Neither is inferred from the other.

## Quality flags

`quality_flags` is the closed set of reasons a sample's environment is
unsuitable:

| Flag | Raised when |
| --- | --- |
| `dirty_tree` | The source tree had uncommitted tracked modifications at build time. |
| `gpu_busy` | Another process was using the target GPU during the run. |
| `gpu_thermal_throttle` | The GPU was thermally throttled. |
| `gpu_clock_unstable` | GPU clocks varied beyond the accepted band. |
| `cpu_affinity_unknown` | The process affinity could not be determined. |
| `cpu_frequency_unstable` | CPU frequency varied beyond the accepted band. |
| `background_load_high` | Competing load was present on the machine. |
| `metadata_incomplete` | Required provenance could not be captured. |

Rules enforced by the schema:

* `environment_valid: true` requires an **empty** `quality_flags` array —
  environment validity is the conjunction, so any raised flag clears it.
* `git_dirty: true` **must** raise `dirty_tree`.
* A run whose `run_status` is not `COMPLETED` is neither numerically nor
  environmentally valid.

## Strict and developer modes

* **Authoritative campaigns** run in strict mode: samples that are not
  `environment_valid` are refused, not merely annotated.
* **Developer runs** remain usable with flags present. A flagged sample is
  still a record; it simply may not back a headline claim.

Detection of these conditions is WP-04. WP-01 fixes only the vocabulary and
the consistency rules.

## Precision semantics

Precision is **not** inferred from CMake option names. The contract therefore
separates four things:

| Field | Meaning |
| --- | --- |
| `requested_precision` | What the build configuration asked for: `SINGLE`, `DOUBLE` or `MIXED`. |
| `effective_cpu_precision` | What the CPU numerical path actually uses, as established by audit. |
| `effective_gpu_precision` | What the GPU numerical path actually uses, as established by audit. |
| `comparison_precision_class` | How the record may be used in a cross-backend comparison. |

The effective fields may legitimately disagree for a single requested value.
The contract explicitly supports:

```json
{
  "requested_precision": "SINGLE",
  "effective_cpu_precision": "DOUBLE",
  "effective_gpu_precision": "SINGLE"
}
```

if that is what auditing establishes. `UNKNOWN` is available before an audit
has been performed; `null` means not applicable (for example the GPU precision
of a CPU-only run).

### Comparison classes

* **`PRECISION_MATCHED`** — the relevant CPU and GPU numerical paths use
  genuinely corresponding, audited precision. A GPU record may only claim this
  when `effective_cpu_precision` and `effective_gpu_precision` agree.
* **`PRODUCTION_CONFIGURATION`** — a comparison between real supported
  production modes that are *not* precision matched, e.g. CPU double against
  CUDA single. Operationally useful, and must never be labelled
  precision-matched.
* **`UNAUDITED`** — the precision audit has not been performed; the record may
  not back a precision claim.
* **`null`** — the record cannot enter a precision comparison at all.

### MIXED is unsupported, not broken

Mixed precision is deliberately not implemented in production. It is therefore
representable as a **first-class unsupported state** rather than as a generic
failed build:

```json
{
  "requested_precision": "MIXED",
  "precision_support_state": "unsupported",
  "run_status": "SKIPPED",
  "numerical_valid": false,
  "comparison_precision_class": null
}
```

`precision_support_state` takes `supported`, `unsupported` or `unaudited`. The
schema requires `MIXED` to be recorded as `unsupported`, and an unsupported
mode can neither be numerically valid nor enter a precision comparison.
