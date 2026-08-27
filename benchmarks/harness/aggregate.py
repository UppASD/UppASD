"""Sample-statistics aggregation (WP-04 section F).

Raw sample records are never overwritten or replaced -- every one of them is
already an individual file on disk by the time this module runs. What this
module adds is the *aggregate* record kind: statistics computed strictly
within one cell (blueprint TERMINOLOGY.md: one campaign/case/variant/size/
build/run-configuration combination), naming every contributing `run_id`.

Two entry points, matching the two ways a `steady_step_seconds`/
`setup_seconds` value differs from `process_wall_seconds`:

* `build_process_wall_aggregate` aggregates the `process_wall_seconds` field
  directly present on completed raw sample records.
* `build_fit_aggregate` aggregates across several independent
  `steady_state.FitResult`s (each itself already derived from several
  `process_wall_seconds` samples at different `nstep`) -- a raw sample record
  never carries a non-null `steady_step_seconds`/`setup_seconds` itself (the
  schema reserves those fields), so there is nothing to pull directly off one
  record for those metrics.

Never reports only the fastest sample: both entry points always return the
full `count`/`median`/`mad`/`minimum`/`maximum` block, and callers are not
offered a "just give me the min" shortcut.
"""

from __future__ import annotations

import statistics as _statistics
from collections import namedtuple

from harness import schema_validate

SCHEMA_VERSION = "1.1.0"

_CELL_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "build_id", "backend", "gpu_backend", "omp_threads", "requested_precision",
    "measurement_profile",
)

SampleStatistics = namedtuple("SampleStatistics", "count median mad minimum maximum")


class AggregateError(ValueError):
    """Samples/fits could not be aggregated into one schema-valid cell."""


def compute_sample_statistics(values):
    """`count`/`median`/`mad`/`minimum`/`maximum` over `values`.

    Never collapses to "the minimum" alone (blueprint section 13: "Do not
    report only the fastest sample") -- every field is always returned
    together. Raises :class:`AggregateError` on an empty input rather than
    fabricating statistics for zero samples.
    """
    values = [float(v) for v in values]
    if not values:
        raise AggregateError("cannot compute sample statistics over zero values")
    med = _statistics.median(values)
    mad = _statistics.median(abs(v - med) for v in values)
    return SampleStatistics(
        count=len(values), median=med, mad=mad, minimum=min(values), maximum=max(values)
    )


def _require_homogeneous_cell(records):
    if not records:
        raise AggregateError("no records to aggregate")
    reference = {field: records[0].get(field) for field in _CELL_FIELDS}
    for record in records[1:]:
        mismatched = {
            field for field in _CELL_FIELDS if record.get(field) != reference[field]
        }
        if mismatched:
            raise AggregateError(
                f"records do not share one cell; differing fields: {sorted(mismatched)} "
                f"(run_ids {records[0].get('run_id')!r} vs {record.get('run_id')!r})"
            )
    return reference


def _base_aggregate_fields(
    cell, *, aggregate_id, metric, run_ids,
    all_numerical_valid, all_environment_valid, quality_flags_union,
):
    authoritative = bool(all_numerical_valid and all_environment_valid and not quality_flags_union)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_kind": "aggregate",
        "aggregate_id": aggregate_id,
        "campaign_id": cell["campaign_id"],
        "case_id": cell["case_id"],
        "variant_id": cell["variant_id"],
        "size_id": cell["size_id"],
        "machine_id": cell["machine_id"],
        "build_id": cell["build_id"],
        "backend": cell["backend"],
        "gpu_backend": cell["gpu_backend"],
        "omp_threads": cell["omp_threads"],
        "requested_precision": cell["requested_precision"],
        "measurement_profile": cell["measurement_profile"],
        "metric": metric,
        "unit": "s",
        "run_ids": sorted(run_ids),
        "all_samples_numerical_valid": all_numerical_valid,
        "all_samples_environment_valid": all_environment_valid,
        "quality_flags_union": sorted(quality_flags_union),
        "authoritative": authoritative,
        "notes": None,
    }


def build_process_wall_aggregate(records, *, aggregate_id):
    """Aggregate `process_wall_seconds` over the `COMPLETED` records of one cell.

    `records` must be the complete set of raw sample records for one
    (campaign, case, variant, size, build, run configuration, nstep) cell --
    every one of them, not a pre-filtered "best of" subset. Records whose
    `run_status` is not `COMPLETED` carry no timing (schema rule) and are
    excluded from the statistics themselves, but every completed one
    contributes its `run_id`; validity/quality-flag unions are likewise taken
    over the completed contributors, since an excluded, non-completed sample
    contributes no timing evidence to summarize.
    """
    records = list(records)
    completed = [r for r in records if r.get("run_status") == "COMPLETED"]
    if not completed:
        raise AggregateError("no COMPLETED records to aggregate process_wall_seconds from")
    cell = _require_homogeneous_cell(completed)

    nsteps = {r.get("nstep") for r in completed}
    if len(nsteps) != 1:
        raise AggregateError(
            f"process_wall_seconds aggregation needs one shared nstep, got {sorted(nsteps)}"
        )

    values = [r["process_wall_seconds"] for r in completed]
    stats = compute_sample_statistics(values)

    quality_flags_union = set()
    for record in completed:
        quality_flags_union.update(record.get("quality_flags", []))

    fields = _base_aggregate_fields(
        cell, aggregate_id=aggregate_id, metric="process_wall_seconds",
        run_ids=[r["run_id"] for r in completed],
        all_numerical_valid=all(r.get("numerical_valid") for r in completed),
        all_environment_valid=all(r.get("environment_valid") for r in completed),
        quality_flags_union=quality_flags_union,
    )
    fields.update({
        "nstep": next(iter(nsteps)),
        "nstep_fit_points": None,
        "fit_method": None,
        "sample_count": stats.count,
        "median": stats.median,
        "mad": stats.mad,
        "minimum": stats.minimum,
        "maximum": stats.maximum,
    })

    schema_validate.validate_record(fields)
    return fields


def build_fit_aggregate(
    fits, *, metric, run_ids, cell, aggregate_id,
    all_numerical_valid, all_environment_valid, quality_flags_union=(),
):
    """Aggregate several independent `steady_state.FitResult`s of the same cell.

    `metric` is `"steady_step_seconds"` or `"setup_seconds"`; the
    corresponding attribute (`steady_step_seconds` or
    `setup_or_fixed_intercept`) is pulled from every fit. All fits must agree
    on `fit_method` and `nstep_fit_points` -- they are repeats of the same
    measurement protocol, not different protocols pooled together.

    `run_ids` must give exactly one identifier per fit, in the same order --
    the schema's `sample_count == len(run_ids)` invariant (enforced by
    `schema_validate`) ties `run_ids` to whatever the statistics were
    actually computed over, and here that is one fit per repetition of the
    multi-nstep protocol, not the flattened set of every raw
    `process_wall_seconds` sample that fed into all of them. A stable id for
    each repetition (e.g. its lowest-`nstep` raw `run_id`) is enough to trace
    a fit back to its raw samples via that raw record's own fields; there is
    no way to recover the full raw-sample set from a `FitResult` alone, so
    it is not attempted here. `cell` supplies the shared
    campaign/case/.../measurement_profile identity fields, since fits carry
    none of that either.
    """
    fits = list(fits)
    if not fits:
        raise AggregateError("cannot build a fit aggregate from zero fits")
    run_ids = list(run_ids)
    if len(run_ids) != len(fits):
        raise AggregateError(
            f"run_ids must give exactly one id per fit ({len(fits)} fits, {len(run_ids)} run_ids)"
        )
    if metric == "steady_step_seconds":
        values = [f.steady_step_seconds for f in fits]
    elif metric == "setup_seconds":
        values = [f.setup_or_fixed_intercept for f in fits]
    else:
        raise AggregateError(f"build_fit_aggregate does not support metric {metric!r}")

    fit_methods = {f.fit_method for f in fits}
    if len(fit_methods) != 1:
        raise AggregateError(f"fits do not share one fit_method: {sorted(fit_methods)}")
    nstep_points = {tuple(f.nstep_fit_points) for f in fits}
    if len(nstep_points) != 1:
        raise AggregateError(f"fits do not share one nstep_fit_points set: {sorted(nstep_points)}")

    stats = compute_sample_statistics(values)
    if len(set(run_ids)) != len(run_ids):
        raise AggregateError(f"run_ids must be unique, got {run_ids}")

    fields = _base_aggregate_fields(
        cell, aggregate_id=aggregate_id, metric=metric, run_ids=sorted(run_ids),
        all_numerical_valid=all_numerical_valid, all_environment_valid=all_environment_valid,
        quality_flags_union=set(quality_flags_union),
    )
    fields.update({
        "nstep": None,
        "nstep_fit_points": sorted(next(iter(nstep_points))),
        "fit_method": next(iter(fit_methods)),
        "sample_count": stats.count,
        "median": stats.median,
        "mad": stats.mad,
        "minimum": stats.minimum,
        "maximum": stats.maximum,
    })

    schema_validate.validate_record(fields)
    return fields
