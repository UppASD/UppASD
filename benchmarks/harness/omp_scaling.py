"""OpenMP scaling metrics and CPU baseline selection (WP-05 sections E-F).

Operates on `aggregate` records (`harness.aggregate`), never on individual
raw samples: blueprint section 8's `S_OMP(p)`/`E_OMP(p)` and `p_best` are
defined over the per-cell statistics WP-04 already computed, not a single
noisy measurement. Every function here takes "the aggregates of one OpenMP
scaling family" -- one (campaign, case, variant, size, build, backend,
requested_precision, measurement_profile) cell varying only `omp_threads` --
and raises rather than silently pooling a heterogeneous mix.

`p_best` is never assumed to be the largest measured thread count: it is
whichever aggregate actually has the lowest median, determined empirically
from the measured aggregates given.
"""

from __future__ import annotations

import json
import pathlib

CPU_1T_THREADS = 1
CPU_1T_ROLE = "CPU-1T"
CPU_BEST_ROLE = "CPU-BEST"

# Identity fields that must agree across every aggregate entering one OpenMP
# scaling family; omp_threads is deliberately excluded -- it is the axis
# being swept -- and metric/unit are included so a caller cannot mix e.g.
# process_wall_seconds and steady_step_seconds statistics in one table.
_SCALING_CELL_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "build_id", "backend", "requested_precision", "measurement_profile", "metric",
)


class OmpScalingError(ValueError):
    """OpenMP scaling metrics could not be computed from the given aggregates."""


def _filter_numerically_valid(aggregates, *, require_environment_valid=False):
    """Section E: 'measured valid samples'. Numerically invalid aggregates
    (a claimed cell that never actually completed usable physics) are never
    considered; environment validity is stricter and only required when the
    caller asks for authoritative-only evidence.
    """
    valid = [a for a in aggregates if a.get("all_samples_numerical_valid")]
    if require_environment_valid:
        valid = [a for a in valid if a.get("all_samples_environment_valid")]
    if not valid:
        raise OmpScalingError("no numerically valid aggregates given")
    return valid


def _require_homogeneous_scaling_family(aggregates):
    if not aggregates:
        raise OmpScalingError("no aggregates given")
    reference = {field: aggregates[0].get(field) for field in _SCALING_CELL_FIELDS}
    for aggregate in aggregates[1:]:
        mismatched = {
            field for field in _SCALING_CELL_FIELDS if aggregate.get(field) != reference[field]
        }
        if mismatched:
            raise OmpScalingError(
                f"aggregates do not share one OpenMP-scaling family (must differ only by "
                f"omp_threads); differing fields: {sorted(mismatched)}"
            )
    if reference["backend"] != "CPU":
        raise OmpScalingError("OpenMP scaling metrics apply to CPU aggregates only")
    threads_seen = [aggregate["omp_threads"] for aggregate in aggregates]
    if len(set(threads_seen)) != len(threads_seen):
        raise OmpScalingError(f"duplicate omp_threads among aggregates: {sorted(threads_seen)}")
    return reference


def compute_omp_speedup_table(aggregates, *, require_environment_valid=False):
    """Return one row per thread count: `omp_threads`, `t_p` (the aggregate's
    median), `S_OMP(p) = T(N,1)/T(N,p)`, `E_OMP(p) = S_OMP(p)/p`.

    Requires an `omp_threads == 1` aggregate (CPU-1T) among `aggregates` as
    the reference; raises :class:`OmpScalingError` otherwise -- there is no
    meaningful speedup without one.
    """
    valid = _filter_numerically_valid(aggregates, require_environment_valid=require_environment_valid)
    _require_homogeneous_scaling_family(valid)
    by_threads = {aggregate["omp_threads"]: aggregate for aggregate in valid}
    if CPU_1T_THREADS not in by_threads:
        raise OmpScalingError("no omp_threads=1 (CPU-1T) aggregate present; cannot compute S_OMP(p)")
    t1 = by_threads[CPU_1T_THREADS]["median"]
    if not t1 > 0:
        raise OmpScalingError(f"CPU-1T median must be positive, got {t1!r}")

    rows = []
    for threads in sorted(by_threads):
        aggregate = by_threads[threads]
        t_p = aggregate["median"]
        speedup = t1 / t_p
        efficiency = speedup / threads
        rows.append({
            "omp_threads": threads,
            "aggregate_id": aggregate["aggregate_id"],
            "t_p": t_p,
            "speedup": speedup,
            "efficiency": efficiency,
        })
    return rows


def select_cpu_1t(aggregates, *, require_environment_valid=False):
    """Return the `omp_threads == 1` aggregate of this scaling family."""
    valid = _filter_numerically_valid(aggregates, require_environment_valid=require_environment_valid)
    _require_homogeneous_scaling_family(valid)
    by_threads = {aggregate["omp_threads"]: aggregate for aggregate in valid}
    if CPU_1T_THREADS not in by_threads:
        raise OmpScalingError("no omp_threads=1 (CPU-1T) aggregate present")
    return by_threads[CPU_1T_THREADS]


def determine_p_best(aggregates, *, require_environment_valid=False):
    """Return the fastest *measured* aggregate of this scaling family --
    i.e. CPU-BEST. Never assumes the largest thread count wins: this is a
    plain empirical minimum over the medians actually measured.
    """
    valid = _filter_numerically_valid(aggregates, require_environment_valid=require_environment_valid)
    _require_homogeneous_scaling_family(valid)
    return min(valid, key=lambda aggregate: aggregate["median"])


def select_cpu_best(aggregates, *, require_environment_valid=False):
    return determine_p_best(aggregates, require_environment_valid=require_environment_valid)


def build_baseline_pointer(aggregate, role):
    """A small, stable-shaped pointer at one canonical CPU baseline aggregate.

    Deliberately not itself an `aggregate` record (adding `baseline_role`
    would violate `benchmark_aggregate.v1.schema.json`'s
    `additionalProperties: false`): this is a convenience lookup artifact for
    later work (WP-06's GPU crossover needs CPU-BEST by cell, not by
    searching every aggregate), not a new result-contract record kind.
    """
    if role not in (CPU_1T_ROLE, CPU_BEST_ROLE):
        raise OmpScalingError(f"unknown baseline role {role!r}; expected {CPU_1T_ROLE!r} or {CPU_BEST_ROLE!r}")
    return {
        "baseline_role": role,
        "aggregate_id": aggregate["aggregate_id"],
        "campaign_id": aggregate["campaign_id"],
        "case_id": aggregate["case_id"],
        "variant_id": aggregate["variant_id"],
        "size_id": aggregate["size_id"],
        "machine_id": aggregate["machine_id"],
        "build_id": aggregate["build_id"],
        "metric": aggregate["metric"],
        "omp_threads": aggregate["omp_threads"],
        "median_seconds": aggregate["median"],
        "authoritative": aggregate.get("authoritative"),
    }


def _baseline_filename(aggregate, role):
    return (
        f"{aggregate['campaign_id']}__{aggregate['case_id']}__{aggregate['variant_id']}__"
        f"{aggregate['size_id']}__{role}.json"
    )


def persist_cpu_baseline(aggregate, role, output_dir):
    """Write `build_baseline_pointer(aggregate, role)` under a stable,
    predictable filename in `output_dir`. Returns the written path.
    """
    pointer = build_baseline_pointer(aggregate, role)
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / _baseline_filename(aggregate, role)
    path.write_text(json.dumps(pointer, indent=2, sort_keys=True) + "\n")
    return path
