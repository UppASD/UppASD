"""GPU SINGLE/DOUBLE throughput ratio (WP-06 section F).

    R_GPU_32_64 = T_GPU_DOUBLE / T_GPU_SINGLE

computed over one homogeneous GPU precision family -- aggregates that share
everything except `requested_precision` (SINGLE vs DOUBLE) and `build_id`
(they are, legitimately, two different builds; see
`precision_audit.py`'s module docstring). Modelled directly on
`omp_scaling.py`'s homogeneous-family pattern: never pools aggregates that
differ in anything but the axis under study, and raises rather than silently
comparing mismatched cells.

Section F is explicit that an analogous CPU ratio may only be computed "if
both effective CPU precision modes genuinely exist". The precision audit
(`precision_audit.CPU_EFFECTIVE_PRECISION`) establishes that they do not --
`dblprec` is defined unconditionally in `source/Parameters/parameters.f90`,
so CPU is DOUBLE in every buildable configuration. `compute_r_cpu_32_64`
exists to make that refusal explicit and documented, rather than the
function simply not existing (which would look like an oversight rather
than an audited conclusion).
"""

from __future__ import annotations

# gpu_backend and omp_threads are deliberately included: a precision family
# must hold the GPU programming model and the section D host thread
# configuration fixed, varying only requested_precision. build_id is
# deliberately excluded -- SINGLE and DOUBLE are different builds by
# definition (TERMINOLOGY.md: a build is fixed by requested precision).
_PRECISION_FAMILY_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "backend", "gpu_backend", "omp_threads", "measurement_profile", "metric",
)


class PrecisionMetricsError(ValueError):
    """A precision ratio could not be computed from the given aggregates."""


def _filter_numerically_valid(aggregates, *, require_environment_valid=False):
    valid = [a for a in aggregates if a.get("all_samples_numerical_valid")]
    if require_environment_valid:
        valid = [a for a in valid if a.get("all_samples_environment_valid")]
    if not valid:
        raise PrecisionMetricsError("no numerically valid aggregates given")
    return valid


def _require_homogeneous_precision_family(aggregates):
    if not aggregates:
        raise PrecisionMetricsError("no aggregates given")
    reference = {field: aggregates[0].get(field) for field in _PRECISION_FAMILY_FIELDS}
    for aggregate in aggregates[1:]:
        mismatched = {
            field for field in _PRECISION_FAMILY_FIELDS if aggregate.get(field) != reference[field]
        }
        if mismatched:
            raise PrecisionMetricsError(
                f"aggregates do not share one GPU precision family (must differ only by "
                f"requested_precision/build_id); differing fields: {sorted(mismatched)}"
            )
    if reference["backend"] != "GPU":
        raise PrecisionMetricsError("R_GPU_32_64 applies to GPU aggregates only")
    precisions_seen = [aggregate["requested_precision"] for aggregate in aggregates]
    if len(set(precisions_seen)) != len(precisions_seen):
        raise PrecisionMetricsError(f"duplicate requested_precision among aggregates: {sorted(precisions_seen)}")
    return reference


def compute_r_gpu_32_64(aggregates, *, require_environment_valid=False):
    """`R_GPU_32_64 = T_GPU_DOUBLE / T_GPU_SINGLE` over one GPU precision family.

    Requires exactly one `SINGLE` and one `DOUBLE` aggregate among
    `aggregates` (after numeric-validity filtering); raises
    :class:`PrecisionMetricsError` otherwise, including when a `MIXED`
    aggregate is present -- `MIXED` never has real timing to compare (section
    C: it is recorded `unsupported`, never a performance result).
    """
    valid = _filter_numerically_valid(aggregates, require_environment_valid=require_environment_valid)
    _require_homogeneous_precision_family(valid)
    by_precision = {aggregate["requested_precision"]: aggregate for aggregate in valid}

    unsupported = set(by_precision) - {"SINGLE", "DOUBLE"}
    if unsupported:
        raise PrecisionMetricsError(
            f"R_GPU_32_64 requires only SINGLE/DOUBLE aggregates, got {sorted(unsupported)} as well"
        )
    missing = {"SINGLE", "DOUBLE"} - set(by_precision)
    if missing:
        raise PrecisionMetricsError(f"missing aggregate(s) for requested_precision {sorted(missing)}")

    t_single = by_precision["SINGLE"]["median"]
    t_double = by_precision["DOUBLE"]["median"]
    if not t_single > 0:
        raise PrecisionMetricsError(f"SINGLE median must be positive, got {t_single!r}")

    return {
        "r_gpu_32_64": t_double / t_single,
        "t_single": t_single,
        "t_double": t_double,
        "single_aggregate_id": by_precision["SINGLE"]["aggregate_id"],
        "double_aggregate_id": by_precision["DOUBLE"]["aggregate_id"],
    }


def compute_r_cpu_32_64(*_args, **_kwargs):
    """Deliberately unimplemented: see module docstring. There is no
    effective CPU SINGLE precision mode in this codebase to ratio against
    effective CPU DOUBLE -- `UPPASD_PRECISION` never gates `dblprec`
    (`source/Parameters/parameters.f90:8`), so a "CPU SINGLE" build is
    numerically identical to CPU DOUBLE. Computing a ratio here would imply
    a real precision effect that source inspection shows does not exist.
    """
    raise PrecisionMetricsError(
        "no analogous CPU precision ratio exists: effective_cpu_precision is DOUBLE "
        "in every buildable UPPASD_PRECISION configuration (see precision_audit.py); "
        "R_CPU_32_64 would compare a build against itself"
    )
