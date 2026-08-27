"""GPU/CPU crossover detection (WP-07 sections A-B, blueprint section 10).

Two things live here:

1. `compute_gpu_speedup` -- the ratio itself, `S = T_cpu.median / T_gpu.median`,
   for one matched (campaign, case, variant, size, machine, measurement_profile,
   metric) cell. Used for both the headline `S_GPU_BESTCPU` (`cpu_aggregate` is
   the case/size's CPU-BEST aggregate) and the secondary `S_GPU_1T`
   (`cpu_aggregate` is CPU-1T) -- the caller picks which CPU aggregate to pass;
   this module does not itself select CPU-BEST/CPU-1T (that is
   `harness.omp_scaling.determine_p_best`/`select_cpu_1t`).

2. `find_crossover` -- given a speedup curve (one point per tested size, sorted
   by `natom`), determines the smallest `natom` at which the curve reaches a
   threshold such as 1.0x/1.25x/2.0x/5.0x. Interpolation is allowed only
   between two measured neighbouring points that bracket the threshold
   (blueprint section 10: "Interpolation is allowed only between valid
   measured neighbouring points. Never extrapolate outside the measured range
   and present the result as a measured crossover"). When the threshold is
   already exceeded at the smallest tested size, or never reached at the
   largest, the true crossover lies outside the tested range and no
   interpolated `natom` is reported -- only a `below_tested_range` /
   `above_tested_range` classification.

Interpolation is done in log(natom) vs log(speedup) space: both throughput and
GPU/CPU speedup are typically close to power-law in problem size over the
scale a size ladder spans, so a straight line in log-log space is the more
defensible model than linear interpolation in the raw values. This is a
documented modelling choice, not a measured fact (compare
`omp_topology.classify_numa_placement`'s documented-model precedent).
"""

from __future__ import annotations

import math

BELOW_TESTED_RANGE = "below_tested_range"
WITHIN_TESTED_RANGE = "within_tested_range"
ABOVE_TESTED_RANGE = "above_tested_range"

DEFAULT_CROSSOVER_THRESHOLDS = (1.0, 1.25, 2.0, 5.0)

# Interpolated crossover natoms are rounded to this many significant figures
# before being reported, so an interpolated estimate is never displayed with
# more precision than the two measured points bracketing it actually support
# (blueprint section F: "Do not display interpolated crossover with spurious
# numerical precision").
_CROSSOVER_DISPLAY_SIG_FIGS = 2

_SHARED_CELL_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "measurement_profile", "metric",
)


class CrossoverError(ValueError):
    """A speedup or crossover could not be computed from the given inputs."""


def _require_positive(value, name):
    if value is None or not value > 0:
        raise CrossoverError(f"{name} must be a positive number, got {value!r}")


def round_significant(value, digits=_CROSSOVER_DISPLAY_SIG_FIGS):
    """Round `value` to `digits` significant figures. `value == 0` passes through."""
    if value == 0:
        return 0.0
    magnitude = math.floor(math.log10(abs(value)))
    factor = 10 ** (digits - 1 - magnitude)
    return round(value * factor) / factor


def compute_gpu_speedup(cpu_aggregate, gpu_aggregate, *, label):
    """`S = T_cpu.median / T_gpu.median` for one matched cell.

    `label` (e.g. `"S_GPU_BESTCPU"` or `"S_GPU_1T"`) is carried through into
    the result and into error messages only -- it does not change behaviour.
    Both aggregates must share campaign/case/variant/size/machine/
    measurement_profile/metric (`metric` is normally `"steady_step_seconds"`,
    the blueprint's primary compute-performance quantity); `cpu_aggregate`
    must be a `CPU` aggregate and `gpu_aggregate` a `GPU` aggregate.
    """
    if cpu_aggregate.get("backend") != "CPU":
        raise CrossoverError(f"{label}: cpu_aggregate backend must be CPU, got {cpu_aggregate.get('backend')!r}")
    if gpu_aggregate.get("backend") != "GPU":
        raise CrossoverError(f"{label}: gpu_aggregate backend must be GPU, got {gpu_aggregate.get('backend')!r}")
    mismatched = {
        field for field in _SHARED_CELL_FIELDS
        if cpu_aggregate.get(field) != gpu_aggregate.get(field)
    }
    if mismatched:
        raise CrossoverError(
            f"{label}: cpu/gpu aggregates do not share one cell; differing fields: {sorted(mismatched)}"
        )
    t_cpu = cpu_aggregate["median"]
    t_gpu = gpu_aggregate["median"]
    _require_positive(t_cpu, f"{label}: cpu median")
    _require_positive(t_gpu, f"{label}: gpu median")

    return {
        "label": label,
        "size_id": gpu_aggregate["size_id"],
        "speedup": t_cpu / t_gpu,
        "t_cpu": t_cpu,
        "t_gpu": t_gpu,
        "cpu_aggregate_id": cpu_aggregate["aggregate_id"],
        "gpu_aggregate_id": gpu_aggregate["aggregate_id"],
    }


def build_speedup_curve(points):
    """Sort `points` (each carrying at least `natom` and `speedup`) by `natom`.

    Raises if `points` is empty or `natom` is not unique across points -- a
    crossover curve is one point per tested size, not several repeats
    silently pooled (repeats belong in the aggregate's own `sample_count`).
    """
    points = list(points)
    if not points:
        raise CrossoverError("no points given to build a speedup curve from")
    points = sorted(points, key=lambda p: p["natom"])
    natoms = [p["natom"] for p in points]
    if len(set(natoms)) != len(natoms):
        raise CrossoverError(f"duplicate natom among curve points: {natoms}")
    return points


def _interpolate_crossover_natom(lo, hi, threshold):
    lo_n, lo_s = lo["natom"], lo["speedup"]
    hi_n, hi_s = hi["natom"], hi["speedup"]
    if not (lo_n > 0 and hi_n > 0):
        raise CrossoverError("interpolation requires positive natom at both bracketing points")
    if not (lo_s > 0 and hi_s > 0):
        raise CrossoverError("interpolation requires positive speedup at both bracketing points")
    if hi_s <= lo_s:
        raise CrossoverError(
            f"interpolation requires strictly increasing speedup between neighbouring points "
            f"(natom={lo_n} speedup={lo_s} -> natom={hi_n} speedup={hi_s})"
        )
    log_lo_n, log_hi_n = math.log(lo_n), math.log(hi_n)
    log_lo_s, log_hi_s = math.log(lo_s), math.log(hi_s)
    fraction = (math.log(threshold) - log_lo_s) / (log_hi_s - log_lo_s)
    return math.exp(log_lo_n + fraction * (log_hi_n - log_lo_n))


def find_crossover(curve, threshold):
    """Smallest `natom` at which the speedup curve reaches `threshold`.

    `curve` must already be `build_speedup_curve`'s sorted-by-natom output.
    Returns a dict with `status` in `{below_tested_range, within_tested_range,
    above_tested_range}`, `crossover_natom` (rounded per
    `round_significant` when interpolated; the exact measured value when
    not; `None` when the crossover falls outside the tested range) and
    `interpolated` (bool). Never extrapolates: `below_tested_range` and
    `above_tested_range` carry no `crossover_natom`, only the tested range
    itself, exactly the blueprint's required states.
    """
    curve = list(curve)
    if not curve:
        raise CrossoverError("no points in curve")
    tested_range = [curve[0]["natom"], curve[-1]["natom"]]

    if curve[0]["speedup"] >= threshold:
        return {
            "threshold": threshold, "status": BELOW_TESTED_RANGE,
            "crossover_natom": None, "interpolated": False, "tested_range": tested_range,
        }
    if curve[-1]["speedup"] < threshold:
        return {
            "threshold": threshold, "status": ABOVE_TESTED_RANGE,
            "crossover_natom": None, "interpolated": False, "tested_range": tested_range,
        }

    for i in range(1, len(curve)):
        if curve[i]["speedup"] < threshold:
            continue
        if curve[i]["speedup"] == threshold:
            return {
                "threshold": threshold, "status": WITHIN_TESTED_RANGE,
                "crossover_natom": curve[i]["natom"], "interpolated": False,
                "tested_range": tested_range, "bracket": [curve[i]["natom"], curve[i]["natom"]],
            }
        lo, hi = curve[i - 1], curve[i]
        crossover_natom = _interpolate_crossover_natom(lo, hi, threshold)
        return {
            "threshold": threshold, "status": WITHIN_TESTED_RANGE,
            "crossover_natom": round_significant(crossover_natom), "interpolated": True,
            "tested_range": tested_range, "bracket": [lo["natom"], hi["natom"]],
        }

    raise CrossoverError("unreachable: curve is neither below, above, nor bracketing the threshold")


def find_all_crossovers(curve, thresholds=DEFAULT_CROSSOVER_THRESHOLDS):
    """`find_crossover` for every threshold in `thresholds`, keyed by threshold."""
    curve = build_speedup_curve(curve)
    return {threshold: find_crossover(curve, threshold) for threshold in thresholds}
