"""Tests for the WP-07 GPU/CPU crossover metrics (sections A-B)."""

from __future__ import annotations

import pytest

from analysis import crossover


def _aggregate(backend, median, size_id="32x32x32", **overrides):
    base = {
        "aggregate_id": f"agg-{backend}-{size_id}",
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": size_id,
        "machine_id": "host1",
        "measurement_profile": "DYNAMICS_ONLY",
        "metric": "steady_step_seconds",
        "backend": backend,
        "median": median,
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# compute_gpu_speedup
# ---------------------------------------------------------------------------

def test_speedup_is_cpu_over_gpu():
    cpu = _aggregate("CPU", 10.0)
    gpu = _aggregate("GPU", 2.0)
    result = crossover.compute_gpu_speedup(cpu, gpu, label="S_GPU_BESTCPU")
    assert result["speedup"] == pytest.approx(5.0)
    assert result["t_cpu"] == pytest.approx(10.0)
    assert result["t_gpu"] == pytest.approx(2.0)
    assert result["label"] == "S_GPU_BESTCPU"


def test_speedup_requires_cpu_and_gpu_backends():
    with pytest.raises(crossover.CrossoverError, match="cpu_aggregate backend must be CPU"):
        crossover.compute_gpu_speedup(_aggregate("GPU", 10.0), _aggregate("GPU", 2.0), label="x")
    with pytest.raises(crossover.CrossoverError, match="gpu_aggregate backend must be GPU"):
        crossover.compute_gpu_speedup(_aggregate("CPU", 10.0), _aggregate("CPU", 2.0), label="x")


def test_speedup_requires_matching_cell():
    cpu = _aggregate("CPU", 10.0)
    gpu = _aggregate("GPU", 2.0, size_id="64x64x64")
    with pytest.raises(crossover.CrossoverError, match="do not share one cell"):
        crossover.compute_gpu_speedup(cpu, gpu, label="x")


def test_speedup_requires_positive_medians():
    with pytest.raises(crossover.CrossoverError, match="positive"):
        crossover.compute_gpu_speedup(_aggregate("CPU", 0.0), _aggregate("GPU", 2.0), label="x")


# ---------------------------------------------------------------------------
# build_speedup_curve
# ---------------------------------------------------------------------------

def test_build_speedup_curve_sorts_by_natom():
    points = [
        {"natom": 65536, "speedup": 3.0},
        {"natom": 4096, "speedup": 0.5},
        {"natom": 16384, "speedup": 1.5},
    ]
    curve = crossover.build_speedup_curve(points)
    assert [p["natom"] for p in curve] == [4096, 16384, 65536]


def test_build_speedup_curve_rejects_empty():
    with pytest.raises(crossover.CrossoverError, match="no points"):
        crossover.build_speedup_curve([])


def test_build_speedup_curve_rejects_duplicate_natom():
    points = [{"natom": 4096, "speedup": 1.0}, {"natom": 4096, "speedup": 2.0}]
    with pytest.raises(crossover.CrossoverError, match="duplicate natom"):
        crossover.build_speedup_curve(points)


# ---------------------------------------------------------------------------
# find_crossover
# ---------------------------------------------------------------------------

def _curve(*natom_speedup_pairs):
    return crossover.build_speedup_curve(
        [{"natom": n, "speedup": s} for n, s in natom_speedup_pairs]
    )


def test_crossover_below_tested_range_when_smallest_point_already_exceeds_threshold():
    curve = _curve((4096, 2.0), (16384, 4.0), (65536, 8.0))
    result = crossover.find_crossover(curve, 1.0)
    assert result["status"] == crossover.BELOW_TESTED_RANGE
    assert result["crossover_natom"] is None
    assert result["interpolated"] is False
    assert result["tested_range"] == [4096, 65536]


def test_crossover_above_tested_range_when_largest_point_never_reaches_threshold():
    curve = _curve((4096, 0.1), (16384, 0.3), (65536, 0.8))
    result = crossover.find_crossover(curve, 1.0)
    assert result["status"] == crossover.ABOVE_TESTED_RANGE
    assert result["crossover_natom"] is None
    assert result["interpolated"] is False


def test_crossover_interpolates_between_bracketing_points():
    # speedup 0.5 at 1000, 2.0 at 100000: log-log interpolation for threshold 1.0.
    curve = _curve((1000, 0.5), (100000, 2.0))
    result = crossover.find_crossover(curve, 1.0)
    assert result["status"] == crossover.WITHIN_TESTED_RANGE
    assert result["interpolated"] is True
    assert result["bracket"] == [1000, 100000]
    # geometric-mean-like midpoint in log-log space; well inside the bracket.
    assert 1000 < result["crossover_natom"] < 100000


def test_crossover_reports_exact_measured_point_without_interpolation():
    curve = _curve((4096, 0.5), (16384, 1.0), (65536, 3.0))
    result = crossover.find_crossover(curve, 1.0)
    assert result["status"] == crossover.WITHIN_TESTED_RANGE
    assert result["interpolated"] is False
    assert result["crossover_natom"] == 16384


def test_crossover_display_precision_is_rounded():
    curve = _curve((1000, 0.9999), (123457, 2.0))
    result = crossover.find_crossover(curve, 1.0)
    assert result["interpolated"] is True
    # rounded to 2 significant figures -- no long decimal tail.
    text = repr(result["crossover_natom"])
    significant_digits = text.replace(".", "").replace("-", "").lstrip("0")
    assert len(significant_digits.rstrip("0") or "0") <= 2 or significant_digits == "0"


def test_find_all_crossovers_covers_default_thresholds():
    curve = [
        {"natom": 1000, "speedup": 0.2},
        {"natom": 10000, "speedup": 0.8},
        {"natom": 100000, "speedup": 3.0},
        {"natom": 1000000, "speedup": 6.0},
    ]
    results = crossover.find_all_crossovers(curve)
    assert set(results) == set(crossover.DEFAULT_CROSSOVER_THRESHOLDS)
    assert results[1.0]["status"] == crossover.WITHIN_TESTED_RANGE
    assert results[5.0]["status"] == crossover.WITHIN_TESTED_RANGE


def test_round_significant_handles_zero():
    assert crossover.round_significant(0.0) == 0.0
