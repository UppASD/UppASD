"""Tests for the WP-05 OpenMP scaling metrics and CPU baseline selection.

Operates on plain aggregate-shaped dicts -- only the fields
`omp_scaling.py` actually reads matter here, not full schema validity (that
is `test_aggregate.py`'s job).
"""

from __future__ import annotations

import json

import pytest

from harness import omp_scaling


def _aggregate(omp_threads, median, **overrides):
    base = {
        "aggregate_id": f"agg-{omp_threads}",
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": "32x32x32",
        "machine_id": "host1",
        "build_id": "build1",
        "backend": "CPU",
        "requested_precision": "DOUBLE",
        "measurement_profile": "DYNAMICS_ONLY",
        "metric": "steady_step_seconds",
        "omp_threads": omp_threads,
        "median": median,
        "all_samples_numerical_valid": True,
        "all_samples_environment_valid": True,
    }
    base.update(overrides)
    return base


def _family(*thread_median_pairs):
    return [_aggregate(threads, median) for threads, median in thread_median_pairs]


# ---------------------------------------------------------------------------
# compute_omp_speedup_table
# ---------------------------------------------------------------------------

def test_speedup_and_efficiency_are_computed_correctly():
    aggregates = _family((1, 10.0), (2, 5.0), (4, 3.0))
    rows = omp_scaling.compute_omp_speedup_table(aggregates)
    by_threads = {row["omp_threads"]: row for row in rows}

    assert by_threads[1]["speedup"] == pytest.approx(1.0)
    assert by_threads[1]["efficiency"] == pytest.approx(1.0)
    assert by_threads[2]["speedup"] == pytest.approx(2.0)
    assert by_threads[2]["efficiency"] == pytest.approx(1.0)
    assert by_threads[4]["speedup"] == pytest.approx(10.0 / 3.0)
    assert by_threads[4]["efficiency"] == pytest.approx((10.0 / 3.0) / 4)


def test_speedup_table_requires_cpu_1t_reference():
    aggregates = _family((2, 5.0), (4, 3.0))
    with pytest.raises(omp_scaling.OmpScalingError, match="omp_threads=1"):
        omp_scaling.compute_omp_speedup_table(aggregates)


def test_speedup_table_rejects_heterogeneous_family():
    aggregates = _family((1, 10.0), (2, 5.0))
    aggregates[1]["case_id"] = "B02_skyrmion2D"
    with pytest.raises(omp_scaling.OmpScalingError, match="one OpenMP-scaling family"):
        omp_scaling.compute_omp_speedup_table(aggregates)


def test_speedup_table_rejects_duplicate_thread_counts():
    aggregates = _family((1, 10.0), (1, 9.0))
    with pytest.raises(omp_scaling.OmpScalingError, match="duplicate"):
        omp_scaling.compute_omp_speedup_table(aggregates)


def test_speedup_table_rejects_gpu_backend():
    aggregates = _family((1, 10.0))
    aggregates[0]["backend"] = "GPU"
    with pytest.raises(omp_scaling.OmpScalingError, match="CPU aggregates only"):
        omp_scaling.compute_omp_speedup_table(aggregates)


def test_speedup_table_excludes_numerically_invalid_aggregates():
    aggregates = _family((1, 10.0), (2, 5.0))
    aggregates.append(_aggregate(4, 0.1, all_samples_numerical_valid=False))
    rows = omp_scaling.compute_omp_speedup_table(aggregates)
    assert {row["omp_threads"] for row in rows} == {1, 2}


def test_speedup_table_can_require_environment_valid():
    aggregates = _family((1, 10.0), (2, 5.0))
    aggregates[1]["all_samples_environment_valid"] = False
    rows = omp_scaling.compute_omp_speedup_table(aggregates, require_environment_valid=True)
    assert {row["omp_threads"] for row in rows} == {1}


# ---------------------------------------------------------------------------
# select_cpu_1t / determine_p_best / select_cpu_best
# ---------------------------------------------------------------------------

def test_select_cpu_1t_returns_the_one_thread_aggregate():
    aggregates = _family((1, 10.0), (2, 5.0))
    assert omp_scaling.select_cpu_1t(aggregates)["aggregate_id"] == "agg-1"


def test_determine_p_best_does_not_assume_largest_p_wins():
    # p=4 is fastest even though p=8 and p=16 are also measured -- a naive
    # "largest p wins" selection would pick 16 and be wrong.
    aggregates = _family((1, 10.0), (2, 5.0), (4, 2.0), (8, 2.5), (16, 3.0))
    best = omp_scaling.determine_p_best(aggregates)
    assert best["omp_threads"] == 4
    assert omp_scaling.select_cpu_best(aggregates)["omp_threads"] == 4


def test_determine_p_best_ignores_invalid_faster_looking_samples():
    aggregates = _family((1, 10.0), (2, 5.0))
    aggregates.append(_aggregate(4, 0.001, all_samples_numerical_valid=False))
    best = omp_scaling.determine_p_best(aggregates)
    assert best["omp_threads"] == 2


def test_determine_p_best_raises_on_all_invalid():
    aggregates = [_aggregate(1, 10.0, all_samples_numerical_valid=False)]
    with pytest.raises(omp_scaling.OmpScalingError, match="numerically valid"):
        omp_scaling.determine_p_best(aggregates)


# ---------------------------------------------------------------------------
# Baseline pointer persistence
# ---------------------------------------------------------------------------

def test_build_baseline_pointer_shape():
    aggregate = _aggregate(4, 2.0, authoritative=True)
    pointer = omp_scaling.build_baseline_pointer(aggregate, omp_scaling.CPU_BEST_ROLE)
    assert pointer["baseline_role"] == "CPU-BEST"
    assert pointer["omp_threads"] == 4
    assert pointer["median_seconds"] == 2.0
    assert pointer["authoritative"] is True


def test_build_baseline_pointer_rejects_unknown_role():
    aggregate = _aggregate(1, 10.0)
    with pytest.raises(omp_scaling.OmpScalingError):
        omp_scaling.build_baseline_pointer(aggregate, "CPU-MEDIAN")


def test_persist_cpu_baseline_writes_a_stable_named_pointer_file(tmp_path):
    cpu_1t = _aggregate(1, 10.0)
    cpu_best = _aggregate(4, 2.0)

    path_1t = omp_scaling.persist_cpu_baseline(cpu_1t, omp_scaling.CPU_1T_ROLE, tmp_path)
    path_best = omp_scaling.persist_cpu_baseline(cpu_best, omp_scaling.CPU_BEST_ROLE, tmp_path)

    assert path_1t.name == "camp1__B01_bccFe__T0__32x32x32__CPU-1T.json"
    assert path_best.name == "camp1__B01_bccFe__T0__32x32x32__CPU-BEST.json"

    written = json.loads(path_best.read_text())
    assert written["baseline_role"] == "CPU-BEST"
    assert written["omp_threads"] == 4
