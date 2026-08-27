"""Tests for the WP-05 OpenMP scaling report (machine-readable summary + plots).

Uses the Agg backend (set in the module itself) so this runs headless; plots
are checked for existence/non-emptiness, not pixel content -- rendering
correctness is `matplotlib`'s job, not this harness's.
"""

from __future__ import annotations

import json

import pytest

from analysis import omp_scaling_report


def _aggregate(size_id, omp_threads, median):
    return {
        "aggregate_id": f"agg-{size_id}-{omp_threads}",
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": size_id,
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


def _size_families():
    return [
        {
            "size_id": "16x16x16",
            "natom": 4096,
            "aggregates": [
                _aggregate("16x16x16", 1, 10.0),
                _aggregate("16x16x16", 2, 5.5),
                _aggregate("16x16x16", 4, 3.0),
            ],
        },
        {
            "size_id": "32x32x32",
            "natom": 32768,
            "aggregates": [
                _aggregate("32x32x32", 1, 80.0),
                _aggregate("32x32x32", 2, 42.0),
                _aggregate("32x32x32", 4, 21.0),
                _aggregate("32x32x32", 8, 22.0),  # p_best is 4, not the largest p
            ],
        },
    ]


def test_build_scaling_summary_shape_and_p_best():
    summary = omp_scaling_report.build_scaling_summary(_size_families())
    by_size = {size["size_id"]: size for size in summary["sizes"]}

    assert by_size["16x16x16"]["p_best"] == 4
    assert by_size["32x32x32"]["p_best"] == 4  # not 8, despite 8 being measured

    assert [row["natom"] for row in summary["p_best_vs_natom"]] == [4096, 32768]
    assert [row["p_best"] for row in summary["p_best_vs_natom"]] == [4, 4]


def test_build_scaling_summary_sorts_p_best_by_natom_regardless_of_input_order():
    families = list(reversed(_size_families()))
    summary = omp_scaling_report.build_scaling_summary(families)
    assert [row["natom"] for row in summary["p_best_vs_natom"]] == [4096, 32768]


def test_build_scaling_summary_rejects_empty_input():
    with pytest.raises(omp_scaling_report.OmpScalingReportError):
        omp_scaling_report.build_scaling_summary([])


def test_write_scaling_summary_writes_valid_json(tmp_path):
    summary = omp_scaling_report.build_scaling_summary(_size_families())
    path = omp_scaling_report.write_scaling_summary(summary, tmp_path / "out" / "summary.json")
    assert path.is_file()
    reloaded = json.loads(path.read_text())
    assert reloaded == summary


def test_plot_omp_scaling_writes_one_png_per_size_plus_p_best(tmp_path):
    summary = omp_scaling_report.build_scaling_summary(_size_families())
    written = omp_scaling_report.plot_omp_scaling(summary, tmp_path / "plots")

    names = {path.name for path in written}
    assert names == {
        "omp_scaling__16x16x16.png",
        "omp_scaling__32x32x32.png",
        "p_best_vs_natom.png",
    }
    for path in written:
        assert path.is_file()
        assert path.stat().st_size > 0
