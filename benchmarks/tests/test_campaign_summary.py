"""Tests for the WP-07 machine-readable campaign summary (section G)."""

from __future__ import annotations

import csv
import json

import pytest

from analysis import campaign_summary


def _aggregate(backend, median, *, size_id, aggregate_id, omp_threads=1, mad=0.0,
               authoritative=True, quality_flags_union=None, **overrides):
    base = {
        "aggregate_id": aggregate_id,
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": size_id,
        "machine_id": "host1",
        "measurement_profile": "DYNAMICS_ONLY",
        "metric": "steady_step_seconds",
        "backend": backend,
        "omp_threads": omp_threads,
        "median": median,
        "mad": mad,
        "authoritative": authoritative,
        "quality_flags_union": quality_flags_union or [],
        "all_samples_numerical_valid": True,
        "all_samples_environment_valid": True,
    }
    base.update(overrides)
    return base


def _entry(size_id, natom, requested_precision, *, t_gpu, t_cpu_best, t_cpu_1t,
           directed_interactions=None, fft_grid_points=None,
           workload_class="NEIGHBOR_LIST", cpu_best_threads=8):
    return {
        "size_id": size_id,
        "natom": natom,
        "directed_interactions": directed_interactions,
        "fft_grid_points": fft_grid_points,
        "workload_class": workload_class,
        "requested_precision": requested_precision,
        "gpu_aggregate": _aggregate(
            "GPU", t_gpu, size_id=size_id, aggregate_id=f"gpu-{requested_precision}-{size_id}",
            gpu_backend="CUDA", requested_precision=requested_precision,
        ),
        "cpu_best_aggregate": _aggregate(
            "CPU", t_cpu_best, size_id=size_id, aggregate_id=f"cpubest-{size_id}",
            omp_threads=cpu_best_threads, requested_precision="DOUBLE",
        ),
        "cpu_1t_aggregate": _aggregate(
            "CPU", t_cpu_1t, size_id=size_id, aggregate_id=f"cpu1t-{size_id}",
            omp_threads=1, requested_precision="DOUBLE",
        ),
    }


def _entries():
    return [
        _entry("4096", 4096, "DOUBLE", t_gpu=1.0, t_cpu_best=0.5, t_cpu_1t=2.0,
               directed_interactions=200000),
        _entry("32768", 32768, "DOUBLE", t_gpu=2.0, t_cpu_best=4.0, t_cpu_1t=16.0,
               directed_interactions=1600000),
        _entry("4096", 4096, "SINGLE", t_gpu=0.6, t_cpu_best=0.5, t_cpu_1t=2.0,
               directed_interactions=200000),
        _entry("32768", 32768, "SINGLE", t_gpu=1.0, t_cpu_best=4.0, t_cpu_1t=16.0,
               directed_interactions=1600000),
    ]


# ---------------------------------------------------------------------------
# build_case_summary
# ---------------------------------------------------------------------------

def test_build_case_summary_shape():
    summary = campaign_summary.build_case_summary(_entries())
    assert summary["identity"] == {
        "campaign_id": "camp1", "case_id": "B01_bccFe", "variant_id": "T0", "machine_id": "host1",
    }
    assert set(summary["by_precision"]) == {"SINGLE", "DOUBLE"}

    double = summary["by_precision"]["DOUBLE"]
    assert double["measured_range"] == {"min_natom": 4096, "max_natom": 32768}
    rows = double["rows"]
    assert [row["natom"] for row in rows] == [4096, 32768]
    assert rows[0]["s_gpu_bestcpu"] == pytest.approx(0.5)
    assert rows[1]["s_gpu_bestcpu"] == pytest.approx(2.0)
    assert rows[0]["s_gpu_1t"] == pytest.approx(2.0)


def test_build_case_summary_computes_crossover_per_precision():
    summary = campaign_summary.build_case_summary(_entries())
    double_1x = summary["by_precision"]["DOUBLE"]["crossover"]["1.0"]
    assert double_1x["status"] in {"within_tested_range", "below_tested_range", "above_tested_range"}


def test_build_case_summary_computes_precision_penalty():
    summary = campaign_summary.build_case_summary(_entries())
    penalties = {p["size_id"]: p for p in summary["precision_penalty"]}
    assert set(penalties) == {"4096", "32768"}
    # DOUBLE at 4096 is 1.0s, SINGLE is 0.6s -> R = 1.0/0.6
    assert penalties["4096"]["r_gpu_32_64"] == pytest.approx(1.0 / 0.6)


def test_build_case_summary_rejects_empty():
    with pytest.raises(campaign_summary.CampaignSummaryError, match="no entries"):
        campaign_summary.build_case_summary([])


def test_build_case_summary_rejects_mixed_case_identity():
    entries = _entries()
    entries[0]["gpu_aggregate"]["case_id"] = "B04_dhcpNd"
    with pytest.raises(campaign_summary.CampaignSummaryError, match="do not share one"):
        campaign_summary.build_case_summary(entries)


def test_build_case_summary_quality_status_reflects_flags():
    entries = _entries()
    entries[0]["gpu_aggregate"]["authoritative"] = False
    entries[0]["gpu_aggregate"]["quality_flags_union"] = ["gpu_busy"]
    summary = campaign_summary.build_case_summary(entries)
    status = summary["by_precision"]["DOUBLE"]["quality_status"]
    assert status["authoritative"] is False
    assert "gpu_busy" in status["quality_flags_union"]


# ---------------------------------------------------------------------------
# write_summary_json / write_summary_csv / write_crossover_csv / write_precision_penalty_csv
# ---------------------------------------------------------------------------

def test_write_summary_json_round_trips(tmp_path):
    summary = campaign_summary.build_case_summary(_entries())
    path = campaign_summary.write_summary_json(summary, tmp_path / "out" / "summary.json")
    assert path.is_file()
    assert json.loads(path.read_text()) == summary


def test_write_summary_csv_has_one_row_per_measured_cell(tmp_path):
    summary = campaign_summary.build_case_summary(_entries())
    path = campaign_summary.write_summary_csv(summary, tmp_path / "summary.csv")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert {row["requested_precision"] for row in rows} == {"SINGLE", "DOUBLE"}


def test_write_crossover_csv_has_one_row_per_precision_threshold(tmp_path):
    summary = campaign_summary.build_case_summary(_entries())
    path = campaign_summary.write_crossover_csv(summary, tmp_path / "crossover.csv")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2 * 4  # 2 precisions x 4 default thresholds


def test_write_precision_penalty_csv(tmp_path):
    summary = campaign_summary.build_case_summary(_entries())
    path = campaign_summary.write_precision_penalty_csv(summary, tmp_path / "penalty.csv")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
