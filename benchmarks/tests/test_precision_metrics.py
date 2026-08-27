"""Tests for the WP-06 R_GPU_32_64 precision ratio (section F)."""

from __future__ import annotations

import pytest

from harness import precision_metrics


def _aggregate(requested_precision, median, **overrides):
    base = {
        "aggregate_id": f"agg-{requested_precision}",
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": "32x32x32",
        "machine_id": "host1",
        "build_id": f"build-{requested_precision.lower()}",
        "backend": "GPU",
        "gpu_backend": "CUDA",
        "omp_threads": 1,
        "requested_precision": requested_precision,
        "measurement_profile": "DYNAMICS_ONLY",
        "metric": "steady_step_seconds",
        "median": median,
        "all_samples_numerical_valid": True,
        "all_samples_environment_valid": True,
    }
    base.update(overrides)
    return base


def _family(single_median, double_median):
    return [_aggregate("SINGLE", single_median), _aggregate("DOUBLE", double_median)]


# ---------------------------------------------------------------------------
# compute_r_gpu_32_64
# ---------------------------------------------------------------------------

def test_ratio_is_double_over_single():
    aggregates = _family(single_median=2.0, double_median=5.0)
    result = precision_metrics.compute_r_gpu_32_64(aggregates)
    assert result["r_gpu_32_64"] == pytest.approx(2.5)
    assert result["t_single"] == pytest.approx(2.0)
    assert result["t_double"] == pytest.approx(5.0)
    assert result["single_aggregate_id"] == "agg-SINGLE"
    assert result["double_aggregate_id"] == "agg-DOUBLE"


def test_ratio_requires_both_precisions():
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="missing"):
        precision_metrics.compute_r_gpu_32_64([_aggregate("SINGLE", 2.0)])
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="missing"):
        precision_metrics.compute_r_gpu_32_64([_aggregate("DOUBLE", 5.0)])


def test_ratio_rejects_mixed_aggregate():
    aggregates = _family(2.0, 5.0) + [_aggregate("MIXED", 0.1)]
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="SINGLE/DOUBLE"):
        precision_metrics.compute_r_gpu_32_64(aggregates)


def test_ratio_requires_gpu_backend():
    aggregates = _family(2.0, 5.0)
    for aggregate in aggregates:
        aggregate["backend"] = "CPU"
        aggregate["gpu_backend"] = None
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="GPU aggregates only"):
        precision_metrics.compute_r_gpu_32_64(aggregates)


def test_ratio_rejects_heterogeneous_family():
    aggregates = _family(2.0, 5.0)
    aggregates[1]["size_id"] = "64x64x64"
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="one GPU precision family"):
        precision_metrics.compute_r_gpu_32_64(aggregates)


def test_ratio_excludes_numerically_invalid_aggregates():
    aggregates = _family(2.0, 5.0)
    aggregates.append(_aggregate("SINGLE", 0.001, aggregate_id="agg-bad", all_samples_numerical_valid=False))
    result = precision_metrics.compute_r_gpu_32_64(aggregates)
    assert result["t_single"] == pytest.approx(2.0)


def test_ratio_can_require_environment_valid():
    aggregates = _family(2.0, 5.0)
    aggregates[0]["all_samples_environment_valid"] = False
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="missing"):
        precision_metrics.compute_r_gpu_32_64(aggregates, require_environment_valid=True)


def test_ratio_rejects_duplicate_precision():
    aggregates = [_aggregate("SINGLE", 2.0), _aggregate("SINGLE", 3.0, aggregate_id="agg-SINGLE-2")]
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="duplicate"):
        precision_metrics.compute_r_gpu_32_64(aggregates)


# ---------------------------------------------------------------------------
# compute_r_cpu_32_64 -- deliberately unsupported
# ---------------------------------------------------------------------------

def test_cpu_ratio_is_deliberately_unimplemented():
    with pytest.raises(precision_metrics.PrecisionMetricsError, match="no analogous CPU precision ratio"):
        precision_metrics.compute_r_cpu_32_64([])
