"""Tests for the WP-06 precision audit (section A/B)."""

from __future__ import annotations

import pytest

from harness import precision_audit


# ---------------------------------------------------------------------------
# effective_cpu_precision
# ---------------------------------------------------------------------------

def test_cpu_precision_is_always_double_regardless_of_request():
    assert precision_audit.effective_cpu_precision("DOUBLE") == "DOUBLE"
    assert precision_audit.effective_cpu_precision("SINGLE") == "DOUBLE"


def test_cpu_precision_unknown_for_mixed():
    assert precision_audit.effective_cpu_precision("MIXED") == "UNKNOWN"


def test_cpu_precision_rejects_bogus_request():
    with pytest.raises(ValueError):
        precision_audit.effective_cpu_precision("BOGUS")


# ---------------------------------------------------------------------------
# effective_gpu_precision
# ---------------------------------------------------------------------------

def test_gpu_precision_none_for_cpu_only_run():
    assert precision_audit.effective_gpu_precision("DOUBLE", None) is None


def test_gpu_precision_mirrors_request_for_cuda():
    assert precision_audit.effective_gpu_precision("SINGLE", "CUDA") == "SINGLE"
    assert precision_audit.effective_gpu_precision("DOUBLE", "CUDA") == "DOUBLE"


def test_gpu_precision_unknown_for_unaudited_hip():
    assert precision_audit.effective_gpu_precision("SINGLE", "HIP") == "UNKNOWN"
    assert precision_audit.effective_gpu_precision("DOUBLE", "HIP") == "UNKNOWN"


def test_gpu_precision_unknown_for_mixed():
    assert precision_audit.effective_gpu_precision("MIXED", "CUDA") == "UNKNOWN"


def test_gpu_precision_rejects_bogus_backend():
    with pytest.raises(ValueError):
        precision_audit.effective_gpu_precision("SINGLE", "OPENCL")


# ---------------------------------------------------------------------------
# classify_comparison_precision_class -- mirrors the four shipped fixtures
# ---------------------------------------------------------------------------

def test_classify_cpu_double_completed_is_production_configuration():
    # cpu_bccfe_dynamics_only.json
    result = precision_audit.classify_comparison_precision_class(
        backend="CPU", gpu_backend=None, precision_support_state="supported",
        effective_cpu_precision="DOUBLE", effective_gpu_precision=None, run_status="COMPLETED",
    )
    assert result == "PRODUCTION_CONFIGURATION"


def test_classify_gpu_single_vs_cpu_double_is_production_configuration():
    # cuda_bccfe_single_precision.json
    result = precision_audit.classify_comparison_precision_class(
        backend="GPU", gpu_backend="CUDA", precision_support_state="supported",
        effective_cpu_precision="DOUBLE", effective_gpu_precision="SINGLE", run_status="COMPLETED",
    )
    assert result == "PRODUCTION_CONFIGURATION"


def test_classify_gpu_double_matches_cpu_double():
    # cuda_dipole_fft_plus_neighbours.json
    result = precision_audit.classify_comparison_precision_class(
        backend="GPU", gpu_backend="CUDA", precision_support_state="supported",
        effective_cpu_precision="DOUBLE", effective_gpu_precision="DOUBLE", run_status="COMPLETED",
    )
    assert result == "PRECISION_MATCHED"


def test_classify_mixed_unsupported_is_always_null():
    # cuda_mixed_precision_unsupported.json
    result = precision_audit.classify_comparison_precision_class(
        backend="GPU", gpu_backend="CUDA", precision_support_state="unsupported",
        effective_cpu_precision="UNKNOWN", effective_gpu_precision="UNKNOWN", run_status="SKIPPED",
    )
    assert result is None


def test_classify_unaudited_skipped_hip_is_null():
    # hip_unmeasured_backend.json -- unaudited AND never executed
    result = precision_audit.classify_comparison_precision_class(
        backend="GPU", gpu_backend="HIP", precision_support_state="unaudited",
        effective_cpu_precision="DOUBLE", effective_gpu_precision="UNKNOWN", run_status="SKIPPED",
    )
    assert result is None


def test_classify_unaudited_completed_reports_unaudited():
    result = precision_audit.classify_comparison_precision_class(
        backend="CPU", gpu_backend=None, precision_support_state="unaudited",
        effective_cpu_precision="UNKNOWN", effective_gpu_precision=None, run_status="COMPLETED",
    )
    assert result == "UNAUDITED"


def test_classify_supported_but_not_completed_is_null():
    result = precision_audit.classify_comparison_precision_class(
        backend="GPU", gpu_backend="CUDA", precision_support_state="supported",
        effective_cpu_precision="DOUBLE", effective_gpu_precision="SINGLE", run_status="SKIPPED",
    )
    assert result is None


def test_classify_rejects_unknown_precision_support_state():
    with pytest.raises(ValueError):
        precision_audit.classify_comparison_precision_class(
            backend="CPU", gpu_backend=None, precision_support_state="bogus",
            effective_cpu_precision="DOUBLE", effective_gpu_precision=None, run_status="COMPLETED",
        )


def test_classify_rejects_unknown_backend():
    with pytest.raises(ValueError):
        precision_audit.classify_comparison_precision_class(
            backend="TPU", gpu_backend=None, precision_support_state="supported",
            effective_cpu_precision="DOUBLE", effective_gpu_precision=None, run_status="COMPLETED",
        )


# ---------------------------------------------------------------------------
# COMPONENT_AUDIT stays internally consistent with the functions above
# ---------------------------------------------------------------------------

def test_component_audit_table_covers_the_required_eight_components():
    assert len(precision_audit.COMPONENT_AUDIT) == 8
    for row in precision_audit.COMPONENT_AUDIT:
        assert row["citation"]
        assert row["double_build"]
        assert row["single_build"]


def test_component_audit_cpu_rows_agree_with_effective_cpu_precision():
    cpu_rows = [row for row in precision_audit.COMPONENT_AUDIT if row["component"].startswith("CPU")]
    assert len(cpu_rows) == 3
    for row in cpu_rows:
        assert row["double_build"] == "DOUBLE"
        assert row["single_build"] == "DOUBLE"  # unaffected by UPPASD_PRECISION=SINGLE
