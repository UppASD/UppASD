"""Tests for the WP-06 GPU memory-limit classification (section E)."""

from __future__ import annotations

import pytest

from harness import gpu_memory


# ---------------------------------------------------------------------------
# estimate_gpu_memory_required_mib
# ---------------------------------------------------------------------------

def test_estimate_scales_with_natom():
    small = gpu_memory.estimate_gpu_memory_required_mib(1000, 50000, "SINGLE")
    large = gpu_memory.estimate_gpu_memory_required_mib(1_000_000, 50_000_000, "SINGLE")
    assert large > small
    assert large == pytest.approx(small * 1000, rel=0.05)


def test_estimate_single_is_smaller_than_double():
    single = gpu_memory.estimate_gpu_memory_required_mib(65536, 3_801_088, "SINGLE")
    double = gpu_memory.estimate_gpu_memory_required_mib(65536, 3_801_088, "DOUBLE")
    assert single < double


def test_estimate_handles_null_directed_interactions_for_fft_workloads():
    estimate = gpu_memory.estimate_gpu_memory_required_mib(65536, None, "DOUBLE")
    assert estimate > 0


def test_estimate_rejects_non_positive_natom():
    with pytest.raises(gpu_memory.GpuMemoryError):
        gpu_memory.estimate_gpu_memory_required_mib(0, 100, "DOUBLE")


def test_estimate_rejects_unaudited_precision():
    with pytest.raises(gpu_memory.GpuMemoryError):
        gpu_memory.estimate_gpu_memory_required_mib(1000, 100, "UNKNOWN")


# ---------------------------------------------------------------------------
# classify_gpu_memory_availability
# ---------------------------------------------------------------------------

def test_classify_available_when_comfortably_under_device_memory():
    result = gpu_memory.classify_gpu_memory_availability(
        natom=65536, directed_interactions=3_801_088, effective_gpu_precision="SINGLE",
        device_memory_total_mib=16376, device_memory_used_mib=100,
    )
    assert result["available"] is True
    assert result["classification"] == "available"


def test_classify_unavailable_when_size_exceeds_device_memory():
    result = gpu_memory.classify_gpu_memory_availability(
        natom=500_000_000, directed_interactions=30_000_000_000, effective_gpu_precision="DOUBLE",
        device_memory_total_mib=16376, device_memory_used_mib=100,
    )
    assert result["available"] is False
    assert result["classification"] == gpu_memory.UNAVAILABLE_MEMORY


def test_classify_accounts_for_already_used_memory():
    kwargs = dict(
        natom=65536, directed_interactions=3_801_088, effective_gpu_precision="SINGLE",
        device_memory_total_mib=16376,
    )
    idle = gpu_memory.classify_gpu_memory_availability(device_memory_used_mib=0, **kwargs)
    busy = gpu_memory.classify_gpu_memory_availability(device_memory_used_mib=16350, **kwargs)
    assert idle["available"] is True
    assert busy["available"] is False


def test_classify_records_requested_size_and_available_memory():
    result = gpu_memory.classify_gpu_memory_availability(
        natom=65536, directed_interactions=3_801_088, effective_gpu_precision="SINGLE",
        device_memory_total_mib=16376, device_memory_used_mib=376,
    )
    assert result["available_mib"] == pytest.approx(16000)
    assert result["estimated_required_mib"] > 0


def test_classify_requires_device_memory_total():
    with pytest.raises(gpu_memory.GpuMemoryError):
        gpu_memory.classify_gpu_memory_availability(
            natom=65536, directed_interactions=3_801_088, effective_gpu_precision="SINGLE",
            device_memory_total_mib=None,
        )


# ---------------------------------------------------------------------------
# detect_gpu_oom
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "CUDA error: out of memory",
    "terminate called after throwing an instance of 'thrust::system::detail::bad_alloc' cudaErrorMemoryAllocation",
    "CUFFT_ALLOC_FAILED",
    "std::bad_alloc: cannot allocate memory",
])
def test_detect_gpu_oom_matches_known_signatures(text):
    assert gpu_memory.detect_gpu_oom("", text) is True
    assert gpu_memory.detect_gpu_oom(text, "") is True


def test_detect_gpu_oom_false_on_normal_output():
    assert gpu_memory.detect_gpu_oom("Simulation finished\n", "") is False


def test_detect_gpu_oom_handles_none_inputs():
    assert gpu_memory.detect_gpu_oom(None, None) is False
