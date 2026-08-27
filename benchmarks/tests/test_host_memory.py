"""Tests for the WP-08D host (CPU) memory-limit classification."""

from __future__ import annotations

import pytest

from harness import host_memory


# ---------------------------------------------------------------------------
# estimate_host_memory_required_mib
# ---------------------------------------------------------------------------

def test_estimate_matches_real_b04_dhcpnd_measurement_16x16x16():
    # Real do_meminfo peak: 265,807,837 B at natom=16384, directed=21,921,792
    # (benchmarks/cases/B04_dhcpNd/README.md section D).
    estimate = host_memory.estimate_host_memory_required_mib(16384, 21_921_792)
    assert estimate == pytest.approx(265_807_837 / (1024 ** 2), rel=0.01)


def test_estimate_matches_real_b04_dhcpnd_measurement_20x20x20():
    # Real do_meminfo peak: 517,975,005 B at natom=32000, directed=42,816,000.
    estimate = host_memory.estimate_host_memory_required_mib(32000, 42_816_000)
    assert estimate == pytest.approx(517_975_005 / (1024 ** 2), rel=0.01)


def test_estimate_scales_with_directed_interactions():
    # Both sizes are large enough that FIXED_OVERHEAD_BYTES is negligible,
    # so the estimate is expected to scale ~linearly with interaction count.
    small = host_memory.estimate_host_memory_required_mib(1_000_000, 50_000_000)
    large = host_memory.estimate_host_memory_required_mib(1_000_000_000, 50_000_000_000)
    assert large > small
    assert large == pytest.approx(small * 1000, rel=0.05)


def test_estimate_rejects_non_positive_natom():
    with pytest.raises(host_memory.HostMemoryError):
        host_memory.estimate_host_memory_required_mib(0, 100)


def test_estimate_rejects_negative_directed_interactions():
    with pytest.raises(host_memory.HostMemoryError):
        host_memory.estimate_host_memory_required_mib(1000, -1)


# ---------------------------------------------------------------------------
# classify_host_memory_availability
# ---------------------------------------------------------------------------

def test_classify_available_when_comfortably_under_host_memory():
    result = host_memory.classify_host_memory_availability(
        natom=1_048_576, directed_interactions=1_402_994_688, host_memory_available_mib=50_000,
    )
    assert result["available"] is True
    assert result["classification"] == "available"


def test_classify_unavailable_when_size_exceeds_host_memory():
    result = host_memory.classify_host_memory_availability(
        natom=4_000_000, directed_interactions=5_352_000_000, host_memory_available_mib=50_000,
    )
    assert result["available"] is False
    assert result["classification"] == host_memory.UNAVAILABLE_MEMORY


def test_classify_applies_safety_margin():
    kwargs = dict(natom=1_048_576, directed_interactions=1_402_994_688, host_memory_available_mib=17_000)
    lenient = host_memory.classify_host_memory_availability(safety_margin=1.0, **kwargs)
    strict = host_memory.classify_host_memory_availability(safety_margin=1.5, **kwargs)
    assert lenient["available"] is True
    assert strict["available"] is False


def test_classify_records_requested_size_and_available_memory():
    result = host_memory.classify_host_memory_availability(
        natom=131_072, directed_interactions=175_374_336, host_memory_available_mib=50_000,
    )
    assert result["available_mib"] == 50_000
    assert result["estimated_required_mib"] > 0


def test_classify_requires_host_memory_available():
    with pytest.raises(host_memory.HostMemoryError):
        host_memory.classify_host_memory_availability(
            natom=131_072, directed_interactions=175_374_336, host_memory_available_mib=None,
        )


# ---------------------------------------------------------------------------
# query_host_memory_mib
# ---------------------------------------------------------------------------

def test_query_host_memory_parses_real_proc_meminfo():
    result = host_memory.query_host_memory_mib()
    assert result["total_mib"] > 0
    assert 0 < result["available_mib"] <= result["total_mib"]


def test_query_host_memory_parses_synthetic_meminfo(tmp_path):
    meminfo = tmp_path / "meminfo"
    meminfo.write_text(
        "MemTotal:       65634084 kB\n"
        "MemFree:        17018312 kB\n"
        "MemAvailable:   50424124 kB\n"
        "Buffers:          902276 kB\n"
    )
    result = host_memory.query_host_memory_mib(meminfo)
    assert result["total_mib"] == pytest.approx(65634084 / 1024)
    assert result["available_mib"] == pytest.approx(50424124 / 1024)


def test_query_host_memory_rejects_missing_fields(tmp_path):
    meminfo = tmp_path / "meminfo"
    meminfo.write_text("MemFree:        17018312 kB\n")
    with pytest.raises(host_memory.HostMemoryError):
        host_memory.query_host_memory_mib(meminfo)


# ---------------------------------------------------------------------------
# detect_host_oom
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text", [
    "problem of allocation of array",
    "Fortran runtime error: Insufficient virtual memory",
    "std::bad_alloc",
    "Cannot allocate memory",
])
def test_detect_host_oom_matches_known_signatures(text):
    assert host_memory.detect_host_oom(1, "", text) is True
    assert host_memory.detect_host_oom(1, text, "") is True


@pytest.mark.parametrize("returncode", [-9, 137])
def test_detect_host_oom_matches_sigkill_returncode(returncode):
    assert host_memory.detect_host_oom(returncode, "", "") is True


def test_detect_host_oom_false_on_normal_output():
    assert host_memory.detect_host_oom(0, "Simulation finished\n", "") is False


def test_detect_host_oom_handles_none_output():
    assert host_memory.detect_host_oom(0, None, None) is False
