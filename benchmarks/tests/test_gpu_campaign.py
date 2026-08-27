"""Tests for the WP-06 GPU campaign record builders (sections C-E)."""

from __future__ import annotations

import json

import pytest

from harness import gpu_campaign
from harness import gpu_memory
from harness import omp_sweep


class _FakeCase:
    def __init__(self, case_id="B01_bccFe", workload_class="NEIGHBOR_LIST"):
        self.id = case_id
        self.record_workload_class = workload_class


def _static_context(**overrides):
    base = {
        "omp_threads": 1,
        "git_commit": "a" * 40,
        "git_dirty": False,
        "git_dirty_files": [],
        "compiler": "nvcc",
        "compiler_version": "12.4.131",
        "compile_flags": ["-O3"],
        "cmake_options": {"USE_CUDA": True},
        "binary_checksum": "b" * 64,
        "build_id": "build-cuda-single-0123456",
        "cpu_model": "AMD Ryzen 9 5950X 16-Core Processor",
        "cpu_sockets": 1,
        "cpu_physical_cores": 16,
        "cpu_threads": 32,
        "numa_nodes": 1,
        "gpu_runtime": "CUDA 12.4",
    }
    base.update(overrides)
    return base


def _workload_metadata(**overrides):
    base = {
        "natom": 65536,
        "directed_interactions": 3_801_088,
        "mean_neighbors": 58.0,
        "median_neighbors": 58.0,
        "max_neighbors": 58,
        "fft_grid": None,
        "fft_grid_padded": None,
        "fft_grid_points": None,
        "dipole_backend": None,
    }
    base.update(overrides)
    return base


# ---------------------------------------------------------------------------
# resolve_gpu_host_env
# ---------------------------------------------------------------------------

def test_resolve_gpu_host_env_defaults_to_one_thread():
    env = gpu_campaign.resolve_gpu_host_env(base_env={})
    assert env["OMP_NUM_THREADS"] == "1"
    assert env["OMP_PLACES"] == "cores"
    assert env["OMP_DYNAMIC"] == "FALSE"


def test_resolve_gpu_host_env_matches_omp_sweep_semantics():
    env = gpu_campaign.resolve_gpu_host_env(omp_threads=1, base_env={"PATH": "/usr/bin"})
    expected = omp_sweep.build_omp_env(1, proc_bind="close", places="cores", base_env={"PATH": "/usr/bin"})
    assert env == expected


# ---------------------------------------------------------------------------
# build_unsupported_precision_record
# ---------------------------------------------------------------------------

def test_unsupported_precision_record_is_schema_valid_and_skipped():
    record = gpu_campaign.build_unsupported_precision_record(
        run_id="lean01-B01-T0-32x32x32-cuda-mixed-s000", campaign_id="lean01", case=_FakeCase(),
        variant_id="T0", size_id="32x32x32", sample_index=0, machine_id="host1",
        gpu_backend="CUDA", workload_metadata=_workload_metadata(), context=_static_context(),
        temperature=0.0, timestep=1e-16, nstep=10000, measurement_profile="DYNAMICS_ONLY",
    )
    assert record["run_status"] == "SKIPPED"
    assert record["precision_support_state"] == "unsupported"
    assert record["numerical_valid"] is False
    assert record["environment_valid"] is False
    assert record["comparison_precision_class"] is None
    assert record["gpu_model"] is None
    assert record["requested_precision"] == "MIXED"
    assert record["natom"] == 65536
    assert "MIXED" in record["notes"]


def test_unsupported_precision_record_rejects_supported_precision():
    with pytest.raises(gpu_campaign.GpuCampaignError):
        gpu_campaign.build_unsupported_precision_record(
            run_id="run1", campaign_id="lean01", case=_FakeCase(), variant_id="T0", size_id="32x32x32",
            sample_index=0, machine_id="host1", gpu_backend="CUDA", workload_metadata=_workload_metadata(),
            context=_static_context(), temperature=0.0, timestep=1e-16, nstep=10000,
            measurement_profile="DYNAMICS_ONLY", requested_precision="DOUBLE",
        )


# ---------------------------------------------------------------------------
# build_unavailable_memory_record
# ---------------------------------------------------------------------------

def _unavailable_classification():
    return gpu_memory.classify_gpu_memory_availability(
        natom=500_000_000, directed_interactions=30_000_000_000, effective_gpu_precision="DOUBLE",
        device_memory_total_mib=16376, device_memory_used_mib=100,
    )


def test_unavailable_memory_record_is_schema_valid_and_skipped():
    device = {"model": "NVIDIA RTX A4000", "uuid": "GPU-abc", "driver_version": "550.54.14", "index": 0}
    record = gpu_campaign.build_unavailable_memory_record(
        run_id="lean01-B04-huge-cuda-double-s000", campaign_id="lean01", case=_FakeCase(case_id="B04_dhcpNd"),
        variant_id="T0", size_id="256x256x256", sample_index=0, machine_id="host1",
        gpu_backend="CUDA", requested_precision="DOUBLE",
        workload_metadata=_workload_metadata(
            natom=500_000_000, directed_interactions=30_000_000_000,
            mean_neighbors=60.0, median_neighbors=60.0, max_neighbors=60,
        ),
        context=_static_context(), memory_classification=_unavailable_classification(), device=device,
        temperature=0.0, timestep=1e-16, nstep=10000, measurement_profile="DYNAMICS_ONLY",
    )
    assert record["run_status"] == "SKIPPED"
    assert record["precision_support_state"] == "supported"
    assert record["numerical_valid"] is False
    assert record["comparison_precision_class"] is None
    assert record["gpu_model"] == "NVIDIA RTX A4000"
    assert record["effective_gpu_precision"] == "DOUBLE"
    assert record["effective_cpu_precision"] == "DOUBLE"
    assert "unavailable_memory" in record["notes"]
    assert "500000000" in record["notes"]


def test_unavailable_memory_record_rejects_an_available_classification():
    device = {"model": "NVIDIA RTX A4000", "uuid": "GPU-abc", "index": 0}
    available = gpu_memory.classify_gpu_memory_availability(
        natom=65536, directed_interactions=3_801_088, effective_gpu_precision="SINGLE",
        device_memory_total_mib=16376, device_memory_used_mib=100,
    )
    with pytest.raises(gpu_campaign.GpuCampaignError):
        gpu_campaign.build_unavailable_memory_record(
            run_id="run1", campaign_id="lean01", case=_FakeCase(), variant_id="T0", size_id="32x32x32",
            sample_index=0, machine_id="host1", gpu_backend="CUDA", requested_precision="SINGLE",
            workload_metadata=_workload_metadata(), context=_static_context(),
            memory_classification=available, device=device,
            temperature=0.0, timestep=1e-16, nstep=10000, measurement_profile="DYNAMICS_ONLY",
        )


# ---------------------------------------------------------------------------
# write_record
# ---------------------------------------------------------------------------

def test_write_record_uses_run_id_filename(tmp_path):
    record = gpu_campaign.build_unsupported_precision_record(
        run_id="my-run-id", campaign_id="lean01", case=_FakeCase(), variant_id="T0", size_id="32x32x32",
        sample_index=0, machine_id="host1", gpu_backend="CUDA", workload_metadata=_workload_metadata(),
        context=_static_context(), temperature=0.0, timestep=1e-16, nstep=10000,
        measurement_profile="DYNAMICS_ONLY",
    )
    path = gpu_campaign.write_record(record, tmp_path)
    assert path.name == "my-run-id.json"
    assert json.loads(path.read_text())["run_id"] == "my-run-id"
