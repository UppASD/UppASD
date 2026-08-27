"""WP-04 section G: provenance capture tests.

Real hardware (`nvidia-smi`, `/proc/cpuinfo`, `git`) is not assumed to be
available or deterministic in CI, so every external call is exercised
through a mock `run` callable or a synthetic filesystem tree rather than the
real toolchain -- this is what the checklist calls "mock parsers".
"""

from __future__ import annotations

import json

import pytest

from harness import cpu_provenance
from harness import gpu_provenance
from harness import provenance
from harness import source_provenance


class FakeCompleted:
    def __init__(self, returncode=0, stdout=""):
        self.returncode = returncode
        self.stdout = stdout


def _fixed_run(stdout="", returncode=0):
    def run(*_args, **_kwargs):
        return FakeCompleted(returncode=returncode, stdout=stdout)
    return run


# ---------------------------------------------------------------------------
# source_provenance -- git
# ---------------------------------------------------------------------------

def test_capture_git_provenance_clean_tree():
    # git is invoked twice (rev-parse, status); a stateful stub gives each
    # call its matching canned output.
    calls = iter([
        FakeCompleted(0, "d" * 40 + "\n"),
        FakeCompleted(0, ""),
    ])

    def run(_argv, **_kwargs):
        return next(calls)

    result = source_provenance.capture_git_provenance(run=run)
    assert result["git_commit"] == "d" * 40
    assert result["git_dirty"] is False
    assert result["git_dirty_files"] == []
    assert result["metadata_incomplete"] is False


def test_capture_git_provenance_dirty_tree_lists_files():
    calls = iter([
        FakeCompleted(0, "c" * 40 + "\n"),
        FakeCompleted(0, " M benchmarks/harness/runner.py\n?? new_file.py\n"),
    ])

    def run(_argv, **_kwargs):
        return next(calls)

    result = source_provenance.capture_git_provenance(run=run)
    assert result["git_dirty"] is True
    assert result["git_dirty_files"] == ["benchmarks/harness/runner.py", "new_file.py"]
    assert result["metadata_incomplete"] is False


def test_capture_git_provenance_missing_git_is_honestly_incomplete():
    def run(_argv, **_kwargs):
        raise FileNotFoundError("git")

    result = source_provenance.capture_git_provenance(run=run)
    assert result["git_commit"] == "0" * 40
    assert result["git_dirty"] is False
    assert result["metadata_incomplete"] is True


# ---------------------------------------------------------------------------
# source_provenance -- CMakeCache
# ---------------------------------------------------------------------------

_SAMPLE_CACHE = """\
# This is the CMakeCache file.
CMAKE_BUILD_TYPE:STRING=Release
CMAKE_Fortran_COMPILER:STRING=/bin/f95
CMAKE_Fortran_FLAGS:STRING= -ffree-line-length-0 -cpp
CMAKE_Fortran_FLAGS_RELEASE:STRING=-O3 -DNDEBUG
UPPASD_GPU_BACKEND:STRING=CUDA
UPPASD_PRECISION:STRING=SINGLE
USE_CUDA:BOOL=OFF
//ADVANCED property for variable: CMAKE_Fortran_COMPILER
CMAKE_Fortran_COMPILER-ADVANCED:INTERNAL=1
"""


def test_parse_cmake_cache_and_option_selection(tmp_path):
    cache_path = tmp_path / "CMakeCache.txt"
    cache_path.write_text(_SAMPLE_CACHE)
    cache = source_provenance.parse_cmake_cache(cache_path)
    assert cache["UPPASD_PRECISION"] == "SINGLE"
    assert cache["USE_CUDA"] is False  # BOOL decoded
    assert cache["CMAKE_Fortran_COMPILER-ADVANCED"] == "1"

    options = source_provenance.select_cmake_options(cache)
    assert options["UPPASD_GPU_BACKEND"] == "CUDA"
    assert "CMAKE_Fortran_COMPILER-ADVANCED" not in options


def test_resolve_requested_precision_rejects_unknown():
    with pytest.raises(source_provenance.SourceProvenanceError):
        source_provenance.resolve_requested_precision({"UPPASD_PRECISION": "BOGUS"})
    with pytest.raises(source_provenance.SourceProvenanceError):
        source_provenance.resolve_requested_precision({})


def test_resolve_backend_maps_off_and_gpu_backends():
    assert source_provenance.resolve_backend({"UPPASD_GPU_BACKEND": "OFF"}) == ("CPU", None)
    assert source_provenance.resolve_backend({}) == ("CPU", None)
    assert source_provenance.resolve_backend({"UPPASD_GPU_BACKEND": "CUDA"}) == ("GPU", "CUDA")
    with pytest.raises(source_provenance.SourceProvenanceError):
        source_provenance.resolve_backend({"UPPASD_GPU_BACKEND": "OPENCL"})


def test_build_source_context_missing_cache_is_honestly_incomplete(tmp_path):
    binary_path = tmp_path / "sd.fake"
    binary_path.write_bytes(b"not a real binary")
    build_dir = tmp_path / "build_missing"
    build_dir.mkdir()

    def run(_argv, **_kwargs):
        return FakeCompleted(0, "0" * 40 + "\n")

    context = source_provenance.build_source_context(binary_path, build_dir, run=run)
    assert context["metadata_incomplete"] is True
    assert context["compiler"] == "unaudited"
    assert context["requested_precision"] == "DOUBLE"
    assert context["gpu_backend"] is None


def test_build_source_context_full_pipeline(tmp_path):
    build_dir = tmp_path / "build_gpu"
    build_dir.mkdir()
    (build_dir / "CMakeCache.txt").write_text(_SAMPLE_CACHE)
    binary_path = tmp_path / "sd.cuda"
    binary_path.write_bytes(b"binary contents")

    calls = iter([
        FakeCompleted(0, "a" * 40 + "\n"),  # git rev-parse
        FakeCompleted(0, ""),                # git status (clean)
        FakeCompleted(0, "GNU Fortran 13.3.0\n"),  # compiler --version
    ])

    def run(_argv, **_kwargs):
        return next(calls)

    context = source_provenance.build_source_context(binary_path, build_dir, run=run)
    assert context["requested_precision"] == "SINGLE"
    assert context["backend"] == "GPU"
    assert context["gpu_backend"] == "CUDA"
    assert context["precision_support_state"] == "unaudited"
    assert context["compiler"] == "f95"
    assert context["compiler_version"] == "GNU Fortran 13.3.0"
    assert "-O3" in context["compile_flags"]
    assert context["cmake_options"]["UPPASD_GPU_BACKEND"] == "CUDA"
    assert context["build_id"].startswith("build_gpu-")
    assert context["metadata_incomplete"] is False


def test_build_source_context_mixed_precision_is_unsupported(tmp_path):
    build_dir = tmp_path / "build_mixed"
    build_dir.mkdir()
    (build_dir / "CMakeCache.txt").write_text(_SAMPLE_CACHE.replace("SINGLE", "MIXED"))
    binary_path = tmp_path / "sd.mixed"
    binary_path.write_bytes(b"x")

    def run(_argv, **_kwargs):
        return FakeCompleted(0, "b" * 40 + "\n")

    context = source_provenance.build_source_context(binary_path, build_dir, run=run)
    assert context["requested_precision"] == "MIXED"
    assert context["precision_support_state"] == "unsupported"


# ---------------------------------------------------------------------------
# cpu_provenance
# ---------------------------------------------------------------------------

_SAMPLE_CPUINFO = "\n\n".join(
    f"processor\t: {i}\nmodel name\t: Fake CPU\nphysical id\t: {i // 2}\ncore id\t: {i % 2}"
    for i in range(4)
)


def test_parse_cpuinfo_two_sockets_two_cores_each():
    model, sockets, physical_cores, threads = cpu_provenance.parse_cpuinfo(_SAMPLE_CPUINFO)
    assert model == "Fake CPU"
    assert sockets == 2
    assert physical_cores == 4  # (socket, core) pairs are all distinct here
    assert threads == 4


def test_capture_cpu_topology_from_synthetic_tree(tmp_path):
    cpuinfo_path = tmp_path / "cpuinfo"
    cpuinfo_path.write_text(_SAMPLE_CPUINFO)
    node_root = tmp_path / "node"
    node_root.mkdir()
    (node_root / "node0").mkdir()
    (node_root / "node1").mkdir()
    (node_root / "not_a_node").mkdir()

    topology = cpu_provenance.capture_cpu_topology(cpuinfo_path=cpuinfo_path, node_root=node_root)
    assert topology["cpu_model"] == "Fake CPU"
    assert topology["numa_nodes"] == 2
    assert topology["metadata_incomplete"] is False


def test_capture_cpu_topology_missing_paths_is_honestly_incomplete(tmp_path):
    topology = cpu_provenance.capture_cpu_topology(
        cpuinfo_path=tmp_path / "no-such-file", node_root=tmp_path / "no-such-dir"
    )
    assert topology["cpu_model"] == "unaudited"
    assert topology["metadata_incomplete"] is True


def test_capture_omp_environment_parses_and_defaults():
    env = {"OMP_NUM_THREADS": "8", "OMP_PLACES": "cores", "OMP_PROC_BIND": "close", "OMP_DYNAMIC": "TRUE"}
    result = cpu_provenance.capture_omp_environment(env)
    assert result == {"omp_threads": 8, "omp_places": "cores", "omp_proc_bind": "close", "omp_dynamic": True}

    empty = cpu_provenance.capture_omp_environment({}, default_threads=1)
    assert empty["omp_threads"] == 1
    assert empty["omp_places"] is None
    assert empty["omp_dynamic"] is None

    garbage = cpu_provenance.capture_omp_environment({"OMP_DYNAMIC": "maybe"})
    assert garbage["omp_dynamic"] is None


def test_format_affinity_ranges():
    assert cpu_provenance._format_affinity_ranges([0, 1, 2, 3, 5, 7, 8]) == "0-3,5,7-8"
    assert cpu_provenance._format_affinity_ranges([]) == ""
    assert cpu_provenance._format_affinity_ranges([4]) == "4"


def test_capture_cpu_frequency_state_from_synthetic_sysfs(tmp_path):
    for cpu_id, governor, freq in ((0, "performance", 3000000), (1, "powersave", 1200000)):
        cpufreq_dir = tmp_path / f"cpu{cpu_id}" / "cpufreq"
        cpufreq_dir.mkdir(parents=True)
        (cpufreq_dir / "scaling_governor").write_text(governor + "\n")
        (cpufreq_dir / "scaling_cur_freq").write_text(str(freq) + "\n")

    state = cpu_provenance.capture_cpu_frequency_state([0, 1, 2], sysfs_root=tmp_path)
    assert state["governors"] == {0: "performance", 1: "powersave"}
    assert state["current_freq_khz"] == {0: 3000000, 1: 1200000}
    assert 2 not in state["governors"]  # cpu2 has no cpufreq dir -- honestly omitted


# ---------------------------------------------------------------------------
# gpu_provenance -- CUDA (mock nvidia-smi output)
# ---------------------------------------------------------------------------

def _cuda_csv_row(**overrides):
    fields = {
        "index": "0", "uuid": "GPU-abc", "name": "NVIDIA RTX A4000", "driver_version": "610.57.04",
        "memory.total": "16376", "memory.used": "102", "utilization.gpu": "0", "temperature.gpu": "45",
        "clocks.sm": "210", "clocks.max.sm": "2100",
        "clocks_throttle_reasons.sw_thermal_slowdown": "Not Active",
        "clocks_throttle_reasons.hw_thermal_slowdown": "Not Active",
        "clocks_throttle_reasons.hw_slowdown": "Not Active",
        "clocks_throttle_reasons.hw_power_brake_slowdown": "Not Active",
        "clocks_throttle_reasons.sw_power_cap": "Not Active",
    }
    fields.update(overrides)
    return ", ".join(fields[name] for name in gpu_provenance._CUDA_QUERY_FIELDS)


def test_query_cuda_devices_parses_canned_csv(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: True)
    run = _fixed_run(_cuda_csv_row() + "\n")

    devices = gpu_provenance.query_cuda_devices(run=run)
    assert len(devices) == 1
    device = devices[0]
    assert device["model"] == "NVIDIA RTX A4000"
    assert device["uuid"] == "GPU-abc"
    assert device["driver_version"] == "610.57.04"
    assert device["temperature_c"] == 45
    assert device["clocks_sm_mhz"] == 210
    assert device["thermal_throttled"] is False
    assert device["throttled"] is False


def test_query_cuda_devices_detects_thermal_throttle(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: True)
    row = _cuda_csv_row(**{"clocks_throttle_reasons.hw_thermal_slowdown": "Active"})
    devices = gpu_provenance.query_cuda_devices(run=_fixed_run(row + "\n"))
    assert devices[0]["thermal_throttled"] is True
    assert devices[0]["throttled"] is True


def test_query_cuda_devices_missing_nvidia_smi_is_honestly_empty(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: False)
    assert gpu_provenance.query_cuda_devices(run=_fixed_run("should not be called")) == []


def test_query_cuda_devices_command_failure_is_honestly_empty(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: True)

    def run(_argv, **_kwargs):
        raise OSError("nvidia-smi not runnable")

    assert gpu_provenance.query_cuda_devices(run=run) == []


def test_query_cuda_runtime_version_parses_nvcc(monkeypatch):
    monkeypatch.setattr(gpu_provenance.shutil, "which", lambda name: "/usr/bin/nvcc" if name == "nvcc" else None)
    run = _fixed_run("Cuda compilation tools, release 13.3, V13.3.73\n")
    assert gpu_provenance.query_cuda_runtime_version(run=run) == "13.3"


def test_query_cuda_runtime_version_missing_nvcc_is_none(monkeypatch):
    monkeypatch.setattr(gpu_provenance.shutil, "which", lambda _name: None)
    assert gpu_provenance.query_cuda_runtime_version(run=_fixed_run("unused")) is None


def test_classify_cuda_contamination_excludes_own_pid(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: True)
    device = {"thermal_throttled": False}
    own_process = [{"pid": 42, "process_name": "sd.f95.cuda"}]
    assert gpu_provenance.classify_cuda_contamination(device, own_process, own_pid=42) == set()

    competing = [{"pid": 42, "process_name": "sd.f95.cuda"}, {"pid": 999, "process_name": "someone-else"}]
    assert gpu_provenance.classify_cuda_contamination(device, competing, own_pid=42) == {
        gpu_provenance.QUALITY_FLAG_GPU_BUSY
    }


def test_classify_cuda_contamination_thermal_throttle_flag():
    device = {"thermal_throttled": True}
    flags = gpu_provenance.classify_cuda_contamination(device, [], own_pid=1)
    assert flags == {gpu_provenance.QUALITY_FLAG_GPU_THERMAL_THROTTLE}


def test_detect_cuda_clock_instability():
    before = {"clocks_sm_mhz": 1000}
    stable_after = {"clocks_sm_mhz": 1050}
    unstable_after = {"clocks_sm_mhz": 300}
    assert gpu_provenance.detect_cuda_clock_instability(before, stable_after) == set()
    assert gpu_provenance.detect_cuda_clock_instability(before, unstable_after) == {
        gpu_provenance.QUALITY_FLAG_GPU_CLOCK_UNSTABLE
    }
    assert gpu_provenance.detect_cuda_clock_instability(None, unstable_after) == set()


def test_query_cuda_compute_processes_parses_csv(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "nvidia_smi_available", lambda: True)
    run = _fixed_run("12345, other-tool, 512\n")
    processes = gpu_provenance.query_cuda_compute_processes(0, run=run)
    assert processes == [{"pid": 12345, "process_name": "other-tool", "used_memory_mib": 512}]


# ---------------------------------------------------------------------------
# gpu_provenance -- HIP (structural, no AMD hardware required)
# ---------------------------------------------------------------------------

_SAMPLE_ROCM_JSON = json.dumps({
    "card0": {
        "Card series": "Fake MI-Series",
        "Unique ID": "0x1234",
        "Driver version": "6.0.0",
        "Temperature (Sensor edge) (C)": "50.0",
        "GPU use (%)": "0",
        "sclk clock speed": "800Mhz",
        "Throttle Status": "0",
    }
})


def test_query_hip_devices_missing_rocm_smi_is_honestly_empty(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "rocm_smi_available", lambda: False)
    assert gpu_provenance.query_hip_devices(run=_fixed_run("should not be called")) == []


def test_query_hip_devices_parses_canned_json(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "rocm_smi_available", lambda: True)
    devices = gpu_provenance.query_hip_devices(run=_fixed_run(_SAMPLE_ROCM_JSON))
    assert len(devices) == 1
    device = devices[0]
    assert device["model"] == "Fake MI-Series"
    assert device["uuid"] == "0x1234"
    assert device["driver_version"] == "6.0.0"
    assert device["temperature_c"] == 50
    assert device["thermal_throttled"] is False


def test_query_hip_devices_detects_thermal_throttle_string(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "rocm_smi_available", lambda: True)
    payload = json.loads(_SAMPLE_ROCM_JSON)
    payload["card0"]["Throttle Status"] = "thermal"
    devices = gpu_provenance.query_hip_devices(run=_fixed_run(json.dumps(payload)))
    assert devices[0]["thermal_throttled"] is True


def test_query_hip_devices_malformed_json_is_honestly_empty(monkeypatch):
    monkeypatch.setattr(gpu_provenance, "rocm_smi_available", lambda: True)
    assert gpu_provenance.query_hip_devices(run=_fixed_run("not json{{{")) == []


def test_gpu_backend_dispatch_rejects_unknown_backend():
    with pytest.raises(ValueError):
        gpu_provenance.query_devices("OPENCL")
    with pytest.raises(ValueError):
        gpu_provenance.classify_contamination("OPENCL", None, [])


# ---------------------------------------------------------------------------
# provenance.py orchestration
# ---------------------------------------------------------------------------

def test_classify_background_load():
    assert provenance.classify_background_load(20.0, 8) == {"background_load_high"}
    assert provenance.classify_background_load(1.0, 8) == set()
    assert provenance.classify_background_load(None, 8) == set()
    assert provenance.classify_background_load(5.0, None) == set()


def test_classify_cpu_frequency_stability():
    before = {"governors": {0: "performance", 1: "performance"}}
    same = {"governors": {0: "performance", 1: "performance"}}
    changed = {"governors": {0: "performance", 1: "powersave"}}
    assert provenance.classify_cpu_frequency_stability(before, same) == set()
    assert provenance.classify_cpu_frequency_stability(before, changed) == {"cpu_frequency_unstable"}


def test_gate_sample_strict_mode():
    # Default (non-strict): never raises, regardless of validity.
    provenance.gate_sample(environment_valid=False, quality_flags=["dirty_tree"], require_clean_environment=False)

    # Strict + valid: no-op.
    provenance.gate_sample(environment_valid=True, quality_flags=[], require_clean_environment=True)

    # Strict + invalid: refused.
    with pytest.raises(provenance.EnvironmentGateError):
        provenance.gate_sample(environment_valid=False, quality_flags=["gpu_busy"], require_clean_environment=True)


def test_classify_gpu_environment_unions_busy_and_clock_instability():
    before = provenance.GpuSnapshot(
        device={"thermal_throttled": False, "clocks_sm_mhz": 200},
        processes=[{"pid": 555, "process_name": "intruder"}],
        taken_at=0.0,
    )
    after = provenance.GpuSnapshot(
        device={"thermal_throttled": False, "clocks_sm_mhz": 1800},
        processes=[],
        taken_at=1.0,
    )
    flags = provenance.classify_gpu_environment("CUDA", before, after, own_pid=1)
    assert gpu_provenance.QUALITY_FLAG_GPU_BUSY in flags
    assert gpu_provenance.QUALITY_FLAG_GPU_CLOCK_UNSTABLE in flags


def test_build_static_context_dirty_tree_raises_flag(tmp_path):
    build_dir = tmp_path / "build_cpu"
    build_dir.mkdir()
    (build_dir / "CMakeCache.txt").write_text(_SAMPLE_CACHE.replace("CUDA", "OFF"))
    binary_path = tmp_path / "sd.f95"
    binary_path.write_bytes(b"binary")

    calls = iter([
        FakeCompleted(0, "e" * 40 + "\n"),
        FakeCompleted(0, " M source/somefile.f90\n"),  # dirty
        FakeCompleted(0, "GNU Fortran 13.3.0\n"),
    ])

    def run(_argv, **_kwargs):
        return next(calls)

    context = provenance.build_static_context(
        binary_path=binary_path, build_dir=build_dir, machine_id="unit-test", env={}, run=run,
    )
    assert "dirty_tree" in context["quality_flags"]
    assert context["backend"] == "CPU"
    assert context["gpu_backend"] is None


def test_build_static_context_missing_gpu_device_flags_incomplete(tmp_path):
    build_dir = tmp_path / "build_gpu"
    build_dir.mkdir()
    (build_dir / "CMakeCache.txt").write_text(_SAMPLE_CACHE)
    binary_path = tmp_path / "sd.cuda"
    binary_path.write_bytes(b"binary")

    calls = iter([
        FakeCompleted(0, "f" * 40 + "\n"),
        FakeCompleted(0, ""),
        FakeCompleted(0, "GNU Fortran 13.3.0\n"),
    ])

    def run(_argv, **_kwargs):
        try:
            return next(calls)
        except StopIteration:
            # nvidia-smi / nvcc calls: simulate "no device reachable".
            return FakeCompleted(1, "")

    context = provenance.build_static_context(
        binary_path=binary_path, build_dir=build_dir, machine_id="unit-test", env={}, run=run,
    )
    assert context["gpu_model"] is None
    assert "metadata_incomplete" in context["quality_flags"]


# ---------------------------------------------------------------------------
# runner.py strict-mode wiring
# ---------------------------------------------------------------------------

def test_runner_environment_gate_error_is_a_runner_error():
    from harness import runner
    assert issubclass(runner.EnvironmentGateError, runner.RunnerError)
