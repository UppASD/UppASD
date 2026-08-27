"""GPU identity and contention/throttle provenance (WP-04 section C).

CUDA capture goes through `nvidia-smi` (device identity/state) and `nvcc
--version` (toolkit/runtime version); HIP exposes the structurally analogous
functions against `rocm-smi`, per the blueprint: representable without
requiring AMD hardware, but no HIP performance claim may ever be backed by
data this module produces without real execution on real hardware.

Every external call goes through an injectable ``run`` callable (default
:func:`subprocess.run`), and availability is checked before use -- there is
no AMD hardware in this environment, so the HIP path is exercised exclusively
through the mock-parser tests required by section G, never against a real
`rocm-smi`.

Nothing here decides `environment_valid` on its own: these functions capture
and classify; `provenance.py` and the campaign driver decide what to do with
the resulting flags.
"""

from __future__ import annotations

import csv
import io
import json
import re
import shutil
import subprocess

QUALITY_FLAG_GPU_BUSY = "gpu_busy"
QUALITY_FLAG_GPU_THERMAL_THROTTLE = "gpu_thermal_throttle"
QUALITY_FLAG_GPU_CLOCK_UNSTABLE = "gpu_clock_unstable"

# ---------------------------------------------------------------------------
# CUDA
# ---------------------------------------------------------------------------

_CUDA_QUERY_FIELDS = (
    "index",
    "uuid",
    "name",
    "driver_version",
    "memory.total",
    "memory.used",
    "utilization.gpu",
    "temperature.gpu",
    "clocks.sm",
    "clocks.max.sm",
    "clocks_throttle_reasons.sw_thermal_slowdown",
    "clocks_throttle_reasons.hw_thermal_slowdown",
    "clocks_throttle_reasons.hw_slowdown",
    "clocks_throttle_reasons.hw_power_brake_slowdown",
    "clocks_throttle_reasons.sw_power_cap",
)

_ACTIVE_TOKENS = {"active", "true", "1"}


def nvidia_smi_available(*, which=shutil.which):
    return which("nvidia-smi") is not None


def _parse_csv_rows(stdout, field_names):
    reader = csv.reader(io.StringIO(stdout), skipinitialspace=True)
    rows = []
    for raw in reader:
        if not raw or len(raw) != len(field_names):
            continue
        rows.append(dict(zip(field_names, (cell.strip() for cell in raw))))
    return rows


def _int_or_none(token):
    try:
        return int(token)
    except (TypeError, ValueError):
        pass
    try:
        return int(round(float(token)))
    except (TypeError, ValueError):
        return None


def query_cuda_devices(*, run=subprocess.run, timeout=10):
    """List every CUDA device `nvidia-smi` reports, parsed and typed.

    Returns `[]` (not an error) when `nvidia-smi` is missing or fails --
    "missing metadata handled honestly" (section G) means a GPU-backend
    record with no reachable device still gets a record, with
    `metadata_incomplete` raised by the caller rather than a fabricated
    device.
    """
    if not nvidia_smi_available():
        return []
    try:
        result = run(
            ["nvidia-smi", f"--query-gpu={','.join(_CUDA_QUERY_FIELDS)}", "--format=csv,noheader,nounits"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    if result.returncode != 0:
        return []
    return [_parse_cuda_device_row(row) for row in _parse_csv_rows(result.stdout, _CUDA_QUERY_FIELDS)]


def _parse_cuda_device_row(row):
    throttle_active = any(
        row.get(field, "").strip().lower() in _ACTIVE_TOKENS
        for field in _CUDA_QUERY_FIELDS
        if field.startswith("clocks_throttle_reasons.")
    )
    thermal_active = any(
        row.get(field, "").strip().lower() in _ACTIVE_TOKENS
        for field in ("clocks_throttle_reasons.sw_thermal_slowdown", "clocks_throttle_reasons.hw_thermal_slowdown")
    )
    return {
        "index": row.get("index"),
        "uuid": row.get("uuid") or None,
        "model": row.get("name") or None,
        "driver_version": row.get("driver_version") or None,
        "memory_total_mib": _int_or_none(row.get("memory.total")),
        "memory_used_mib": _int_or_none(row.get("memory.used")),
        "utilization_gpu_percent": _int_or_none(row.get("utilization.gpu")),
        "temperature_c": _int_or_none(row.get("temperature.gpu")),
        "clocks_sm_mhz": _int_or_none(row.get("clocks.sm")),
        "clocks_max_sm_mhz": _int_or_none(row.get("clocks.max.sm")),
        "throttled": throttle_active,
        "thermal_throttled": thermal_active,
    }


def query_cuda_runtime_version(*, run=subprocess.run, timeout=10):
    """CUDA toolkit version actually used to build, from `nvcc --version`.

    `None` when `nvcc` is unavailable -- the record then carries no
    fabricated runtime version rather than one guessed from the driver.
    """
    nvcc = shutil.which("nvcc")
    if nvcc is None:
        return None
    try:
        result = run([nvcc, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"release\s+([\d.]+)", result.stdout or "")
    return match.group(1) if match else None


def query_cuda_compute_processes(gpu_index, *, run=subprocess.run, timeout=10):
    """Processes currently using the given CUDA device -- "visible competing processes"."""
    if not nvidia_smi_available():
        return []
    try:
        result = run(
            [
                "nvidia-smi", f"-i={gpu_index}",
                "--query-compute-apps=pid,process_name,used_memory",
                "--format=csv,noheader,nounits",
            ],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    if result.returncode != 0:
        return []
    fields = ("pid", "process_name", "used_memory")
    processes = []
    for row in _parse_csv_rows(result.stdout, fields):
        processes.append({
            "pid": _int_or_none(row.get("pid")),
            "process_name": row.get("process_name") or None,
            "used_memory_mib": _int_or_none(row.get("used_memory")),
        })
    return processes


def classify_cuda_contamination(device, processes, *, own_pid=None):
    """Quality flags from one CUDA device snapshot + its compute-process list.

    `gpu_busy` fires only for a process that is not this benchmark's own
    (`own_pid` excluded) -- a solitary sample process legitimately using the
    GPU is not contamination.
    """
    flags = set()
    if device and device.get("thermal_throttled"):
        flags.add(QUALITY_FLAG_GPU_THERMAL_THROTTLE)
    competing = [p for p in processes if p.get("pid") is not None and p.get("pid") != own_pid]
    if competing:
        flags.add(QUALITY_FLAG_GPU_BUSY)
    return flags


def detect_cuda_clock_instability(before, after, *, relative_band=0.10):
    """`gpu_clock_unstable` if the SM clock moved more than `relative_band` across a run."""
    if not before or not after:
        return set()
    b, a = before.get("clocks_sm_mhz"), after.get("clocks_sm_mhz")
    if b is None or a is None or b == 0:
        return set()
    if abs(a - b) / b > relative_band:
        return {QUALITY_FLAG_GPU_CLOCK_UNSTABLE}
    return set()


# ---------------------------------------------------------------------------
# HIP -- structurally analogous, unverified against real hardware
# ---------------------------------------------------------------------------

def rocm_smi_available(*, which=shutil.which):
    return which("rocm-smi") is not None


def query_hip_devices(*, run=subprocess.run, timeout=10):
    """List every HIP device `rocm-smi --showallinfo --json` reports.

    `rocm-smi`'s JSON key names have varied across ROCm releases; this reads
    the handful of keys documented as stable (`Device Name`/`Card series`,
    `Unique ID`/`GPU ID`, `Driver version`, temperature and GPU-use
    percentages) defensively via `.get()`, and never raises on an unexpected
    shape -- it returns as much as it can parse. There is no AMD hardware to
    validate this against in this environment; it is exercised only by the
    mock-parser tests (section G), and no HIP performance claim may rely on
    it without real execution on real hardware.
    """
    if not rocm_smi_available():
        return []
    try:
        result = run(
            ["rocm-smi", "--showallinfo", "--json"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    if result.returncode != 0:
        return []
    try:
        payload = json.loads(result.stdout or "{}")
    except json.JSONDecodeError:
        return []
    return [_parse_hip_device_entry(index, card, entry) for index, (card, entry) in enumerate(payload.items())]


def _parse_hip_device_entry(index, card, entry):
    def first(*keys):
        for key in keys:
            if key in entry:
                return entry[key]
        return None

    temperature = first("Temperature (Sensor edge) (C)", "Temperature (Sensor junction) (C)")
    util = first("GPU use (%)")
    sm_clock = first("sclk clock speed", "sclk clock speed:")
    throttle = first("Throttle Status", "throttle_status")
    return {
        "index": index,
        "uuid": first("Unique ID", "GPU ID") or f"card{card}",
        "model": first("Card series", "Device Name") or None,
        "driver_version": first("Driver version") or None,
        "memory_total_mib": _int_or_none(first("VRAM Total Memory (B)")),
        "memory_used_mib": _int_or_none(first("VRAM Total Used Memory (B)")),
        "utilization_gpu_percent": _int_or_none(util),
        "temperature_c": _int_or_none(temperature),
        "clocks_sm_mhz": _int_or_none(re.sub(r"[^\d]", "", sm_clock) if isinstance(sm_clock, str) else sm_clock),
        "clocks_max_sm_mhz": None,
        "throttled": bool(throttle) and str(throttle).strip().lower() not in ("", "0", "none", "normal"),
        "thermal_throttled": bool(throttle) and "therm" in str(throttle).lower(),
    }


def query_hip_runtime_version(*, run=subprocess.run, timeout=10):
    """ROCm/HIP toolkit version from `hipcc --version`, or `None`."""
    hipcc = shutil.which("hipcc")
    if hipcc is None:
        return None
    try:
        result = run([hipcc, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"HIP version:\s*([\d.]+)", result.stdout or "")
    return match.group(1) if match else None


def classify_hip_contamination(device, *, own_pid=None):
    """Structurally analogous to `classify_cuda_contamination`.

    `rocm-smi` has no equivalent of `nvidia-smi --query-compute-apps` in
    widespread stable form, so competing-process detection is not attempted
    here -- only device-reported throttle state -- until real HIP hardware
    makes a process-level check verifiable.
    """
    flags = set()
    if device and device.get("thermal_throttled"):
        flags.add(QUALITY_FLAG_GPU_THERMAL_THROTTLE)
    return flags


def detect_hip_clock_instability(before, after, *, relative_band=0.10):
    return detect_cuda_clock_instability(before, after, relative_band=relative_band)


# ---------------------------------------------------------------------------
# Backend-dispatching helpers
# ---------------------------------------------------------------------------

def query_devices(gpu_backend, *, run=subprocess.run):
    if gpu_backend == "CUDA":
        return query_cuda_devices(run=run)
    if gpu_backend == "HIP":
        return query_hip_devices(run=run)
    raise ValueError(f"unknown gpu_backend {gpu_backend!r}")


def query_runtime_version(gpu_backend, *, run=subprocess.run):
    if gpu_backend == "CUDA":
        return query_cuda_runtime_version(run=run)
    if gpu_backend == "HIP":
        return query_hip_runtime_version(run=run)
    raise ValueError(f"unknown gpu_backend {gpu_backend!r}")


def classify_contamination(gpu_backend, device, processes, *, own_pid=None):
    if gpu_backend == "CUDA":
        return classify_cuda_contamination(device, processes, own_pid=own_pid)
    if gpu_backend == "HIP":
        return classify_hip_contamination(device, own_pid=own_pid)
    raise ValueError(f"unknown gpu_backend {gpu_backend!r}")


def detect_clock_instability(gpu_backend, before, after, *, relative_band=0.10):
    if gpu_backend == "CUDA":
        return detect_cuda_clock_instability(before, after, relative_band=relative_band)
    if gpu_backend == "HIP":
        return detect_hip_clock_instability(before, after, relative_band=relative_band)
    raise ValueError(f"unknown gpu_backend {gpu_backend!r}")
