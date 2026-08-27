"""Provenance orchestration and environment-quality classification (WP-04).

Composes `source_provenance.py` (section A), `cpu_provenance.py` (section B)
and `gpu_provenance.py` (section C) into one `context` dict shaped exactly
like what `runner.run_sample` requires (see `runner._REQUIRED_CONTEXT_FIELDS`)
-- this is the authoritative counterpart to `runner.developer_context`, meant
for real campaign use rather than local self-tests.

Static identity (source/build/CPU topology/GPU identity) is resolved once per
(build, machine) via `build_static_context` and reused across samples, the
same way `runner.resolve_workload_metadata` is resolved once per cell. GPU
contention/throttle/clock state is inherently per-sample -- a campaign driver
takes a `GpuSnapshot` immediately before and after each `run_sample` call and
passes both to `classify_gpu_environment`, then folds the resulting flags
into that sample's `context["quality_flags"]` before writing it (or, more
simply, into `notes`/a post-hoc quality pass, since raw sample records are
never rewritten in place).

Section D's quality flags and section E's strict mode live here too:
`classify_environment_flags` assembles the closed flag set from real signals,
and `EnvironmentPolicy`/`gate_sample` decide whether a flagged sample may be
accepted as authoritative. Developer runs are never blocked by this module on
their own -- only a caller that opts into strict mode is.
"""

from __future__ import annotations

import subprocess
import time

from harness import cpu_provenance
from harness import gpu_provenance
from harness import source_provenance

BACKGROUND_LOAD_RATIO_THRESHOLD = 1.15


def build_static_context(
    *,
    binary_path,
    build_dir,
    machine_id,
    repo_root=None,
    gpu_index=0,
    env=None,
    run=subprocess.run,
    language="Fortran",
):
    """Resolve everything about a (binary, build, machine) that does not
    change between samples: source/build identity, CPU topology, the OpenMP
    environment and static GPU identity (model/driver/runtime -- not its
    live contention/throttle state).

    Returns a dict directly usable as `runner.run_sample`'s `context`
    (`quality_flags` included), plus the source `backend`/`gpu_backend` this
    build was actually configured for (read from `CMakeCache.txt`, never
    guessed from the binary's path).
    """
    source = source_provenance.build_source_context(binary_path, build_dir, repo_root=repo_root, run=run)
    topology = cpu_provenance.capture_cpu_topology()
    omp = cpu_provenance.capture_omp_environment(env)
    affinity = cpu_provenance.capture_process_affinity()

    flags = set()
    if source["metadata_incomplete"] or topology["metadata_incomplete"]:
        flags.add("metadata_incomplete")
    if affinity is None:
        flags.add("cpu_affinity_unknown")
    if source["git_dirty"]:
        flags.add("dirty_tree")

    gpu_model = gpu_id = gpu_driver = gpu_runtime = None
    backend = source["backend"]
    gpu_backend = source["gpu_backend"]
    if gpu_backend is not None:
        devices = gpu_provenance.query_devices(gpu_backend, run=run)
        device = next((d for d in devices if str(d.get("index")) == str(gpu_index)), None)
        if device is None:
            flags.add("metadata_incomplete")
        else:
            gpu_model = device.get("model")
            gpu_id = device.get("uuid") or (str(device["index"]) if device.get("index") is not None else None)
            gpu_driver = device.get("driver_version")
        gpu_runtime = gpu_provenance.query_runtime_version(gpu_backend, run=run)
        if gpu_runtime is None:
            flags.add("metadata_incomplete")

    context = {
        "machine_id": machine_id,
        "backend": backend,
        "gpu_backend": gpu_backend,
        "omp_threads": omp["omp_threads"],
        "omp_places": omp["omp_places"],
        "omp_proc_bind": omp["omp_proc_bind"],
        "omp_dynamic": omp["omp_dynamic"],
        "process_affinity": affinity,
        "requested_precision": source["requested_precision"],
        "effective_cpu_precision": "UNKNOWN",
        "effective_gpu_precision": "UNKNOWN" if gpu_backend else None,
        "precision_support_state": source["precision_support_state"],
        "comparison_precision_class": None if source["precision_support_state"] == "unsupported" else "UNAUDITED",
        "git_commit": source["git_commit"],
        "git_dirty": source["git_dirty"],
        "git_dirty_files": source["git_dirty_files"],
        "compiler": source["compiler"],
        "compiler_version": source["compiler_version"],
        "compile_flags": source["compile_flags"],
        "cmake_options": source["cmake_options"],
        "binary_checksum": source["binary_checksum"],
        "build_id": source["build_id"],
        "cpu_model": topology["cpu_model"],
        "cpu_sockets": topology["cpu_sockets"],
        "cpu_physical_cores": topology["cpu_physical_cores"],
        "cpu_threads": topology["cpu_threads"],
        "numa_nodes": topology["numa_nodes"],
        "gpu_model": gpu_model,
        "gpu_id": gpu_id,
        "gpu_driver": gpu_driver,
        "gpu_runtime": gpu_runtime,
        "dipole_backend": None,
        "quality_flags": sorted(flags),
    }
    return context


class GpuSnapshot:
    """One point-in-time GPU device + compute-process observation."""

    def __init__(self, device, processes, taken_at):
        self.device = device
        self.processes = processes
        self.taken_at = taken_at


def capture_gpu_snapshot(gpu_backend, gpu_index=0, *, run=subprocess.run):
    """A `GpuSnapshot` for `gpu_backend`'s device at `gpu_index`, right now.

    Call once immediately before and once immediately after a `run_sample`
    invocation to bracket that sample's execution window; `None` device means
    the device could not be reached at that instant (still a valid, honestly
    incomplete snapshot -- not a crash).
    """
    devices = gpu_provenance.query_devices(gpu_backend, run=run)
    device = next((d for d in devices if str(d.get("index")) == str(gpu_index)), None)
    processes = gpu_provenance.query_cuda_compute_processes(gpu_index, run=run) if gpu_backend == "CUDA" else []
    return GpuSnapshot(device=device, processes=processes, taken_at=time.time())


def classify_gpu_environment(gpu_backend, before, after, *, own_pid=None):
    """Union of contamination + clock-instability flags across one bracketed run.

    `before`/`after` are `GpuSnapshot`s from `capture_gpu_snapshot`. Missing
    devices on either side contribute no flags of their own here (that is
    `metadata_incomplete` territory, decided by the caller) but also cannot
    hide contamination that the reachable snapshot did detect.
    """
    flags = set()
    for snapshot in (before, after):
        flags |= gpu_provenance.classify_contamination(gpu_backend, snapshot.device, snapshot.processes, own_pid=own_pid)
    flags |= gpu_provenance.detect_clock_instability(gpu_backend, before.device, after.device)
    return flags


def classify_background_load(load1, cpu_count, *, threshold=BACKGROUND_LOAD_RATIO_THRESHOLD):
    """`background_load_high` when 1-minute load exceeds `threshold * cpu_count`.

    `None` inputs (load average unavailable on this platform) classify as no
    flag -- absence of evidence is not evidence of a busy machine.
    """
    if load1 is None or not cpu_count:
        return set()
    return {"background_load_high"} if load1 > threshold * cpu_count else set()


def classify_cpu_frequency_stability(before_state, after_state):
    """`cpu_frequency_unstable` if any core's governor or frequency moved between snapshots.

    Cores absent from one snapshot (only ever seen before or after) are
    ignored rather than treated as an implicit change -- `cpu_provenance`
    already reports missing per-core data honestly by omission.
    """
    governors_before = before_state.get("governors", {})
    governors_after = after_state.get("governors", {})
    for cpu_id, governor in governors_before.items():
        if cpu_id in governors_after and governors_after[cpu_id] != governor:
            return {"cpu_frequency_unstable"}
    return set()


class EnvironmentGateError(ValueError):
    """A strict-mode environment requirement was not met."""


def gate_sample(*, environment_valid, quality_flags, require_clean_environment):
    """Section E: refuse a sample as authoritative evidence in strict mode.

    Raises `EnvironmentGateError` when `require_clean_environment` is set and
    the sample is not `environment_valid`; otherwise a no-op. Never called
    implicitly -- a normal developer run (the default, `require_clean_environment=False`)
    is unaffected by a dirty tree or an unreachable GPU query.
    """
    if require_clean_environment and not environment_valid:
        raise EnvironmentGateError(
            f"sample rejected by --require-clean-environment: quality_flags={sorted(quality_flags)}"
        )
