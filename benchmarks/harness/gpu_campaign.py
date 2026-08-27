"""CUDA precision campaign support (WP-06 sections C-E).

Three things a GPU precision campaign needs beyond the WP-03/04/05 CPU
machinery:

1. **A fixed host OpenMP configuration** (section D). Source inspection
   (`source/uppasd.f90:823-826` RNG init, `source/Hamiltonian/dipolecommon.f90`
   / `fftdipole_fftw.f90` dipole setup, `source/Fields/spintorques.f90`
   spin-transfer-torque field construction) shows real `!$omp parallel`
   host-side work runs even in a GPU production run -- at minimum once
   during setup, and per-step whenever the case's variant uses STT. Nothing
   in that inspection showed *every* case needing more than one host thread
   for its per-step GPU loop, so `GPU_HOST_OMP_THREADS_DEFAULT = 1` is used
   as the section D default, exactly as the blueprint prefers -- but it is
   an initial choice, not a proof that 1 thread is never a bottleneck (the
   per-step host/GPU interaction was explicitly out of scope to fully
   characterize here; "Do not optimize host OMP+GPU interaction in this
   task"). Record `omp_threads` on every GPU record regardless, exactly as
   CPU records already do.

2. **MIXED emitted as `unsupported`, not a failure** (section C).
   `UPPASD_PRECISION=MIXED` fatal-errors at CMake configure time
   (`CMakeLists.txt:159-163`) -- there is no binary to run, so
   `build_unsupported_precision_record` writes the schema's first-class
   unsupported-state record directly, the same shape as
   `tests/fixtures/valid/cuda_mixed_precision_unsupported.json`, never
   invoking `runner.run_sample`.

3. **GPU memory-limit handling** (section E). `build_unavailable_memory_record`
   is the `gpu_memory.UNAVAILABLE_MEMORY` counterpart: a real, queried GPU
   device that this campaign has decided *not* to attempt at this size,
   recorded rather than silently reduced. Also written directly, never via
   `runner.run_sample` -- attempting the run is exactly what this exists to
   avoid.

Both record builders take `workload_metadata` as already resolved (the same
contract `runner.run_sample` uses -- see its module docstring): natom is a
property of the case/size, resolved once via a cheap CPU probe, independent
of whether the GPU build under audit can actually hold it. `context` is the
same static provenance/hardware shape `runner.run_sample` requires
(`runner._REQUIRED_CONTEXT_FIELDS` -- machine/build/CPU identity, never
`temperature`/`timestep`/`nstep`/`measurement_profile`, which are passed
explicitly here exactly as they are to `runner.run_sample`). For
`build_unsupported_precision_record` specifically, no build of the requested
(unsupported) precision exists to have produced a `context` via
`provenance.build_static_context` -- pass the static context of any real
build already made on this machine (e.g. this campaign's CUDA SINGLE or
DOUBLE build): source/machine identity (`git_commit`, `cpu_model`, ...) is a
property of the machine and repository state, not of the unbuildable
configuration.
"""

from __future__ import annotations

import datetime
import json
import pathlib

from harness import omp_sweep
from harness import precision_audit
from harness import schema_validate
from harness.gpu_memory import UNAVAILABLE_MEMORY

SCHEMA_VERSION = "1.0.0"

GPU_HOST_OMP_THREADS_DEFAULT = 1

SUPPORTED_CUDA_PRECISIONS = ("SINGLE", "DOUBLE")
UNSUPPORTED_PRECISIONS = ("MIXED",)


class GpuCampaignError(ValueError):
    """A GPU campaign record could not be built from the given inputs."""


def resolve_gpu_host_env(*, omp_threads=GPU_HOST_OMP_THREADS_DEFAULT, base_env=None):
    """The fixed, explicit host OpenMP environment for a GPU campaign sample.

    Reuses `omp_sweep.build_omp_env` so a GPU run's host environment is set
    with the same rigor as a CPU sweep's (`OMP_NUM_THREADS`/`OMP_PLACES`/
    `OMP_PROC_BIND`/`OMP_DYNAMIC` all explicit, never inherited ambiently);
    `proc_bind="close"` is an arbitrary but harmless choice at
    `omp_threads == 1`, since there is nothing to bind across.
    """
    return omp_sweep.build_omp_env(omp_threads, proc_bind="close", places="cores", base_env=base_env)


def _base_skipped_record(
    *, run_id, campaign_id, case_id, variant_id, size_id, sample_index, machine_id,
    workflow, workload_class, workload_metadata, backend, gpu_backend, omp_threads,
    requested_precision, effective_cpu_precision, effective_gpu_precision, precision_support_state,
    context, gpu_model, gpu_id, gpu_driver, gpu_runtime, temperature, timestep, nstep, measurement_profile,
    notes, timestamp_utc=None,
):
    comparison_precision_class = precision_audit.classify_comparison_precision_class(
        backend=backend, gpu_backend=gpu_backend, precision_support_state=precision_support_state,
        effective_cpu_precision=effective_cpu_precision, effective_gpu_precision=effective_gpu_precision,
        run_status="SKIPPED",
    )
    stamp = timestamp_utc or datetime.datetime.now(datetime.timezone.utc)

    record = {
        "schema_version": SCHEMA_VERSION,
        "record_kind": "raw_sample",
        "run_id": run_id,
        "campaign_id": campaign_id,
        "case_id": case_id,
        "variant_id": variant_id,
        "size_id": size_id,
        "sample_index": sample_index,
        "machine_id": machine_id,
        "timestamp_utc": stamp.strftime("%Y-%m-%dT%H:%M:%SZ") if hasattr(stamp, "strftime") else stamp,
        "workflow": workflow,
        "workload_class": workload_class,
        "natom": workload_metadata["natom"],
        "directed_interactions": workload_metadata.get("directed_interactions"),
        "mean_neighbors": workload_metadata.get("mean_neighbors"),
        "median_neighbors": workload_metadata.get("median_neighbors"),
        "max_neighbors": workload_metadata.get("max_neighbors"),
        "fft_grid": workload_metadata.get("fft_grid"),
        "fft_grid_padded": workload_metadata.get("fft_grid_padded"),
        "fft_grid_points": workload_metadata.get("fft_grid_points"),
        "dipole_backend": workload_metadata.get("dipole_backend"),
        "backend": backend,
        "gpu_backend": gpu_backend,
        "omp_threads": omp_threads,
        "omp_places": "cores",
        "omp_proc_bind": "close",
        "omp_dynamic": False,
        "process_affinity": None,
        "requested_precision": requested_precision,
        "effective_cpu_precision": effective_cpu_precision,
        "effective_gpu_precision": effective_gpu_precision,
        "precision_support_state": precision_support_state,
        "comparison_precision_class": comparison_precision_class,
        "git_commit": context["git_commit"],
        "git_dirty": context["git_dirty"],
        "git_dirty_files": context.get("git_dirty_files", []),
        "compiler": context["compiler"],
        "compiler_version": context["compiler_version"],
        "compile_flags": context.get("compile_flags", []),
        "cmake_options": context.get("cmake_options", {}),
        "binary_checksum": context["binary_checksum"],
        "build_id": context["build_id"],
        "cpu_model": context["cpu_model"],
        "cpu_sockets": context["cpu_sockets"],
        "cpu_physical_cores": context["cpu_physical_cores"],
        "cpu_threads": context["cpu_threads"],
        "numa_nodes": context["numa_nodes"],
        "gpu_model": gpu_model,
        "gpu_id": gpu_id,
        "gpu_driver": gpu_driver,
        "gpu_runtime": gpu_runtime,
        "temperature": temperature,
        "timestep": timestep,
        "nstep": nstep,
        "measurement_profile": measurement_profile,
        "run_status": "SKIPPED",
        "timing_method": "NOT_MEASURED",
        "process_wall_seconds": None,
        "setup_seconds": None,
        "steady_step_seconds": None,
        "numerical_valid": False,
        "environment_valid": False,
        "quality_flags": [],
        "notes": notes,
    }
    schema_validate.validate_record(record)
    return record


def build_unsupported_precision_record(
    *, run_id, campaign_id, case, variant_id, size_id, sample_index, machine_id,
    gpu_backend, workload_metadata, context, temperature, timestep, nstep, measurement_profile,
    requested_precision="MIXED", workflow="ORDINARY_ASD", notes=None,
):
    """Section C: a `MIXED` (or otherwise `UNSUPPORTED_PRECISIONS`) request,
    recorded as the schema's first-class unsupported state -- never as a
    build/run failure, and never actually attempted (there is nothing to
    build). No GPU device is queried, matching an unsupported request never
    reaching a binary that could report one.
    """
    if requested_precision not in UNSUPPORTED_PRECISIONS:
        raise GpuCampaignError(
            f"requested_precision {requested_precision!r} is not one of "
            f"UNSUPPORTED_PRECISIONS {UNSUPPORTED_PRECISIONS!r}"
        )
    return _base_skipped_record(
        run_id=run_id, campaign_id=campaign_id, case_id=case.id, variant_id=variant_id, size_id=size_id,
        sample_index=sample_index, machine_id=machine_id,
        workflow=workflow, workload_class=case.record_workload_class,
        workload_metadata=workload_metadata, backend="GPU", gpu_backend=gpu_backend,
        omp_threads=context["omp_threads"], requested_precision=requested_precision,
        effective_cpu_precision=precision_audit.effective_cpu_precision(requested_precision),
        effective_gpu_precision=precision_audit.effective_gpu_precision(requested_precision, gpu_backend),
        precision_support_state="unsupported", context=context,
        gpu_model=None, gpu_id=None, gpu_driver=None, gpu_runtime=None,
        temperature=temperature, timestep=timestep, nstep=nstep, measurement_profile=measurement_profile,
        notes=notes or (
            f"{requested_precision} is deliberately not implemented in production "
            "(CMakeLists.txt fatal-errors at configure time); this is an unsupported "
            "mode, not a build failure."
        ),
    )


def build_unavailable_memory_record(
    *, run_id, campaign_id, case, variant_id, size_id, sample_index, machine_id,
    gpu_backend, requested_precision, workload_metadata, context, memory_classification, device,
    temperature, timestep, nstep, measurement_profile, workflow="ORDINARY_ASD", notes=None,
):
    """Section E: a build that genuinely supports `requested_precision`, on a
    size this campaign decided not to attempt because
    `gpu_memory.classify_gpu_memory_availability` predicted it would not fit
    the queried device. `memory_classification` is that function's return
    dict; `device` is the real `gpu_provenance` device query the estimate was
    compared against -- both numbers this section requires ("record
    requested system size and available GPU memory") are threaded straight
    into `notes` rather than silently reducing `size_id`.
    """
    if memory_classification["classification"] != UNAVAILABLE_MEMORY:
        raise GpuCampaignError(
            f"memory_classification does not report {UNAVAILABLE_MEMORY!r} "
            f"(got {memory_classification['classification']!r}); nothing to skip"
        )
    effective_gpu = precision_audit.effective_gpu_precision(requested_precision, gpu_backend)
    return _base_skipped_record(
        run_id=run_id, campaign_id=campaign_id, case_id=case.id, variant_id=variant_id, size_id=size_id,
        sample_index=sample_index, machine_id=machine_id,
        workflow=workflow, workload_class=case.record_workload_class,
        workload_metadata=workload_metadata, backend="GPU", gpu_backend=gpu_backend,
        omp_threads=context["omp_threads"], requested_precision=requested_precision,
        effective_cpu_precision=precision_audit.effective_cpu_precision(requested_precision),
        effective_gpu_precision=effective_gpu,
        precision_support_state="supported", context=context,
        gpu_model=device.get("model"), gpu_id=device.get("uuid") or str(device.get("index")),
        gpu_driver=device.get("driver_version"), gpu_runtime=context.get("gpu_runtime"),
        temperature=temperature, timestep=timestep, nstep=nstep, measurement_profile=measurement_profile,
        notes=notes or (
            f"{UNAVAILABLE_MEMORY}: requested natom={workload_metadata['natom']} "
            f"directed_interactions={workload_metadata.get('directed_interactions')} "
            f"estimated_required_mib={memory_classification['estimated_required_mib']:.1f} "
            f"available_mib={memory_classification['available_mib']:.1f} "
            f"safety_margin={memory_classification['safety_margin']}"
        ),
    )


def write_record(record, results_dir):
    """Write a GPU-campaign-built record with the same naming convention
    `runner.run_sample` uses (`{run_id}.json` under `results_dir`)."""
    results_dir = pathlib.Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    path = results_dir / f"{record['run_id']}.json"
    path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return path
