"""Production executable timing harness (WP-03).

Times the real UppASD production executable end to end -- never a kernel
microbenchmark -- per the blueprint's core principle (docs/UppASD Production
Benchmark Framework.md section 2). Each sample:

1. generates an isolated work directory from an immutable case template
   (``cases.generate_run_directory``, WP-02);
2. invokes the real executable in that directory, capturing full wall-clock
   time and stdout/stderr separately;
3. checks the run against the section F validity rules;
4. attaches this cell's workload metadata (resolved once, see
   ``resolve_workload_metadata`` below -- it describes the workload that was
   *configured*, not something a possibly-crashed process reported, so it is
   present even on a FAILED sample);
5. writes one schema-valid raw sample record.

Provenance (git/build/compiler), hardware (CPU/GPU identity) and precision
audit are WP-04's job, not this module's: callers supply them as ``context``.
``developer_context`` below is a best-effort stand-in for local/self-test
use; it always raises the ``metadata_incomplete`` quality flag, so records it
produces can never be ``environment_valid`` -- exactly VALIDITY.md's
"developer runs remain usable with flags present" rule, never silently
promoted to authoritative evidence.
"""

from __future__ import annotations

import datetime
import hashlib
import json
import os
import pathlib
import re
import subprocess
import time

from harness import cases as cases_mod
from harness import measurement_profiles
from harness import schema_validate
from harness import workload_metadata as workload_metadata_mod

SCHEMA_VERSION = "1.0.0"

_COMPLETION_MARKER = "Simulation finished"
_FATAL_ERROR_RE = re.compile(r"^\s*ERROR:", re.MULTILINE)
_NAN_INF_RE = re.compile(r"\b(NaN|Infinity|-Infinity|Inf)\b")
_REPORTED_NATOM_RE = re.compile(r"Number of atoms\s+(\d+)")

_REQUIRED_CONTEXT_FIELDS = (
    "machine_id",
    "backend",
    "omp_threads",
    "requested_precision",
    "precision_support_state",
    "git_commit",
    "git_dirty",
    "compiler",
    "compiler_version",
    "binary_checksum",
    "build_id",
    "cpu_model",
    "cpu_sockets",
    "cpu_physical_cores",
    "cpu_threads",
    "numa_nodes",
)


class RunnerError(ValueError):
    """A harness usage/configuration error -- distinct from a sample that
    ran but failed validity (that is recorded, not raised)."""


class ExecutionResult:
    """Raw outcome of invoking the executable once."""

    def __init__(self, returncode, stdout, stderr, wall_seconds, timed_out):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr
        self.wall_seconds = wall_seconds
        self.timed_out = timed_out


class ValidityAssessment:
    def __init__(self, run_status, numerical_valid, reasons):
        self.run_status = run_status
        self.numerical_valid = numerical_valid
        self.reasons = reasons


class SampleResult:
    """One executed sample: its run directory and the record written for it."""

    def __init__(self, run_directory, record, record_path, execution):
        self.run_directory = run_directory
        self.record = record
        self.record_path = record_path
        self.execution = execution


def _parse_fortran_float(token):
    return float(token.strip().replace("D", "E").replace("d", "e"))


def _execute_binary(binary_path, cwd, *, env=None, timeout_seconds=None):
    """Invoke the real executable once. Never used to derive timing from a
    subset of the run -- this is the complete process, start to finish."""
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            [str(binary_path)],
            cwd=str(cwd),
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - started
        stdout = exc.stdout if isinstance(exc.stdout, str) else (exc.stdout or b"").decode(errors="replace")
        stderr = exc.stderr if isinstance(exc.stderr, str) else (exc.stderr or b"").decode(errors="replace")
        return ExecutionResult(returncode=None, stdout=stdout, stderr=stderr, wall_seconds=elapsed, timed_out=True)
    elapsed = time.perf_counter() - started
    return ExecutionResult(
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        wall_seconds=elapsed,
        timed_out=False,
    )


def _write_raw_artifacts(run_dir, execution, run_status):
    (run_dir / "stdout.log").write_text(execution.stdout)
    (run_dir / "stderr.log").write_text(execution.stderr)
    timing_record = {
        "run_status": run_status,
        "process_wall_seconds": execution.wall_seconds,
        "returncode": execution.returncode,
        "timed_out": execution.timed_out,
    }
    (run_dir / "timing.json").write_text(json.dumps(timing_record, indent=2, sort_keys=True) + "\n")


def assess_validity(execution, run_dir, *, expected_natom, essential_output_required=True):
    """Section F: decide run_status/numerical_valid from real run evidence.

    Every check is against something the run itself produced -- exit code,
    stdout/stderr text, the files it wrote -- never assumed.
    """
    reasons = []

    if execution.timed_out:
        return ValidityAssessment(run_status="TIMEOUT", numerical_valid=False, reasons=["process exceeded timeout"])

    if execution.returncode != 0:
        reasons.append(f"nonzero exit code {execution.returncode}")

    combined = execution.stdout + "\n" + execution.stderr
    fatal = _FATAL_ERROR_RE.search(combined)
    if fatal:
        reasons.append(f"fatal error detected: {fatal.group(0).strip()!r}")

    if _COMPLETION_MARKER not in execution.stdout:
        reasons.append(f"expected completion evidence ({_COMPLETION_MARKER!r}) missing from stdout")

    nan_inf = _NAN_INF_RE.search(combined)
    if nan_inf:
        reasons.append(f"NaN/Inf detected in process output: {nan_inf.group(0)!r}")

    reported_natom_match = _REPORTED_NATOM_RE.search(execution.stdout)
    if reported_natom_match is not None:
        reported_natom = int(reported_natom_match.group(1))
        if reported_natom != expected_natom:
            reasons.append(
                f"reported atom count ({reported_natom}) does not match the "
                f"expected generated case ({expected_natom})"
            )

    if essential_output_required and not reasons:
        try:
            simid = cases_mod.read_simid(run_dir)
        except cases_mod.CaseManifestError:
            reasons.append("essential output cannot be checked: inpsd.dat has no simid")
        else:
            restart_path = run_dir / f"restart.{simid}.out"
            if not restart_path.is_file():
                reasons.append(f"essential output absent: {restart_path.name} was not written")

    if reasons:
        return ValidityAssessment(run_status="FAILED", numerical_valid=False, reasons=reasons)
    return ValidityAssessment(run_status="COMPLETED", numerical_valid=True, reasons=[])


def resolve_workload_metadata(case, variant_id, size_id, work_root, binary_path, *, env=None,
                               timeout_seconds=120, probe_nstep=1, run_id=None):
    """Run one minimal probe to determine a (case, variant, size) cell's real
    workload metadata (WP-03 section A.8 / B).

    This is a workload property, not a sample outcome: it does not change
    between repeated samples of the same cell, so it is resolved once here
    and then attached to every sample's record for that cell -- including a
    later sample that fails, per the result contract (a raw sample record
    always carries real `natom`/workload fields, regardless of `run_status`).
    """
    size = case.resolve_size(size_id)
    overrides = {"Nstep": probe_nstep}
    method = case.manifest["workload_metadata_method"]
    if method == "neighbor_list_from_struct_output":
        if "do_prnstruct" not in case.manifest["allowed_input_overrides"]:
            raise RunnerError(
                f"case {case.id!r} uses workload_metadata_method "
                f"{method!r} but does not allow overriding do_prnstruct"
            )
        overrides["do_prnstruct"] = 1

    probe_run_id = run_id or f"probe__{case.id}__{variant_id}__{size_id}"
    run_directory = cases_mod.generate_run_directory(
        case, variant_id, size_id, work_root, run_id=probe_run_id, extra_overrides=overrides
    )
    execution = _execute_binary(binary_path, run_directory.path, env=env, timeout_seconds=timeout_seconds)
    _write_raw_artifacts(run_directory.path, execution, "PROBE")
    if execution.timed_out or execution.returncode != 0 or _COMPLETION_MARKER not in execution.stdout:
        raise RunnerError(
            f"workload-metadata probe failed for case={case.id!r} "
            f"variant={variant_id!r} size={size_id!r}; see {run_directory.path}"
        )
    return workload_metadata_mod.compute_workload_metadata(case, size, run_directory.path)


def developer_context(binary_path, *, machine_id="dev-local", backend="CPU", omp_threads=None):
    """Best-effort context for local/self-test runs (never authoritative).

    Every field WP-04 would actually audit (build toolchain, CPU/GPU
    identity, precision audit) is filled with an honest "not audited"
    placeholder rather than a guess, and `metadata_incomplete` is always
    raised so `environment_valid` can never come out true for a record built
    from this context.
    """
    binary_path = str(binary_path)
    digest = hashlib.sha256()
    with open(binary_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    binary_checksum = digest.hexdigest()

    git_commit = "0" * 40
    git_dirty = False
    git_dirty_files = []
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=10
        )
        if commit.returncode == 0:
            git_commit = commit.stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain"], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=10
        )
        if status.returncode == 0:
            git_dirty_files = [line[3:] for line in status.stdout.splitlines() if line.strip()]
            git_dirty = bool(git_dirty_files)
    except (OSError, subprocess.TimeoutExpired):
        pass

    threads = omp_threads or int(os.environ.get("OMP_NUM_THREADS", "1"))
    cpu_count = os.cpu_count() or 1

    return {
        "machine_id": machine_id,
        "backend": backend,
        "gpu_backend": None,
        "omp_threads": threads,
        "omp_places": os.environ.get("OMP_PLACES"),
        "omp_proc_bind": os.environ.get("OMP_PROC_BIND"),
        "omp_dynamic": None,
        "process_affinity": None,
        "requested_precision": "DOUBLE",
        "effective_cpu_precision": "UNKNOWN",
        "effective_gpu_precision": None,
        "precision_support_state": "unaudited",
        "comparison_precision_class": "UNAUDITED",
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "git_dirty_files": git_dirty_files,
        "compiler": "unaudited",
        "compiler_version": "unaudited",
        "compile_flags": [],
        "cmake_options": {},
        "binary_checksum": binary_checksum,
        "build_id": f"dev-{binary_checksum[:12]}",
        "cpu_model": "unaudited",
        "cpu_sockets": 1,
        "cpu_physical_cores": cpu_count,
        "cpu_threads": cpu_count,
        "numa_nodes": 1,
        "gpu_model": None,
        "gpu_id": None,
        "gpu_driver": None,
        "gpu_runtime": None,
        "dipole_backend": None,
        "quality_flags": ["metadata_incomplete"],
    }


def _check_context(context):
    missing = [field for field in _REQUIRED_CONTEXT_FIELDS if field not in context]
    if missing:
        raise RunnerError(f"context is missing required fields: {sorted(missing)}")


def run_sample(
    case,
    variant_id,
    size_id,
    *,
    nstep,
    run_id,
    campaign_id,
    sample_index,
    work_root,
    results_dir,
    binary_path,
    measurement_profile,
    workload_metadata,
    context,
    env=None,
    timeout_seconds=None,
    extra_overrides=None,
):
    """Execute one complete production sample and write its raw record.

    ``workload_metadata`` must come from :func:`resolve_workload_metadata`
    for this ``(case, variant_id, size_id)`` cell -- it is not recomputed
    per sample. ``context`` supplies the provenance/hardware/precision/
    run-configuration fields this WP does not itself audit (see module
    docstring); :func:`developer_context` is a non-authoritative stand-in.
    """
    _check_context(context)
    if measurement_profile not in measurement_profiles.PROFILES:
        raise RunnerError(f"unknown measurement_profile {measurement_profile!r}")
    extra_overrides = dict(extra_overrides or {})
    if "Nstep" in extra_overrides:
        raise RunnerError("pass nstep explicitly; do not also set it via extra_overrides")

    overrides = {"Nstep": nstep}
    overrides.update(measurement_profiles.profile_overrides(measurement_profile, nstep))
    overrides.update(extra_overrides)

    size = case.resolve_size(size_id)
    run_directory = cases_mod.generate_run_directory(
        case, variant_id, size_id, work_root, run_id=run_id, extra_overrides=overrides
    )
    run_dir = run_directory.path

    expected_natom = workload_metadata_mod.expected_natom(case, size, run_dir)
    if workload_metadata["natom"] != expected_natom:
        raise RunnerError(
            f"workload_metadata natom ({workload_metadata['natom']}) does not match "
            f"this cell's expected atom count ({expected_natom}); it was not resolved "
            f"for this (case, variant_id, size_id)"
        )

    started_utc = datetime.datetime.now(datetime.timezone.utc)
    execution = _execute_binary(binary_path, run_dir, env=env, timeout_seconds=timeout_seconds)

    assessment = assess_validity(execution, run_dir, expected_natom=expected_natom)
    _write_raw_artifacts(run_dir, execution, assessment.run_status)

    temp_token = cases_mod.read_keyword(run_dir, "temp")
    if temp_token is None:
        raise RunnerError(f"case {case.id!r}'s generated input has no 'temp' keyword to report")
    temperature = _parse_fortran_float(temp_token)
    timestep_token = cases_mod.read_keyword(run_dir, "timestep")
    if timestep_token is None:
        raise RunnerError(
            f"case {case.id!r}'s template does not declare a 'timestep' keyword; "
            "the result contract requires a real, non-null timestep for every record"
        )
    timestep = _parse_fortran_float(timestep_token)
    recorded_nstep = int(cases_mod.read_keyword(run_dir, "Nstep"))

    completed = assessment.run_status == "COMPLETED"
    quality_flags = set(context.get("quality_flags", []))
    if context["git_dirty"]:
        quality_flags.add("dirty_tree")
    environment_valid = completed and not quality_flags

    record = {
        "schema_version": SCHEMA_VERSION,
        "record_kind": "raw_sample",
        "run_id": run_id,
        "campaign_id": campaign_id,
        "case_id": case.id,
        "variant_id": variant_id,
        "size_id": size_id,
        "sample_index": sample_index,
        "machine_id": context["machine_id"],
        "timestamp_utc": started_utc.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "workflow": context.get("workflow", "ORDINARY_ASD"),
        "workload_class": case.record_workload_class,
        "natom": workload_metadata["natom"],
        "directed_interactions": workload_metadata.get("directed_interactions"),
        "mean_neighbors": workload_metadata.get("mean_neighbors"),
        "median_neighbors": workload_metadata.get("median_neighbors"),
        "max_neighbors": workload_metadata.get("max_neighbors"),
        "fft_grid": workload_metadata.get("fft_grid"),
        "fft_grid_padded": workload_metadata.get("fft_grid_padded"),
        "fft_grid_points": workload_metadata.get("fft_grid_points"),
        "dipole_backend": context.get("dipole_backend"),
        "backend": context["backend"],
        "gpu_backend": context.get("gpu_backend"),
        "omp_threads": context["omp_threads"],
        "omp_places": context.get("omp_places"),
        "omp_proc_bind": context.get("omp_proc_bind"),
        "omp_dynamic": context.get("omp_dynamic"),
        "process_affinity": context.get("process_affinity"),
        "requested_precision": context["requested_precision"],
        "effective_cpu_precision": context.get("effective_cpu_precision"),
        "effective_gpu_precision": context.get("effective_gpu_precision"),
        "precision_support_state": context["precision_support_state"],
        "comparison_precision_class": context.get("comparison_precision_class"),
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
        "gpu_model": context.get("gpu_model"),
        "gpu_id": context.get("gpu_id"),
        "gpu_driver": context.get("gpu_driver"),
        "gpu_runtime": context.get("gpu_runtime"),
        "temperature": temperature,
        "timestep": timestep,
        "nstep": recorded_nstep,
        "measurement_profile": measurement_profile,
        "run_status": assessment.run_status,
        "timing_method": "EXTERNAL_PROCESS_WALL_CLOCK" if completed else "NOT_MEASURED",
        "process_wall_seconds": execution.wall_seconds if completed else None,
        "setup_seconds": None,
        "steady_step_seconds": None,
        "numerical_valid": assessment.numerical_valid,
        "environment_valid": environment_valid,
        "quality_flags": sorted(quality_flags),
        "notes": "; ".join(assessment.reasons) if assessment.reasons else None,
    }

    schema_validate.validate_record(record)

    results_dir = pathlib.Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    record_path = results_dir / f"{run_id}.json"
    record_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")

    return SampleResult(run_directory=run_directory, record=record, record_path=record_path, execution=execution)
