#!/usr/bin/env python3
"""CPU-HAM-05 production backend crossover campaign.

This driver intentionally stays independent of the optional benchmark-harness
Python dependencies. It runs the real CPU executable in isolated copies of
the admitted production templates, varying only backend-dispatch keywords.
The output is a raw JSONL stream and a derived JSON summary; both are useful
for an external run and neither is part of the tracked benchmark contract.

The timed quantity is the complete production process wall time. Fixed cost
and steady-state cost are fitted from two or more complete runs at different
``nstep`` values. The production executable does not export a pair-only
timer, so this driver records that limitation and does not manufacture one.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import platform
import re
import shutil
import statistics
import subprocess
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]

BACKEND_OVERRIDES = {
    "DIRECT": {"do_sparse": "N", "do_convolution": "N", "do_reduced": "N"},
    "REDUCED-DIRECT": {"do_sparse": "N", "do_convolution": "N", "do_reduced": "Y"},
    "SPARSE": {"do_sparse": "Y", "do_convolution": "N", "do_reduced": "N"},
    "CONVOLUTION": {"do_sparse": "N", "do_convolution": "Y", "do_reduced": "Y"},
}

# CPU-HAM-05 has a narrower scalar-J scope than the general CPU campaign and
# must not accidentally grow when a new case is admitted elsewhere.
CASES = {
    "B01_bccFe": {
        "template": ROOT / "benchmarks/cases/B01_bccFe/template",
        "variant": "bcc_fe_t0",
        "sizes": {
            "13x13x13": (13, 13, 13),
            "20x20x20": (20, 20, 20),
            "32x32x32": (32, 32, 32),
        },
        "basis": 2,
        "directed_per_atom": 96,
        "mean_neighbors": 96.0,
        "max_neighbors": 96,
        "ensembles": 10,
        # B01's anisotropy and symmetry-expanded input are a deliberate
        # negative control for the reduced/convolution eligibility gate.
        "backends": ("DIRECT", "SPARSE"),
        "ineligible": {"REDUCED-DIRECT": "B01 is not an eligible reduced scalar-J Hamiltonian"},
    },
    "B04_dhcpNd": {
        "template": ROOT / "benchmarks/cases/B04_dhcpNd/template",
        "variant": "dhcpNd_t0",
        "sizes": {
            "16x16x16": (16, 16, 16),
            "20x20x20": (20, 20, 20),
            "25x25x25": (25, 25, 25),
        },
        "basis": 4,
        "directed_per_atom": 1338,
        "mean_neighbors": 1338.0,
        "max_neighbors": 1340,
        "ensembles": 1,
        "backends": ("DIRECT", "REDUCED-DIRECT", "SPARSE", "CONVOLUTION"),
        "ineligible": {},
    },
    "B06_shortRangeScalarJ": {
        "template": ROOT / "benchmarks/cpu_ham05/short_range_scalar_control",
        "variant": "short_scalar_t0",
        "sizes": {
            "16x16x16": (16, 16, 16),
            "32x32x32": (32, 32, 32),
            "64x64x64": (64, 64, 64),
        },
        "basis": 1,
        "directed_per_atom": 6,
        "mean_neighbors": 6.0,
        "max_neighbors": 6,
        "ensembles": 1,
        "backends": ("DIRECT", "REDUCED-DIRECT", "SPARSE", "CONVOLUTION"),
        "ineligible": {},
    },
}

_TIME_RE = re.compile(
    r"Time for (initialization|meas\. phase|one meas\. iter)\s*:\s*([0-9.Ee+-]+)"
)


class CampaignError(RuntimeError):
    pass


def _format_value(value):
    if isinstance(value, (list, tuple)):
        return " ".join(str(item) for item in value)
    return str(value)


def _keyword(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    return stripped.split(None, 1)[0]


def apply_overrides(input_path, overrides):
    """Replace input keywords case-insensitively in one generated run only."""
    input_path = pathlib.Path(input_path)
    lines = input_path.read_text().splitlines()
    remaining = {key.lower(): value for key, value in overrides.items()}
    output = []
    for line in lines:
        key = _keyword(line)
        lower = key.lower() if key else None
        if lower in remaining:
            output.append(f"{key} {_format_value(remaining.pop(lower))}")
        else:
            output.append(line)
    for key, value in remaining.items():
        output.append(f"{key} {_format_value(value)}")
    input_path.write_text("\n".join(output) + "\n")


def _parse_phase_times(stdout):
    values = {}
    for name, value in _TIME_RE.findall(stdout):
        values[name] = float(value)
    return {
        "initialization_seconds_reported": values.get("initialization"),
        "measurement_phase_seconds_reported": values.get("meas. phase"),
        "one_measurement_iteration_seconds_reported": values.get("one meas. iter"),
    }


def _fit(points):
    if len(points) < 2:
        raise CampaignError("at least two complete nstep points are required for a fit")
    xs = [float(point[0]) for point in points]
    ys = [float(point[1]) for point in points]
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    denominator = sum((x - mean_x) ** 2 for x in xs)
    if denominator <= 0.0:
        raise CampaignError(f"nstep points are not distinct: {points!r}")
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denominator
    intercept = mean_y - slope * mean_x
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    rmse = math.sqrt(statistics.fmean([residual * residual for residual in residuals]))
    return {
        "setup_seconds": intercept,
        "steady_step_seconds": slope,
        "fit_rmse_seconds": rmse,
        "fit_points": [[int(x), y] for x, y in points],
    }


def _parse_threads(value):
    threads = [int(token) for token in value.split(",") if token.strip()]
    if not threads or any(thread < 1 for thread in threads):
        raise argparse.ArgumentTypeError("threads must be a comma-separated list of positive integers")
    return sorted(set(threads))


def _parse_steps(value):
    steps = [int(token) for token in value.split(",") if token.strip()]
    if len(steps) < 2 or any(step < 1 for step in steps):
        raise argparse.ArgumentTypeError("nsteps must contain at least two positive integers")
    return sorted(set(steps))


def _metadata(case, replication):
    spec = CASES[case]
    natom = spec["basis"] * replication[0] * replication[1] * replication[2]
    return {
        "natom": natom,
        "directed_interactions": natom * spec["directed_per_atom"],
        "mean_neighbors": spec["mean_neighbors"],
        "max_neighbors": spec["max_neighbors"],
        "basis_atoms": spec["basis"],
        "ensembles": spec["ensembles"],
    }


def _active_backend(backend, stdout):
    if backend == "DIRECT":
        return "DIRECT", True, None
    if backend == "REDUCED-DIRECT":
        active = "Validated reduced scalar-J stencil available" in stdout
        return backend if active else "DIRECT", active, None if active else "reduced stencil declined"
    if backend == "SPARSE":
        active = "Persistent scalar-J sparse backend ready:" in stdout
        return backend if active else "DIRECT", active, None if active else "sparse backend declined"
    if backend == "CONVOLUTION":
        active = "Persistent scalar-J CPU convolution ready:" in stdout
        return backend if active else "DIRECT", active, None if active else "convolution backend declined"
    raise CampaignError(f"unknown backend {backend!r}")


def _run_one(binary, run_dir, env, timeout):
    started = time.perf_counter()
    try:
        completed = subprocess.run(
            [str(binary)], cwd=run_dir, env=env, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "returncode": None,
            "timed_out": True,
            "wall_seconds": time.perf_counter() - started,
            "stdout": exc.stdout or "",
            "stderr": exc.stderr or "",
        }
    return {
        "returncode": completed.returncode,
        "timed_out": False,
        "wall_seconds": time.perf_counter() - started,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def _write_log(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def run_campaign(args):
    binary = pathlib.Path(args.binary).resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise CampaignError(f"CPU executable is not executable: {binary}")
    work_root = pathlib.Path(args.work_root).resolve()
    results_dir = pathlib.Path(args.results_dir).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    raw_path = results_dir / "raw.jsonl"
    summary_path = results_dir / "summary.json"
    if raw_path.exists():
        raw_path.unlink()

    selected_cases = args.cases or list(CASES)
    for case in selected_cases:
        if case not in CASES:
            raise CampaignError(f"unknown case {case!r}; choose from {sorted(CASES)}")
    selected_sizes = set(args.sizes or [])
    selected_backends = set(args.backends or [])
    for backend in selected_backends:
        if backend not in BACKEND_OVERRIDES:
            raise CampaignError(f"unknown backend {backend!r}")

    common_env = os.environ.copy()
    common_env.update({"OMP_DYNAMIC": "FALSE", "OMP_PLACES": "cores", "OMP_PROC_BIND": "close"})

    records = []
    with raw_path.open("a") as raw_file:
        for case in selected_cases:
            spec = CASES[case]
            for size_id, replication in spec["sizes"].items():
                if selected_sizes and size_id not in selected_sizes:
                    continue
                metadata = _metadata(case, replication)
                for backend in spec["backends"]:
                    if selected_backends and backend not in selected_backends:
                        continue
                    for threads in args.threads:
                        points_by_sample = []
                        point_records = []
                        for sample in range(args.samples):
                            sample_points = []
                            for nstep in args.nsteps:
                                run_id = f"{case}__{size_id}__{backend}__t{threads}__s{sample}__n{nstep}"
                                run_dir = work_root / run_id
                                if run_dir.exists():
                                    shutil.rmtree(run_dir)
                                shutil.copytree(spec["template"], run_dir)
                                apply_overrides(run_dir / "inpsd.dat", {
                                    "ncell": list(replication),
                                    "temp": 0.0,
                                    "tseed": 1,
                                    "do_prnstruct": 0,
                                    "nstep": nstep,
                                    "Nstep": nstep,
                                    "avrg_step": nstep + 1,
                                    "cumu_step": nstep + 1,
                                    "tottraj_step": nstep + 1,
                                    "ene_step": nstep + 1,
                                    **BACKEND_OVERRIDES[backend],
                                })
                                env = dict(common_env)
                                env["OMP_NUM_THREADS"] = str(threads)
                                result = _run_one(binary, run_dir, env, args.timeout)
                                stdout = result["stdout"]
                                stderr = result["stderr"]
                                active_backend, active, activation_reason = _active_backend(backend, stdout)
                                record = {
                                    "case": case,
                                    "size_id": size_id,
                                    "replication": list(replication),
                                    "backend_requested": backend,
                                    "backend_active": active_backend,
                                    "backend_active_verified": active,
                                    "backend_activation_reason": activation_reason,
                                    "omp_threads": threads,
                                    "sample": sample,
                                    "nstep": nstep,
                                    **metadata,
                                    "fft_grid": list(replication) if backend == "CONVOLUTION" else None,
                                    "fft_provider": "FFTW" if backend == "CONVOLUTION" else None,
                                    "returncode": result["returncode"],
                                    "timed_out": result["timed_out"],
                                    "completed": result["returncode"] == 0 and not result["timed_out"] and "Simulation finished" in stdout,
                                    "process_wall_seconds": result["wall_seconds"],
                                    "pair_field_time_seconds": None,
                                    "pair_field_time_status": "not exported by production executable",
                                    **_parse_phase_times(stdout),
                                }
                                if args.keep_logs:
                                    _write_log(results_dir / "logs" / f"{run_id}.stdout", stdout)
                                    _write_log(results_dir / "logs" / f"{run_id}.stderr", stderr)
                                raw_file.write(json.dumps(record, sort_keys=True) + "\n")
                                raw_file.flush()
                                records.append(record)
                                point_records.append(record)
                                if not record["completed"] or not active:
                                    raise CampaignError(
                                        f"run failed or backend activation was not verified: {run_id}; "
                                        f"returncode={record['returncode']} reason={activation_reason!r}"
                                    )
                                sample_points.append((nstep, result["wall_seconds"]))
                                shutil.rmtree(run_dir, ignore_errors=True)
                            points_by_sample.append(sample_points)

                        sample_fits = [_fit(points) for points in points_by_sample]
                        setup = statistics.median(item["setup_seconds"] for item in sample_fits)
                        steady = statistics.median(item["steady_step_seconds"] for item in sample_fits)
                        rmse = statistics.median(item["fit_rmse_seconds"] for item in sample_fits)
                        if steady <= 0.0:
                            raise CampaignError(
                                f"non-positive steady-state slope for {case}/{size_id}/{backend}/t{threads}: {sample_fits}"
                            )
                        summary = {
                            "case": case,
                            "size_id": size_id,
                            "replication": list(replication),
                            "backend": backend,
                            "backend_active": all(item["backend_active_verified"] for item in point_records),
                            "omp_threads": threads,
                            "sample_count": args.samples,
                            **metadata,
                            "fft_grid": list(replication) if backend == "CONVOLUTION" else None,
                            "fft_provider": "FFTW" if backend == "CONVOLUTION" else None,
                            "setup_seconds": setup,
                            "steady_step_seconds": steady,
                            "fit_rmse_seconds": rmse,
                            "full_effective_field_time_estimate_seconds": steady / 2.0,
                            "pair_field_time_seconds": None,
                            "pair_field_time_status": "not exported by production executable; see full-step/2 estimate",
                            "full_effective_field_time_status": "estimate: two effective-field calls per production step",
                            "spin_steps_per_second": 1.0 / steady,
                            "interaction_million_per_second": (
                                2.0 * metadata["directed_interactions"] * metadata["ensembles"] / steady
                            ) / 1.0e6,
                            "sample_fits": sample_fits,
                        }
                        records.append({"kind": "summary", **summary})

    summaries = [record for record in records if record.get("kind") == "summary"]
    payload = {
        "campaign": "CPU-HAM-05",
        "source_revision": os.environ.get("CPU_HAM05_SOURCE_REVISION", "unknown-at-runtime"),
        "host": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "processor": platform.processor(),
        },
        "binary": str(binary),
        "measurement": {
            "complete_production_process": True,
            "omp_dynamic": "FALSE",
            "omp_places": "cores",
            "omp_proc_bind": "close",
            "nsteps": args.nsteps,
            "samples": args.samples,
            "pair_field_timer": "unavailable from production executable",
        },
        "summaries": summaries,
        "ineligible": {
            case: spec["ineligible"] for case, spec in CASES.items()
            if case in selected_cases and spec["ineligible"]
        },
    }
    summary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"raw": str(raw_path), "summary": str(summary_path), "summary_count": len(summaries)}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, help="release CPU production executable")
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--threads", type=_parse_threads, default=[1, 2, 4, 8])
    parser.add_argument("--nsteps", type=_parse_steps, default=[5, 10, 20])
    parser.add_argument("--samples", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=900.0)
    parser.add_argument("--case", dest="cases", action="append")
    parser.add_argument("--size", dest="sizes", action="append")
    parser.add_argument("--backend", dest="backends", action="append")
    parser.add_argument("--keep-logs", action="store_true")
    args = parser.parse_args(argv)
    if args.samples < 1:
        parser.error("--samples must be positive")
    try:
        run_campaign(args)
    except CampaignError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
