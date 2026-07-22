#!/usr/bin/env python3
"""Benchmark GPU real-space and FFT-convolution Hamiltonian backends.

Each benchmark gets a private fixture copy and is run once to warm up followed
by several timed repetitions.  The convolution run must print the device-backend
activation line, so a requested FFT benchmark cannot silently time a fallback.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[2]
INPUT_NAME = "inpsd.dat"
CONVOLUTION_ACTIVE = "Gpu: device lattice convolution active"
MEASUREMENT_PHASE_TIME = re.compile(
    r"^\s*Time for meas\. phase\s*:\s*([0-9.+\-EeDd]+)", re.MULTILINE
)
GPU_TIMING_ROW = re.compile(
    r"^\s*(?P<name>Total|hamiltonian)\s+"
    r"(?P<wall>[0-9.+\-EeDd]+)\s+ms wall\s+"
    r"(?P<gpu>[0-9.+\-EeDd]+)\s+ms GPU",
    re.MULTILINE,
)

# The BCC cases expose the small/large-system crossover; FeCo also exercises a
# two-atom basis. dhcpNd adds a four-site basis with long-range exchange
# (1,340 neighbours/atom). The SC cases use the RKKY generator in their fixture
# for a controlled system-size sweep at fixed neighbour count.
DEFAULT_MATRIX = (
    {"name": "bcc_small", "fixture": "tests/bccFe_cuda", "ncell": (16, 16, 16), "basis": 2, "steps": 1_000},
    {"name": "bcc_medium", "fixture": "tests/bccFe_cuda", "ncell": (32, 32, 32), "basis": 2, "steps": 1_000},
    {"name": "bcc_medium_long", "fixture": "tests/bccFe_cuda", "ncell": (32, 32, 32), "basis": 2, "steps": 10_000},
    {"name": "feco_medium", "fixture": "tests/FeCo_cuda", "ncell": (20, 20, 20), "basis": 2, "steps": 1_000},
    {"name": "dhcp_nd_long_range", "fixture": "examples/SimpleSystems/dhcpNd", "ncell": (32, 32, 32), "basis": 4, "steps": 1000},
    # Hold RKKY range fixed so this is a pure system-size sweep of the sparse
    # and convolution backends (rmax=5 gives 514 neighbours per atom).
    {"name": "sc_rkky_16", "fixture": "tests/SC_conv", "ncell": (16, 16, 16), "basis": 1,
     "steps": 5_000, "rkky_rmax": 5, "neighbours_per_atom": 514},
    {"name": "sc_rkky_40", "fixture": "tests/SC_conv", "ncell": (40, 40, 40), "basis": 1,
     "steps": 5_000, "rkky_rmax": 5, "neighbours_per_atom": 514},
    {"name": "sc_rkky_100", "fixture": "tests/SC_conv", "ncell": (100, 100, 100), "basis": 1,
     "steps": 5_000, "rkky_rmax": 5, "neighbours_per_atom": 514},
    # At fixed size, vary rmax to expose the neighbour-count crossover. The
    # rmax=5 point is sc_rkky_100 above and is shared by both sweeps.
    {"name": "sc_rkky_100_r3", "fixture": "tests/SC_conv", "ncell": (100, 100, 100), "basis": 1,
     "steps": 5_000, "rkky_rmax": 3, "neighbours_per_atom": 122},
    {"name": "sc_rkky_100_r8", "fixture": "tests/SC_conv", "ncell": (100, 100, 100), "basis": 1,
     "steps": 5_000, "rkky_rmax": 8, "neighbours_per_atom": 2108},
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, help="GPU-enabled UppASD executable")
    parser.add_argument("--case", action="append", dest="selected", help="benchmark name; may be repeated")
    parser.add_argument("--repetitions", type=int, default=3, help="timed repetitions per backend (default: 3)")
    parser.add_argument("--steps", type=int, help="override Nstep for every selected benchmark")
    parser.add_argument("--workdir", type=Path, help="directory for retained run files")
    parser.add_argument("--keep", action="store_true", help="retain a temporary work directory")
    parser.add_argument("--output", type=Path, help="write JSON report to this path")
    parser.add_argument("--csv", type=Path, help="write flat CSV samples to this path")
    parser.add_argument("--list", action="store_true", help="list the default benchmark matrix and exit")
    return parser.parse_args()


def set_input_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*)({re.escape(key)})(\s+).*$", re.IGNORECASE | re.MULTILINE)
    text, count = pattern.subn(rf"\g<1>{key}\g<3>{value}", text)
    return text if count else f"{text.rstrip()}\n{key} {value}\n"


def gpu_identity() -> str:
    for command in (("nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"), ("rocminfo",)):
        try:
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, timeout=10)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    return "unavailable"


def prepare_run(case: dict, destination: Path, convolution: bool, steps: int) -> None:
    shutil.copytree(ROOT / case["fixture"], destination)
    if "rkky_rmax" in case:
        generator = destination / "rkky_gen.py"
        generated = subprocess.run(
            [sys.executable, str(generator), "--rmax", str(case["rkky_rmax"])],
            cwd=destination, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        if generated.returncode:
            raise RuntimeError(
                f"{destination}: RKKY generation failed ({generated.returncode}): "
                f"{generated.stderr.strip()}"
            )
        (destination / "jij").write_text(generated.stdout)
    path = destination / INPUT_NAME
    text = path.read_text()
    overrides = {
        "SDEalgh": "5", "do_gpu": "Y", "gpu_mode": "1", "do_gpu_llg": "Y",
        "do_gpu_measurements": "N", "do_gpu_correlations": "N", "temp": "0.0",
        "Nstep": str(steps), "plotenergy": "0", "do_avrg": "N", "do_cumu": "N",
        "do_tottraj": "N", "tseed": "5", "do_gpu_convolution": "Y" if convolution else "N",
        "do_gpu_timings": "Y",
        # FFT fields are defined over the unit-cell basis, so the convolution
        # backend requires nHam (NH) to equal NA rather than the full atom count.
        # Use that same Hamiltonian representation for the real-space baseline.
        "do_reduced": "Y",
        "ncell": " ".join(map(str, case["ncell"])),
    }
    for key, value in overrides.items():
        text = set_input_value(text, key, value)
    # Do not append initial-phase keys: dhcp-Nd uses ip_mcanneal to allocate
    # ipTemp, so an appended scalar ip_temp would allocate it a second time.
    # When a fixture already defines an initial phase, disable it for timing.
    for key in ("ip_nphase", "ip_mcanneal"):
        if re.search(rf"^\s*{re.escape(key)}\s+", text, re.IGNORECASE | re.MULTILINE):
            text = set_input_value(text, key, "0")
    path.write_text(text)


def execute(binary: Path, directory: Path, require_convolution: bool, announce: bool = False) -> dict:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    started = time.perf_counter()
    completed = subprocess.run([str(binary)], cwd=directory, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    (directory / "run.log").write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(f"{directory} failed ({completed.returncode}); see run.log")
    if "Simulation finished" not in completed.stdout:
        raise RuntimeError(f"{directory} ended before simulation completion; see run.log")
    phase_match = MEASUREMENT_PHASE_TIME.search(completed.stdout)
    if phase_match is None:
        raise RuntimeError(f"{directory} did not report measurement-phase timing; see run.log")
    measurement_seconds = float(phase_match.group(1).replace("D", "E").replace("d", "e"))
    if measurement_seconds <= 0.0:
        raise RuntimeError(
            f"{directory} reported a zero measurement-phase time; increase --steps for this benchmark"
        )
    timing_start = completed.stdout.find("GPU measurement phase timing:")
    if timing_start < 0:
        raise RuntimeError(f"{directory} did not report detailed GPU timings; see run.log")
    gpu_times: dict[str, float] = {}
    for match in GPU_TIMING_ROW.finditer(completed.stdout[timing_start:]):
        gpu_times[match.group("name")] = float(match.group("gpu").replace("D", "E").replace("d", "e"))
    if "Total" not in gpu_times or "hamiltonian" not in gpu_times:
        raise RuntimeError(f"{directory} detailed GPU timings lack Total or hamiltonian; see run.log")
    if require_convolution:
        backend_lines = [line for line in completed.stdout.splitlines()
                         if line.startswith("Gpu: device lattice convolution")]
        if CONVOLUTION_ACTIVE not in completed.stdout:
            diagnostic = backend_lines[-1] if backend_lines else "no convolution backend diagnostic was printed"
            raise RuntimeError(f"{directory}: convolution requested but the device backend was not activated: {diagnostic}; see run.log")
        if announce:
            print(f"  {backend_lines[-1]}")
    return {"wall_seconds": elapsed, "measurement_seconds": measurement_seconds,
            "measurement_gpu_ms": gpu_times["Total"], "hamiltonian_gpu_ms": gpu_times["hamiltonian"]}


def summarize(case: dict, backend: str, samples: list[dict], steps: int) -> dict:
    atoms = case["basis"] * case["ncell"][0] * case["ncell"][1] * case["ncell"][2]
    wall_samples = [sample["wall_seconds"] for sample in samples]
    phase_samples = [sample["measurement_seconds"] for sample in samples]
    total_gpu_samples = [sample["measurement_gpu_ms"] for sample in samples]
    hamiltonian_gpu_samples = [sample["hamiltonian_gpu_ms"] for sample in samples]
    median_wall = statistics.median(wall_samples)
    median_phase = statistics.median(phase_samples)
    median_total_gpu = statistics.median(total_gpu_samples)
    median_hamiltonian_gpu = statistics.median(hamiltonian_gpu_samples)
    return {"backend": backend, "wall_seconds": wall_samples,
            "measurement_phase_seconds": phase_samples,
            "measurement_gpu_ms": total_gpu_samples,
            "hamiltonian_gpu_ms": hamiltonian_gpu_samples,
            "median_wall_seconds": median_wall,
            "median_measurement_phase_seconds": median_phase,
            "median_measurement_gpu_ms": median_total_gpu,
            "median_hamiltonian_gpu_ms": median_hamiltonian_gpu,
            "min_measurement_phase_seconds": min(phase_samples),
            "max_measurement_phase_seconds": max(phase_samples),
            "measurement_gpu_ms_per_step": median_total_gpu / steps,
            "hamiltonian_gpu_ms_per_step": median_hamiltonian_gpu / steps,
            "steps_per_second": 1_000.0 * steps / median_total_gpu,
            "atom_steps_per_second": 1_000.0 * atoms * steps / median_total_gpu}


def write_csv(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("case", "backend", "sample", "wall_seconds",
                                                    "measurement_phase_seconds", "measurement_gpu_ms",
                                                    "hamiltonian_gpu_ms"))
        writer.writeheader()
        for case in report["cases"]:
            for result in case["results"]:
                for index, (wall_seconds, phase_seconds, total_gpu_ms, hamiltonian_gpu_ms) in enumerate(
                        zip(result["wall_seconds"], result["measurement_phase_seconds"],
                            result["measurement_gpu_ms"], result["hamiltonian_gpu_ms"]), start=1):
                    writer.writerow({"case": case["name"], "backend": result["backend"],
                                     "sample": index, "wall_seconds": f"{wall_seconds:.9f}",
                                     "measurement_phase_seconds": f"{phase_seconds:.9f}",
                                     "measurement_gpu_ms": f"{total_gpu_ms:.6f}",
                                     "hamiltonian_gpu_ms": f"{hamiltonian_gpu_ms:.6f}"})


def main() -> int:
    args = parse_args()
    if args.list:
        for case in DEFAULT_MATRIX:
            details = ""
            if "neighbours_per_atom" in case:
                details = (f" rmax={case['rkky_rmax']} "
                           f"mnn={case['neighbours_per_atom']}")
            print(f"{case['name']}: {case['fixture']} ncell={'x'.join(map(str, case['ncell']))} "
                  f"steps={case['steps']}{details}")
        return 0
    if args.binary is None:
        raise SystemExit("--binary is required unless --list is used")
    if args.repetitions < 3:
        raise SystemExit("--repetitions must be at least 3")
    if args.steps is not None and args.steps <= 0:
        raise SystemExit("--steps must be positive")
    selected = set(args.selected or [case["name"] for case in DEFAULT_MATRIX])
    unknown = selected - {case["name"] for case in DEFAULT_MATRIX}
    if unknown:
        raise SystemExit(f"unknown case(s): {', '.join(sorted(unknown))}")
    binary = args.binary.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"not an executable: {binary}")
    root = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="uppasd-gpu-bench-"))
    root.mkdir(parents=True, exist_ok=True)
    report = {"binary": str(binary), "gpu": gpu_identity(), "repetitions": args.repetitions, "cases": []}
    failed = False
    try:
        for case in DEFAULT_MATRIX:
            if case["name"] not in selected:
                continue
            steps = args.steps or case["steps"]
            entry = {"name": case["name"], "fixture": case["fixture"], "ncell": case["ncell"], "basis": case["basis"],
                     "steps": steps, "results": []}
            if "neighbours_per_atom" in case:
                entry["rkky_rmax"] = case["rkky_rmax"]
                entry["neighbours_per_atom"] = case["neighbours_per_atom"]
            details = ""
            if "neighbours_per_atom" in case:
                details = (f" rmax={case['rkky_rmax']} "
                           f"mnn={case['neighbours_per_atom']}")
            print(f"{case['name']}: ncell={'x'.join(map(str, case['ncell']))} "
                  f"steps={steps}{details}")
            for convolution, backend in ((False, "real_space"), (True, "convolution")):
                warmup = root / case["name"] / f"{backend}-warmup"
                prepare_run(case, warmup, convolution, steps)
                execute(binary, warmup, convolution, announce=convolution)
                samples = []
                for repetition in range(args.repetitions):
                    directory = root / case["name"] / f"{backend}-{repetition + 1}"
                    prepare_run(case, directory, convolution, steps)
                    samples.append(execute(binary, directory, convolution))
                result = summarize(case, backend, samples, steps)
                entry["results"].append(result)
                print(f"  {backend:11s} wall={result['median_wall_seconds']:.6f}s "
                      f"GPU total={result['median_measurement_gpu_ms']:.3f}ms "
                      f"Hamiltonian={result['median_hamiltonian_gpu_ms']:.3f}ms "
                      f"steps/s={result['steps_per_second']:.1f}")
            real_space, convolution = entry["results"]
            entry["convolution_speedup"] = (
                real_space["median_measurement_gpu_ms"] /
                convolution["median_measurement_gpu_ms"]
            )
            entry["convolution_hamiltonian_speedup"] = (
                real_space["median_hamiltonian_gpu_ms"] /
                convolution["median_hamiltonian_gpu_ms"]
            )
            print(f"  convolution speedup: total={entry['convolution_speedup']:.3f}x "
                  f"hamiltonian={entry['convolution_hamiltonian_speedup']:.3f}x")
            report["cases"].append(entry)
        print(f"GPU: {report['gpu']}")
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
        if args.csv:
            write_csv(args.csv, report)
        return 0
    except Exception as error:
        failed = True
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    finally:
        if args.workdir or args.keep or failed:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
