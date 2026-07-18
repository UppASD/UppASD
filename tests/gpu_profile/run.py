#!/usr/bin/env python3
"""Benchmark the deterministic GPU spin-dynamics BCC fixture.

Every warm-up and measured repetition is a separate copy of tests/bccFe_cuda,
so no generated UppASD output is ever written into the repository fixture.
"""
from __future__ import annotations

import argparse
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
FIXTURE = ROOT / "tests" / "bccFe_cuda"
INPUT_NAME = "inpsd.dat"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="GPU-enabled UppASD executable")
    parser.add_argument("--steps", type=int, default=10_000, help="SD steps per run (default: 10000)")
    parser.add_argument("--repetitions", type=int, default=3, help="measured warm repetitions per mode")
    parser.add_argument("--ncell", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                        help="override BCC cell count")
    parser.add_argument("--convolution", nargs="?", const="Y", choices=("Y", "N", "both"), default="both",
                        help="benchmark convolution Y, ordinary N, or both (default: both; bare flag means Y)")
    parser.add_argument("--workdir", type=Path, help="directory in which to retain per-run fixture copies")
    parser.add_argument("--keep", action="store_true", help="retain the temporary work directory")
    parser.add_argument("--output", type=Path, help="write profile-summary.json here")
    return parser.parse_args()


def set_value(text: str, key: str, value: str) -> str:
    pattern = re.compile(rf"^(\s*)({re.escape(key)})(\s+).*$", re.IGNORECASE | re.MULTILINE)
    text, count = pattern.subn(rf"\g<1>{key}\g<3>{value}", text)
    return text if count else f"{text.rstrip()}\n{key} {value}\n"


def gpu_identity() -> str:
    for command in (("nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"),
                    ("rocminfo",)):
        try:
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL, timeout=10)
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    return "unavailable"


def prepare(directory: Path, args: argparse.Namespace, convolution: bool) -> None:
    shutil.copytree(FIXTURE, directory)
    path = directory / INPUT_NAME
    text = path.read_text()
    overrides = {
        "SDEalgh": "5", "do_gpu": "Y", "do_gpu_llg": "Y",
        "do_gpu_measurements": "N", "do_gpu_correlations": "N",
        "temp": "0.0", "ip_temp": "0.0", "Nstep": str(args.steps),
        "plotenergy": "0", "do_avrg": "N", "do_cumu": "N",
        "do_tottraj": "N", "tseed": "5",
        "do_gpu_convolution": "Y" if convolution else "N",
    }
    if args.ncell:
        overrides["ncell"] = " ".join(map(str, args.ncell))
    for key, value in overrides.items():
        text = set_value(text, key, value)
    path.write_text(text)


def execute(binary: Path, directory: Path) -> float:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    started = time.perf_counter()
    result = subprocess.run([str(binary)], cwd=directory, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - started
    (directory / "run.log").write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"{directory} failed ({result.returncode}); see run.log")
    return elapsed


def main() -> int:
    args = parse_args()
    if args.steps <= 0 or args.repetitions < 3:
        raise SystemExit("--steps must be positive and --repetitions must be at least 3")
    binary = args.binary.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"not an executable: {binary}")
    root = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="uppasd-gpu-profile-"))
    root.mkdir(parents=True, exist_ok=True)
    modes = [False, True] if args.convolution == "both" else [args.convolution == "Y"]
    report = {"command": [str(binary)], "binary": str(binary), "steps": args.steps,
              "repetitions": args.repetitions, "ncell": args.ncell,
              "gpu": gpu_identity(), "modes": []}
    print(f"Configuration: binary={binary} steps={args.steps} repetitions={args.repetitions} "
          f"ncell={args.ncell or 'fixture default'} modes={args.convolution}")
    try:
        for convolution in modes:
            label = "convolution" if convolution else "ordinary"
            warmup = root / f"{label}-warmup"
            prepare(warmup, args, convolution)
            execute(binary, warmup)
            samples = []
            for repetition in range(args.repetitions):
                directory = root / f"{label}-{repetition + 1}"
                prepare(directory, args, convolution)
                samples.append(execute(binary, directory))
            median = statistics.median(samples)
            spread = max(samples) - min(samples)
            result = {"convolution": convolution, "seconds": samples,
                      "median_seconds": median, "spread_seconds": spread}
            report["modes"].append(result)
            print(f"{label}: median={median:.6f}s spread={spread:.6f}s samples={samples}")
        print(f"GPU: {report['gpu']}")
        if args.output:
            args.output.mkdir(parents=True, exist_ok=True)
            (args.output / "profile-summary.json").write_text(json.dumps(report, indent=2) + "\n")
        return 0
    finally:
        if args.workdir or args.keep:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
