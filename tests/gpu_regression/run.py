#!/usr/bin/env python3
"""Compare deterministic CPU and GPU UppASD runs on a GPU-capable machine.

The runner deliberately uses the GPU-enabled executable for both sides: `do_gpu N`
exercises the Fortran implementation and `do_gpu Y` exercises the device path.
Each run receives a private copy of its fixture, so generated files never modify the
repository's examples or regular regression fixtures.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CASES = Path(__file__).with_name("cases.json")
INPUT_NAME = "inpsd.dat"


@dataclass
class Difference:
    label: str
    max_abs: float
    max_rel: float
    count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True,
                        help="GPU-enabled UppASD executable")
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--case", action="append", dest="selected",
                        help="case name to run; may be repeated")
    parser.add_argument("--workdir", type=Path,
                        help="directory for retained run files (default: temporary)")
    parser.add_argument("--keep", action="store_true", help="retain a temporary work directory")
    parser.add_argument("--list", action="store_true", help="list cases and exit")
    return parser.parse_args()


def load_cases(path: Path) -> tuple[dict, list[dict]]:
    data = json.loads(path.read_text())
    return data.get("defaults", {}), data["cases"]


def set_input_value(text: str, key: str, value: str) -> str:
    """Replace active input records for key, or append an input record."""
    pattern = re.compile(rf"^(\s*)({re.escape(key)})(\s+).*$", re.IGNORECASE | re.MULTILINE)
    replacement = rf"\1{key}\3{value}"
    text, count = pattern.subn(replacement, text)
    return text if count else f"{text.rstrip()}\n{key} {value}\n"


def prepare_run(case: dict, defaults: dict, gpu: bool, destination: Path) -> None:
    fixture = ROOT / case["fixture"]
    shutil.copytree(fixture, destination)
    input_path = destination / INPUT_NAME
    if not input_path.exists():
        raise FileNotFoundError(f"{case['name']}: fixture has no {INPUT_NAME}: {fixture}")
    text = input_path.read_text()
    overrides = dict(case.get("overrides", {}))
    steps = str(case.get("steps", defaults.get("steps", 100)))
    # T=0 and a short fixed trajectory make the CPU/GPU comparison reproducible.
    controlled = {
        "temp": "0.0",
        "ip_temp": "0.0",
        "Nstep": steps,
        "mcnstep": steps,
        "plotenergy": "1",
        "do_gpu": "Y" if gpu else "N",
        "do_gpu_llg": "Y",
        "do_gpu_measurements": "N",
        "do_gpu_correlations": "N",
    }
    if gpu and overrides.get("mode", "").upper() == "M":
        controlled["do_gpu_mc"] = "Y"
    for key, value in {**controlled, **overrides}.items():
        text = set_input_value(text, key, str(value))
    input_path.write_text(text)


def run(binary: Path, directory: Path) -> None:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    completed = subprocess.run([str(binary)], cwd=directory, env=env,
                               text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (directory / "run.log").write_text(completed.stdout)
    if completed.returncode:
        raise RuntimeError(f"{directory.name} failed ({completed.returncode}); see {directory / 'run.log'}")


def numeric_rows(path: Path, columns: slice | None = None) -> list[float]:
    values: list[float] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            row = [float(token.replace("D", "E").replace("d", "e")) for token in line.split()]
        except ValueError:
            continue
        values.extend(row[columns] if columns else row)
    return values


def compare_values(label: str, cpu: Iterable[float], gpu: Iterable[float], rtol: float, atol: float) -> Difference:
    left, right = list(cpu), list(gpu)
    if len(left) != len(right):
        raise AssertionError(f"{label}: value-count mismatch: CPU={len(left)}, GPU={len(right)}")
    max_abs = max_rel = 0.0
    for index, (expected, actual) in enumerate(zip(left, right)):
        abs_error = abs(expected - actual)
        rel_error = abs_error / max(abs(expected), atol)
        max_abs, max_rel = max(max_abs, abs_error), max(max_rel, rel_error)
        if abs_error > atol + rtol * abs(expected):
            raise AssertionError(
                f"{label}: element {index}: CPU={expected:.17e}, GPU={actual:.17e}, "
                f"abs={abs_error:.3e}, rel={rel_error:.3e}")
    return Difference(label, max_abs, max_rel, len(left))


def single_file(directory: Path, pattern: str, label: str) -> Path:
    files = sorted(directory.glob(pattern))
    if len(files) != 1:
        raise FileNotFoundError(f"{label}: expected one {pattern} in {directory}, found {len(files)}")
    return files[0]


def compare_outputs(case: dict, defaults: dict, cpu_dir: Path, gpu_dir: Path) -> list[Difference]:
    rtol = float(case.get("rtol", defaults.get("rtol", 1e-10)))
    atol = float(case.get("atol", defaults.get("atol", 1e-12)))
    comparisons = case.get("compare", defaults.get("compare", ["restart", "energy"]))
    result: list[Difference] = []
    if "restart" in comparisons:
        # The first three fields are identical metadata; compare |m| and mx,my,mz only.
        cpu = numeric_rows(single_file(cpu_dir, "restart.*.out", case["name"]), slice(3, 7))
        gpu = numeric_rows(single_file(gpu_dir, "restart.*.out", case["name"]), slice(3, 7))
        result.append(compare_values("restart", cpu, gpu, rtol, atol))
    if "energy" in comparisons:
        cpu = numeric_rows(single_file(cpu_dir, "totenergy.*.out", case["name"]))
        gpu = numeric_rows(single_file(gpu_dir, "totenergy.*.out", case["name"]))
        result.append(compare_values("energy", cpu, gpu, rtol, atol))
    return result


def main() -> int:
    args = parse_args()
    defaults, cases = load_cases(args.cases)
    if args.list:
        print("\n".join(case["name"] for case in cases))
        return 0
    selected = set(args.selected or [case["name"] for case in cases])
    unknown = selected - {case["name"] for case in cases}
    if unknown:
        raise SystemExit(f"unknown case(s): {', '.join(sorted(unknown))}")
    binary = args.binary.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise SystemExit(f"not an executable: {binary}")
    root = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="uppasd-gpu-regression-"))
    root.mkdir(parents=True, exist_ok=True)
    failed = 0
    try:
        for case in cases:
            if case["name"] not in selected:
                continue
            cpu_dir, gpu_dir = root / f"{case['name']}-cpu", root / f"{case['name']}-gpu"
            try:
                prepare_run(case, defaults, False, cpu_dir)
                prepare_run(case, defaults, True, gpu_dir)
                run(binary, cpu_dir)
                run(binary, gpu_dir)
                diffs = compare_outputs(case, defaults, cpu_dir, gpu_dir)
                summary = "; ".join(f"{d.label}: n={d.count}, abs={d.max_abs:.2e}, rel={d.max_rel:.2e}" for d in diffs)
                print(f"PASS {case['name']}: {summary}")
            except Exception as error:  # Keep running so a GPU session gets a complete table.
                failed += 1
                print(f"FAIL {case['name']}: {error}")
        return 1 if failed else 0
    finally:
        if args.workdir or args.keep:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
