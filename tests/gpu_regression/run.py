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


def set_input_value(text: str, key: str, value: str, before_key: str | None = None) -> str:
    """Replace active input records for key, or insert/append an input record."""
    pattern = re.compile(rf"^(\s*)({re.escape(key)})(\s+).*$", re.IGNORECASE | re.MULTILINE)
    replacement = rf"\g<1>{key}\g<3>{value}"
    text, count = pattern.subn(replacement, text)
    if count:
        return text
    if before_key:
        before_pattern = re.compile(rf"^(\s*)({re.escape(before_key)})(\s+).*$", re.IGNORECASE | re.MULTILINE)
        match = before_pattern.search(text)
        if match:
            return f"{text[:match.start()]}{key} {value}\n{text[match.start():]}"
    return f"{text.rstrip()}\n{key} {value}\n"


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
        "plotenergy": "1",
        "do_gpu": "Y" if gpu else "N",
        "do_gpu_llg": "Y",
        "do_gpu_measurements": "N",
        "do_gpu_correlations": "N",
    }
    if overrides.get("mode", "").upper() == "M":
        controlled["mcnstep"] = steps
        if gpu:
            controlled["do_gpu_mc"] = "Y"
    for key, value in {**controlled, **overrides}.items():
        text = set_input_value(text, key, str(value), before_key="Nstep" if key.lower() == "mcnstep" else None)
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


def numeric_table(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        try:
            rows.append([float(token.replace("D", "E").replace("d", "e")) for token in line.split()])
        except ValueError:
            continue
    return rows


def energy_rows_by_iteration(path: Path, label: str) -> dict[int, list[list[float]]]:
    rows: dict[int, list[list[float]]] = {}
    for row in numeric_table(path):
        if len(row) < 2:
            continue
        iteration = int(row[0])
        rows.setdefault(iteration, []).append(row[1:])
    if not rows:
        raise AssertionError(f"{label}: no numeric energy rows")
    return rows


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
        cpu_rows = energy_rows_by_iteration(
            single_file(cpu_dir, "totenergy.*.out", case["name"]), f"{case['name']} CPU")
        gpu_rows = energy_rows_by_iteration(
            single_file(gpu_dir, "totenergy.*.out", case["name"]), f"{case['name']} GPU")
        common_iterations = sorted(cpu_rows.keys() & gpu_rows.keys())
        if not common_iterations:
            raise AssertionError(f"{case['name']}: CPU and GPU have no energy iterations in common")
        cpu = [
            value
            for iteration in common_iterations
            for row in cpu_rows[iteration][:min(len(cpu_rows[iteration]), len(gpu_rows[iteration]))]
            for value in row
        ]
        gpu = [
            value
            for iteration in common_iterations
            for row in gpu_rows[iteration][:min(len(cpu_rows[iteration]), len(gpu_rows[iteration]))]
            for value in row
        ]
        result.append(compare_values("energy (common iterations)", cpu, gpu, rtol, atol))
    return result


def assert_easy_axis(case: dict, directory: Path) -> Difference:
    """Check the final single-spin state has relaxed to either easy-axis pole."""
    threshold = float(case["minimum_abs_mz"])
    rows = numeric_table(single_file(directory, "restart.*.out", case["name"]))
    if not rows or len(rows[-1]) < 7:
        raise AssertionError(f"{case['name']}: restart output lacks a final moment row")
    mz = abs(rows[-1][6])
    if not math.isfinite(mz):
        raise AssertionError(f"{case['name']}: final m_z is non-finite ({rows[-1][6]!r})")
    if mz < threshold:
        raise AssertionError(f"{case['name']}: |m_z|={mz:.6f}, expected at least {threshold:.6f}")
    return Difference("easy-axis |m_z|", 1.0 - mz, 0.0, 1)


def assert_finite_restart(case: dict, directory: Path) -> Difference:
    """Reject NaN/Inf final moments without requiring matching RNG trajectories."""
    rows = numeric_table(single_file(directory, "restart.*.out", case["name"]))
    if not rows:
        raise AssertionError(f"{case['name']}: restart output has no moment rows")
    values = [value for row in rows for value in row[3:]]
    if not values or not all(math.isfinite(value) for value in values):
        raise AssertionError(f"{case['name']}: restart output contains non-finite moment data")
    return Difference("finite restart", 0.0, 0.0, len(values))


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
                if "minimum_abs_mz" in case:
                    diffs.extend((assert_easy_axis(case, cpu_dir), assert_easy_axis(case, gpu_dir)))
                if case.get("restart_finite"):
                    diffs.extend((assert_finite_restart(case, cpu_dir), assert_finite_restart(case, gpu_dir)))
                summary = "; ".join(f"{d.label}: n={d.count}, abs={d.max_abs:.2e}, rel={d.max_rel:.2e}" for d in diffs)
                print(f"PASS {case['name']}: {summary}")
            except Exception as error:  # Keep running so a GPU session gets a complete table.
                failed += 1
                print(f"FAIL {case['name']}: {error}")
        return 1 if failed else 0
    finally:
        if args.workdir or args.keep or failed:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
