#!/usr/bin/env python3
"""Compare the shift_b trajectory in UppASD's C8 frame with pyswatter.

The UppASD trajectory path fixes one global frame at initialization.  The
moment file consumed by pyswatter is in Cartesian coordinates, so this helper
reproduces ``oam_init``'s frame construction, rotates the moments, runs
``spin-oam-balance`` with the matching origin, and compares the shared
origin-referenced observables.

The total-OAM and balance columns are reported by both programs but are not
used for pass/fail: UppASD defines its total OAM with the centroid-referenced
lambda, while pyswatter's balance output uses its origin-referenced lambda.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np


FORTRAN_COLUMNS = (
    "step",
    "lambda_L_origin",
    "lambda_L_centroid",
    "N_m",
    "Lz_tot_hbar",
    "dSz_hbar",
    "balance",
    "R_x",
    "R_y",
    "sigma_psi",
)


def numeric_rows(path: Path):
    """Yield non-comment rows from an UppASD moment/restart file."""

    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 7 and fields[0].lstrip("+-").isdigit():
            yield fields


def c8_frame(initial_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct ex, ey, ez using the same algorithm as ``oam_init``."""

    rows = list(numeric_rows(initial_path))
    if not rows:
        raise RuntimeError(f"no moment rows found in {initial_path}")

    initial_iteration = rows[0][0]
    rows = [row for row in rows if row[0] == initial_iteration]
    moments = np.array(
        [[float(row[4]), float(row[5]), float(row[6])] for row in rows],
        dtype=float,
    )
    mean_m = moments.mean(axis=0)
    mean_norm = np.linalg.norm(mean_m)
    if mean_norm <= 1.0e-14:
        raise RuntimeError("initial average moment is zero")

    ez = mean_m / mean_norm
    if np.min(moments @ ez) < 0.9:
        raise RuntimeError("initial state violates the C8 alignment guard")

    seed = np.array([1.0, 0.0, 0.0])
    if abs(seed @ ez) >= 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    ex = seed - (seed @ ez) * ez
    ex /= np.linalg.norm(ex)
    ey = np.cross(ez, ex)
    return ex, ey, ez


def rotate_moments(source: Path, destination: Path, ex, ey, ez) -> None:
    """Write a pyswatter-readable moment file in the C8 frame."""

    with destination.open("w") as output:
        for line in source.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                output.write(line + "\n")
                continue

            fields = line.split()
            if len(fields) < 7 or not fields[0].lstrip("+-").isdigit():
                output.write(line + "\n")
                continue

            moment = np.array([float(fields[4]), float(fields[5]), float(fields[6])])
            rotated = (moment @ ex, moment @ ey, moment @ ez)
            output.write(
                " ".join(fields[:4] + [f"{value:.16E}" for value in rotated] + fields[7:])
                + "\n"
            )


def read_fortran_row(path: Path) -> dict[str, float]:
    rows = []
    for line in path.read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            rows.append([float(value) for value in line.split()])
    if not rows:
        raise RuntimeError(f"no data rows found in {path}")
    return dict(zip(FORTRAN_COLUMNS, rows[0]))


def read_pyswatter_row(path: Path) -> dict[str, float]:
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream))
    return {key: float(value) for key, value in row.items() if value is not None}


def check_close(label: str, actual: float, expected: float, tolerance: float) -> bool:
    difference = actual - expected
    ok = math.isfinite(actual) and math.isfinite(expected) and abs(difference) <= tolerance
    status = "PASS" if ok else "FAIL"
    print(
        f"[{status}] {label}: pyswatter={actual:.15g} "
        f"Fortran={expected:.15g} diff={difference:+.3e} tol={tolerance:.3e}"
    )
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", nargs="?", default="traj_checks/shift_b")
    parser.add_argument("--pyswatter", default="pyswatter-animate")
    parser.add_argument("--rotated-moments", default="restart.oamtest.c8.out")
    parser.add_argument("--output", default="ref_c8.csv")
    args = parser.parse_args()

    case_dir = Path(args.case_dir).resolve()
    position_file = case_dir / "coord.oamtest.out"
    initial_file = case_dir / "restart.in"
    moment_file = case_dir / "restart.oamtest.out"
    rotated_file = case_dir / args.rotated_moments
    output_file = case_dir / args.output
    for required in (position_file, initial_file, moment_file):
        if not required.is_file():
            parser.error(f"missing required file: {required}")

    executable = shutil.which(args.pyswatter)
    if executable is None:
        parser.error(f"pyswatter executable not found: {args.pyswatter}")

    ex, ey, ez = c8_frame(initial_file)
    print("C8 frame:")
    print(f"  ex = {ex}")
    print(f"  ey = {ey}")
    print(f"  ez = {ez}")
    rotate_moments(moment_file, rotated_file, ex, ey, ez)
    print(f"Rotated moments: {rotated_file}")

    subprocess.run(
        [
            executable,
            "spin-oam-balance",
            position_file.name,
            rotated_file.name,
            "--lz-integration",
            "site",
            "--shift",
            "0",
            "0",
            "0",
            "--output",
            output_file.name,
        ],
        cwd=case_dir,
        check=True,
    )

    pyswatter = read_pyswatter_row(output_file)
    fortran = read_fortran_row(case_dir / "oam_traj.oamtest.out")
    checks = [
        check_close("lambda_L", pyswatter["lambda_L"], fortran["lambda_L_origin"], 1.0e-3),
        check_close("N_m", pyswatter["N_m"], fortran["N_m"], 2.0e-6),
        check_close(
            "delta_Sz_over_hbar",
            pyswatter["delta_Sz_over_hbar"],
            fortran["dSz_hbar"],
            2.0e-6,
        ),
    ]
    print(
        "Not compared: pyswatter Lz_total_over_hbar/balance_over_hbar versus "
        "Fortran Lz_tot_hbar/balance (origin versus centroid convention)."
    )
    return 0 if all(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
