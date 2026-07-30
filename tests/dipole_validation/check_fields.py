#!/usr/bin/env python3
"""Compare UppASD total effective fields with analytic dipole references.

Usage: run a generated finite case with an UppASD executable, then run this
script on that case directory.  It compares `befftot.<simid>.out` (Tesla) to
the independently generated `analytic_reference.json`; it never compares CPU
and GPU output to each other.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def input_value(path: Path, key: str) -> str:
    for line in path.read_text().splitlines():
        tokens = line.split()
        if tokens and tokens[0].lower() == key.lower():
            if len(tokens) < 2:
                raise ValueError(f"{path}: {key} has no value")
            return tokens[1]
    raise ValueError(f"{path}: missing {key}")


def read_field_blocks(path: Path, sites: int) -> dict[tuple[int, int], list[tuple[float, float, float]]]:
    blocks: dict[tuple[int, int], dict[int, tuple[float, float, float]]] = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = line.replace("D", "E").replace("d", "e").split()
        if len(values) != 7:
            continue
        iteration, site, replica = int(float(values[0])), int(values[1]), int(values[2])
        if not 1 <= site <= sites:
            raise AssertionError(f"{path}: invalid site {site}")
        blocks.setdefault((iteration, replica), {})[site] = tuple(map(float, values[3:6]))
    complete: dict[tuple[int, int], list[tuple[float, float, float]]] = {}
    for key, field in blocks.items():
        if len(field) == sites:
            complete[key] = [field[site] for site in range(1, sites + 1)]
    return complete


def read_dipole_energy(path: Path) -> tuple[int, float]:
    """Read the first `Dip` value from a standard UppASD totenergy table."""
    columns: list[str] | None = None
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("#"):
            candidate = line.lstrip()[1:].split()
            if "Dip" in candidate:
                columns = candidate
            continue
        if not line.strip() or columns is None:
            continue
        values = line.replace("D", "E").replace("d", "e").split()
        if len(values) != len(columns):
            continue
        return int(float(values[0])), float(values[columns.index("Dip")])
    raise AssertionError(f"{path}: no numeric Dip column")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="generated case directory")
    parser.add_argument("--rtol", type=float, default=2e-8)
    parser.add_argument("--atol", type=float, default=2e-12)
    parser.add_argument("--energy-rtol", type=float, default=2e-8)
    parser.add_argument("--energy-atol", type=float, default=2e-14,
                        help="absolute tolerance in mRy")
    args = parser.parse_args()

    reference_path = args.case / "analytic_reference.json"
    reference = json.loads(reference_path.read_text())
    if reference.get("reference_model") != "finite open-boundary point-dipole sum":
        raise SystemExit(f"{args.case}: no atomistic analytic field reference")
    expected = [tuple(row) for row in reference["fields_t"]]
    simid = input_value(args.case / "inpsd.dat", "simid")
    field_path = args.case / f"befftot.{simid}.out"
    blocks = read_field_blocks(field_path, len(expected))
    if not blocks:
        raise SystemExit(f"{field_path}: no complete field block")
    # Initial field is the cleanest deterministic observable.  A later block
    # may contain evolved moments after the one integration step.
    key = min(blocks)
    observed = blocks[key]
    max_abs = max_rel = 0.0
    for site, (want, got) in enumerate(zip(expected, observed), start=1):
        for component, (expected_value, actual_value) in enumerate(zip(want, got)):
            absolute = abs(expected_value - actual_value)
            relative = absolute / max(abs(expected_value), args.atol)
            max_abs, max_rel = max(max_abs, absolute), max(max_rel, relative)
            if absolute > args.atol + args.rtol * abs(expected_value):
                raise AssertionError(
                    f"site {site}, component {component}: expected {expected_value:.17e} T, "
                    f"got {actual_value:.17e} T (abs {absolute:.3e}, rel {relative:.3e})")
    if not all(math.isfinite(value) for field in observed for value in field):
        raise AssertionError("field output contains NaN or Inf")
    energy_path = args.case / f"totenergy.{simid}.out"
    energy_iteration, actual_energy = read_dipole_energy(energy_path)
    expected_energy = float(reference["dipole_energy_per_atom_mry"])
    energy_abs = abs(expected_energy - actual_energy)
    energy_rel = energy_abs / max(abs(expected_energy), args.energy_atol)
    if energy_abs > args.energy_atol + args.energy_rtol * abs(expected_energy):
        raise AssertionError(
            f"Dip energy: expected {expected_energy:.17e} mRy/atom, got {actual_energy:.17e} "
            f"at iteration {energy_iteration} (abs {energy_abs:.3e}, rel {energy_rel:.3e})")
    print(f"PASS {args.case.name}: iteration={key[0]} replica={key[1]} "
          f"field_max_abs={max_abs:.3e} T field_max_rel={max_rel:.3e}; "
          f"Dip={actual_energy:.8e} mRy/atom energy_rel={energy_rel:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
