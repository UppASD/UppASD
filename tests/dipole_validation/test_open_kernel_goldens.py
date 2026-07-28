#!/usr/bin/env python3
"""Compare the C++ host open-kernel fields with Luna's frozen WP10.1 cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent


def flatten_moments(record: dict) -> list[float]:
    # C++ active index: component + 3*(basis + NA*(cell + active_cells*ensemble)).
    return [
        component
        for ensemble in record["moments_mu_b"]
        for cell in ensemble
        for basis in cell
        for component in basis
    ]


def flatten_fields(record: dict) -> list[float]:
    return [
        component
        for ensemble in record["fields_dimensionless"]
        for cell in ensemble
        for basis in cell
        for component in basis
    ]


def run_case(driver: Path, record: dict, fft_grid: tuple[int, int, int]) -> tuple[list[float], list[float]]:
    tokens: list[str] = []
    tokens.extend(str(value) for value in record["grid"])
    tokens.extend(str(value) for value in fft_grid)
    tokens.extend(("1", "1", "1"))
    tokens.extend((str(len(record["basis_offsets"])), str(len(record["moments_mu_b"]))))
    tokens.extend(repr(value) for value in record["primitive_vectors"])
    tokens.extend(repr(value) for offset in record["basis_offsets"] for value in offset)
    tokens.extend(repr(value) for value in flatten_moments(record))
    completed = subprocess.run(
        [str(driver)], input=" ".join(tokens), text=True, capture_output=True, check=True
    )
    values = completed.stdout.split()
    if len(values) < 10:
        raise AssertionError(f"{record['name']}: malformed driver output: {completed.stdout!r}")
    field_count = int(values[0])
    diagnostics = [float(value) for value in values[1:10]]
    fields = [float(value) for value in values[10:]]
    if len(fields) != field_count:
        raise AssertionError(
            f"{record['name']}: driver reported {field_count} fields but returned {len(fields)}"
        )
    return fields, diagnostics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver", required=True, type=Path)
    args = parser.parse_args()
    goldens = json.loads((HERE / "open_fft_goldens.json").read_text())
    if goldens["oracle"] != "independent finite open-boundary point-dipole sum":
        raise AssertionError("unexpected WP10.1 golden provenance")

    maximum_absolute_error = 0.0
    maximum_scaled_error = 0.0
    maximum_energy_absolute_error = 0.0
    maximum_energy_scaled_error = 0.0
    maximum_reciprocity_error = 0.0
    minimum_positive_r2 = float("inf")
    comparisons = 0
    cases_run = 0
    for record in goldens["cases"]:
        grid = tuple(record["grid"])
        minimum_padding = tuple(2 * extent - 1 for extent in grid)
        paddings = [minimum_padding]
        # Recheck one non-cubic golden with a legal, deliberately oversized P.
        if record["name"] == "nonuniform_2x3x1":
            paddings.append(tuple(minimum_padding[axis] + (3, 2, 4)[axis] for axis in range(3)))
        expected = flatten_fields(record)
        moments = flatten_moments(record)
        values_per_ensemble = 3 * grid[0] * grid[1] * grid[2] * len(record["basis_offsets"])
        for padding in paddings:
            actual, diagnostics = run_case(args.driver, record, padding)
            if len(actual) != len(expected):
                raise AssertionError(f"{record['name']}: field count mismatch")
            active_cells, fft_cells, batches, nonfinite, all_finite, minimum_r2, reciprocity, self_abs, gap_abs = diagnostics
            expected_active = grid[0] * grid[1] * grid[2]
            expected_fft = padding[0] * padding[1] * padding[2]
            expected_batches = 9 * len(record["basis_offsets"]) ** 2
            if (int(active_cells), int(fft_cells), int(batches)) != (
                expected_active, expected_fft, expected_batches
            ):
                raise AssertionError(f"{record['name']}: count diagnostics are wrong")
            if (
                nonfinite != 0.0
                or all_finite != 1.0
                or reciprocity > 1.0e-11
                or self_abs != 0.0
                or gap_abs != 0.0
            ):
                raise AssertionError(f"{record['name']}: integrity diagnostics failed: {diagnostics}")
            if expected_active * len(record["basis_offsets"]) > 1 and not minimum_r2 > 0.0:
                raise AssertionError(f"{record['name']}: minimum nonself distance is not positive")
            maximum_reciprocity_error = max(maximum_reciprocity_error, reciprocity)
            if minimum_r2 > 0.0:
                minimum_positive_r2 = min(minimum_positive_r2, minimum_r2)
            for got, want in zip(actual, expected):
                error = abs(got - want)
                scaled = error / max(1.0, abs(want))
                maximum_absolute_error = max(maximum_absolute_error, error)
                maximum_scaled_error = max(maximum_scaled_error, scaled)
                comparisons += 1
                if error > 1.0e-12 * max(1.0, abs(want)):
                    raise AssertionError(
                        f"{record['name']} P={padding}: got {got:.17g}, "
                        f"want {want:.17g}, error {error:.3e}"
                    )
            for ensemble in range(len(record["moments_mu_b"])):
                start = ensemble * values_per_ensemble
                stop = start + values_per_ensemble
                actual_energy = -0.5 * sum(
                    moment * field for moment, field in zip(moments[start:stop], actual[start:stop])
                )
                expected_energy = -0.5 * sum(
                    moment * field for moment, field in zip(moments[start:stop], expected[start:stop])
                )
                energy_error = abs(actual_energy - expected_energy)
                energy_scaled = energy_error / max(1.0, abs(expected_energy))
                maximum_energy_absolute_error = max(maximum_energy_absolute_error, energy_error)
                maximum_energy_scaled_error = max(maximum_energy_scaled_error, energy_scaled)
                if energy_error > 1.0e-12 * max(1.0, abs(expected_energy)):
                    raise AssertionError(
                        f"{record['name']} P={padding} ensemble={ensemble}: "
                        f"energy got {actual_energy:.17g}, want {expected_energy:.17g}, "
                        f"error {energy_error:.3e}"
                    )
            cases_run += 1

    print(
        "open-host-goldens "
        f"golden_cases={len(goldens['cases'])} builds={cases_run} comparisons={comparisons} "
        f"absolute_max={maximum_absolute_error:.17g} scaled_max={maximum_scaled_error:.17g} "
        f"energy_absolute_max={maximum_energy_absolute_error:.17g} "
        f"energy_scaled_max={maximum_energy_scaled_error:.17g} "
        f"reciprocity_max={maximum_reciprocity_error:.17g} "
        f"minimum_nonself_r2={minimum_positive_r2:.17g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
