#!/usr/bin/env python3
"""Production OPEN_FFT coarse-block equivalence and rejection checks.

Expected values are relations between a block-one finite run and its exact
uniform-block projection; this script never rewrites Luna's oracle goldens.
"""
from __future__ import annotations

import argparse
import tempfile
from pathlib import Path

import run_wp5_e2e as wp5


OPEN = {"gpu_dipole_mode": "OPEN_FFT", "BC": "0 0 0"}


def projected_fields(fine: list[tuple[float, float, float]], basis: int,
                     ensembles: int) -> list[tuple[float, float, float]]:
    """Average the two x-neighbour fine cells for each basis channel."""
    atoms_per_ensemble = 2 * basis
    result = []
    for ensemble in range(ensembles):
        fields = fine[ensemble * atoms_per_ensemble:(ensemble + 1) * atoms_per_ensemble]
        for cell in range(2):
            for channel in range(basis):
                first = fields[channel]
                second = fields[basis + channel]
                # Coarse scatter writes the conjugate macro field to both
                # member atoms; the fine reference is averaged back to that
                # same macro coordinate.
                result.append(tuple((left + right) / 2.0 for left, right in zip(first, second)))
    return result


def check_projection(binary: Path, parent: Path, name: str, *, basis: int,
                     ensembles: int, posfile: str | None = None, momfile: str | None = None,
                     jfile: str | None = None) -> tuple[float, float]:
    if len(name) > 7:
        raise ValueError("OPEN_FFT E2E simid base must leave room for its fine/coarse suffix")
    updates = {**OPEN, "simid": f"{name}f", "ncell": "2 1 1", "Mensemble": str(ensembles)}
    fine_dir, _ = wp5.launch(binary, parent, f"{name}_fine", updates, posfile=posfile, momfile=momfile,
                             jfile=jfile)
    coarse_updates = {**updates, "simid": f"{name}c", "block_size_x": "2"}
    coarse_dir, coarse_log = wp5.launch(binary, parent, f"{name}_coarse", coarse_updates,
                                        posfile=posfile, momfile=momfile, jfile=jfile)
    if "OPEN_FFT device allocation:" not in coarse_log or "block 2 x 1 x 1" not in coarse_log:
        raise AssertionError("OPEN_FFT coarse run did not report its allocation and enabled block shape")
    fine = wp5.final_fields(fine_dir, f"{name}f", 2 * basis * ensembles)
    coarse = wp5.final_fields(coarse_dir, f"{name}c", 2 * basis * ensembles)
    expected = projected_fields(fine, basis, ensembles)
    maximum_field = 0.0
    if jfile is None:
        for index, (actual, wanted) in enumerate(zip(coarse, expected)):
            for component, (value, reference) in enumerate(zip(actual, wanted)):
                maximum_field = max(maximum_field, wp5.require_close(
                    f"{name} projected field [{index},{component}]", value, reference))
    # Exchange remains atomistic and is deliberately not block-projected, so
    # its local boundary field cannot equal a fine-block average.  The total
    # and dipole energy comparison below is the valid coexistence assertion.
    fine_energy = wp5.final_energy(fine_dir, f"{name}f")
    coarse_energy = wp5.final_energy(coarse_dir, f"{name}c")
    maximum_energy = max(
        wp5.require_close(f"{name} projected total energy", coarse_energy[1], fine_energy[1]),
        wp5.require_close(f"{name} projected dip energy", coarse_energy[8], fine_energy[8]),
    )
    return maximum_field, maximum_energy


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)
    with tempfile.TemporaryDirectory(prefix="uppasd-wp10-open-coarse-") as temporary:
        parent = Path(temporary)
        maximum_field, maximum_energy = check_projection(binary, parent, "o11", basis=1, ensembles=1)
        na2_posfile = "1 1 0.0 0.0 0.0\n2 1 0.23 0.17 0.11\n"
        na2_momfile = "1 1 1.0 0.4 -0.2 0.7\n2 1 1.0 -0.3 0.5 0.1\n"
        field, energy = check_projection(binary, parent, "o24", basis=2, ensembles=4,
                                         posfile=na2_posfile, momfile=na2_momfile)
        maximum_field, maximum_energy = max(maximum_field, field), max(maximum_energy, energy)

        # Exchange remains additive: repeat the exact projected relation with
        # a nonzero exchange Hamiltonian and do not use it to form an oracle.
        exchange = "1 1 3.0 0.0 0.0 1.0e-40\n"
        field, energy = check_projection(binary, parent, "oxe", basis=1, ensembles=1,
                                         jfile=exchange)
        maximum_field, maximum_energy = max(maximum_field, field), max(maximum_energy, energy)
        wp5.check_rejection(binary, parent, "open_reject_partial",
                            {**OPEN, "ncell": "3 1 1", "block_size_x": "2"},
                            "OPEN_FFT requires positive macrocell blocks that divide N1, N2, and N3")
        print("wp10-open-coarse-e2e max_field_tesla_error={:.17g} max_energy_mry_error={:.17g}".format(
            maximum_field, maximum_energy))


if __name__ == "__main__":
    main()
