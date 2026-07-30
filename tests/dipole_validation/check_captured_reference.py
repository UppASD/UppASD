#!/usr/bin/env python3
"""Compare a completed dipole run with a CPU reference capture."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from capture_references import atom_count
from check_fields import input_value, read_dipole_energy, read_field_blocks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--reference", type=Path, help="default: cpu_reference.json in case")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()
    case = args.case.resolve()
    reference = json.loads((args.reference or case / "cpu_reference.json").read_text())
    if reference.get("reference_model") != "captured CPU UppASD dipole output":
        raise SystemExit("reference is not a CPU dipole capture")
    simid = input_value(case / "inpsd.dat", "simid")
    fields = read_field_blocks(case / f"befftot.{simid}.out", sites=atom_count(case))
    expected = reference["field_blocks"]
    max_abs = max_rel = 0.0
    for key, wanted in expected.items():
        step, replica = map(int, key.split(":"))
        if (step, replica) not in fields:
            raise AssertionError(f"missing complete field block {key}")
        for want_site, got_site in zip(wanted, fields[(step, replica)]):
            for want, got in zip(want_site, got_site):
                absolute = abs(want - got)
                relative = absolute / max(abs(want), args.atol)
                max_abs, max_rel = max(max_abs, absolute), max(max_rel, relative)
                if absolute > args.atol + args.rtol * abs(want):
                    raise AssertionError(f"field differs: expected {want:.17e}, got {got:.17e}")
    _, energy = read_dipole_energy(case / f"totenergy.{simid}.out")
    wanted_energy = float(reference["dipole_energy_per_atom_mry"])
    if abs(energy - wanted_energy) > args.atol + args.rtol * abs(wanted_energy):
        raise AssertionError(f"Dip differs: expected {wanted_energy:.17e}, got {energy:.17e}")
    print(f"PASS captured CPU reference: field_abs={max_abs:.3e}, field_rel={max_rel:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
