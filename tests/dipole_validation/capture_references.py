#!/usr/bin/env python3
"""Capture deterministic CPU dipole fields/energies for non-analytic modes.

The result is a *CPU reference capture*, not an analytic oracle.  Use it for
the macrocell (`do_dip=2`) and CPU FFT (`do_dip=3`) regression baselines only;
keep `check_fields.py` as the independent point-dipole validation.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from check_fields import input_value, read_dipole_energy, read_field_blocks


def atom_count(case: Path) -> int:
    """Return Natom from the basis file and the three ncell dimensions."""
    basis = sum(1 for line in (case / "posfile").read_text().splitlines() if line.strip())
    for line in (case / "inpsd.dat").read_text().splitlines():
        words = line.split()
        if words and words[0].lower() == "ncell":
            if len(words) < 4:
                raise ValueError("ncell must contain three dimensions")
            n1, n2, n3 = (int(value) for value in words[1:4])
            break
    else:
        raise ValueError("missing ncell while determining atom count")
    if basis == 0 or min(n1, n2, n3) <= 0:
        raise ValueError("invalid posfile or ncell while determining atom count")
    return basis * n1 * n2 * n3


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="generated and runnable case directory")
    parser.add_argument("--binary", type=Path, required=True, help="CPU UppASD executable")
    parser.add_argument("--output", type=Path, help="capture path (default: cpu_reference.json in case)")
    args = parser.parse_args()

    case, binary = args.case.resolve(), args.binary.resolve()
    if not (case / "inpsd.dat").is_file() or not binary.is_file():
        raise SystemExit("case must contain inpsd.dat and --binary must name an executable")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    completed = subprocess.run([str(binary)], cwd=case, env=env, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (case / "cpu_reference_run.log").write_text(completed.stdout)
    if completed.returncode:
        raise SystemExit(f"UppASD exited {completed.returncode}; see cpu_reference_run.log")

    simid = input_value(case / "inpsd.dat", "simid")
    fields = read_field_blocks(case / f"befftot.{simid}.out", sites=atom_count(case))
    if not fields:
        raise SystemExit("no complete effective-field block was written")
    iteration, energy = read_dipole_energy(case / f"totenergy.{simid}.out")
    capture = {
        "reference_model": "captured CPU UppASD dipole output",
        "do_dip": int(input_value(case / "inpsd.dat", "do_dip")),
        "simid": simid,
        "field_units": "T",
        "dipole_energy_units": "mRy per atom",
        "dipole_energy_iteration": iteration,
        "dipole_energy_per_atom_mry": energy,
        "field_blocks": {
            f"{step}:{replica}": values for (step, replica), values in sorted(fields.items())
        },
    }
    output = (args.output or case / "cpu_reference.json").resolve()
    output.write_text(json.dumps(capture, indent=2) + "\n")
    print(f"captured do_dip={capture['do_dip']} reference to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
