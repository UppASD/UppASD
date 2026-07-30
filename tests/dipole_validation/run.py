#!/usr/bin/env python3
"""Run all finite analytic dipole validations against a CPU UppASD executable."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
ANALYTIC_CASES = (
    "sc3d_single_x",
    "sc3d_pair_axial_x",
    "sc3d_pair_transverse_z",
    "sc2d_finite_slab_x",
    "sc2d_finite_slab_z",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True, help="CPU UppASD executable")
    parser.add_argument("--workdir", type=Path, help="retain generated cases here")
    parser.add_argument("--keep", action="store_true", help="retain a temporary work directory")
    args = parser.parse_args()
    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"missing executable: {binary}")

    # This has no UppASD dependency: it verifies the independent periodic
    # Ewald oracle used by the forthcoming macrocell-PME tests before CPU
    # output is considered at all.
    subprocess.run([sys.executable, str(HERE / "test_periodic_ewald_reference.py")], check=True)

    temporary: tempfile.TemporaryDirectory[str] | None = None
    if args.workdir:
        workdir = args.workdir.resolve()
        if workdir.exists():
            raise SystemExit(f"workdir exists: {workdir}")
    else:
        temporary = tempfile.TemporaryDirectory(prefix="uppasd-dipole-validation-")
        workdir = Path(temporary.name) / "cases"
    subprocess.run([sys.executable, str(HERE / "generate_cases.py"), "--output", str(workdir)], check=True)

    for name in ANALYTIC_CASES:
        case = workdir / name
        print(f"RUN  {name}")
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = "1"
        completed = subprocess.run([str(binary)], cwd=case, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (case / "run.log").write_text(completed.stdout)
        if completed.returncode:
            raise SystemExit(f"FAIL {name}: UppASD exited {completed.returncode}; see {case / 'run.log'}")
        subprocess.run([sys.executable, str(HERE / "check_fields.py"), str(case)], check=True)

    print(f"PASS periodic Ewald self-check + {len(ANALYTIC_CASES)} analytic dipole cases")
    if temporary and args.keep:
        destination = Path.cwd() / workdir.name
        shutil.copytree(workdir, destination)
        print(f"retained {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
