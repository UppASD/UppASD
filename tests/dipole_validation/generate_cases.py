#!/usr/bin/env python3
"""Generate deterministic dipole-validation cases from the SC base fixtures.

The references written beside each generated input are independent point-dipole
calculations in SI units.  They are deliberately *not* CPU benchmark output.
They cover the finite/open-boundary cases needed to validate the tensor sign,
prefactor, self term and energy before a GPU FFT backend is enabled.  Periodic
3D and 2D-slab PME checks need their own reciprocal-space analytic modes and
are intentionally not faked by this finite-sum generator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
MU_B = 9.2740100783e-24  # J/T, CODATA 2018 exact value used for test metadata
MU0_OVER_4PI = 1.0e-7   # T m / A
MRY_J = 2.179872325e-21 # UppASD Constants::mry


def set_record(text: str, key: str, value: str) -> str:
    """Replace an UppASD input record, appending it when absent."""
    lines = text.splitlines()
    for index, line in enumerate(lines):
        words = line.split()
        if words and words[0].lower() == key.lower():
            lines[index] = f"{key} {value}"
            return "\n".join(lines) + "\n"
    return text.rstrip() + f"\n{key} {value}\n"


def points(shape: tuple[int, int, int], alat_m: float) -> list[tuple[float, float, float]]:
    # The one-site SC fixture uses x as the fastest lattice index.
    return [
        (ix * alat_m, iy * alat_m, iz * alat_m)
        for iz in range(shape[2])
        for iy in range(shape[1])
        for ix in range(shape[0])
    ]


def dipole_reference(shape: tuple[int, int, int], moment_mu_b: tuple[float, float, float],
                     alat_m: float) -> dict:
    """Finite point-dipole field and energy for one unit moment at every site."""
    xyz = points(shape, alat_m)
    fields: list[tuple[float, float, float]] = []
    for target in xyz:
        bx = by = bz = 0.0
        for source in xyz:
            rx, ry, rz = target[0] - source[0], target[1] - source[1], target[2] - source[2]
            r2 = rx * rx + ry * ry + rz * rz
            if r2 == 0.0:
                continue
            r = r2 ** 0.5
            mdotr = moment_mu_b[0] * rx + moment_mu_b[1] * ry + moment_mu_b[2] * rz
            pref = MU0_OVER_4PI * MU_B
            bx += pref * (3.0 * rx * mdotr / r**5 - moment_mu_b[0] / r**3)
            by += pref * (3.0 * ry * mdotr / r**5 - moment_mu_b[1] / r**3)
            bz += pref * (3.0 * rz * mdotr / r**5 - moment_mu_b[2] / r**3)
        fields.append((bx, by, bz))
    energy_j = -0.5 * MU_B * sum(
        moment_mu_b[0] * field[0] + moment_mu_b[1] * field[1] + moment_mu_b[2] * field[2]
        for field in fields
    )
    return {
        "reference_model": "finite open-boundary point-dipole sum",
        "field_units": "T",
        "energy_units": "J",
        "mu0_over_4pi": MU0_OVER_4PI,
        "mu_b_j_per_t": MU_B,
        "alat_m": alat_m,
        "moment_mu_b": list(moment_mu_b),
        "shape": list(shape),
        "fields_t": [list(field) for field in fields],
        "dipole_energy_j": energy_j,
        # totenergy reports the ensemble average per atom in mRy.
        "dipole_energy_per_atom_mry": energy_j / len(xyz) / MRY_J,
    }


def write_case(name: str, base: Path, output: Path, shape: tuple[int, int, int],
               moment: tuple[float, float, float], *, do_dip: int,
               bc: tuple[str, str, str], macro_block: tuple[int, int, int] | None = None,
               analytic: bool = True, simid: str, ensembles: int = 1) -> None:
    if len(simid) > 8:
        raise ValueError(f"simid must fit UppASD's 8-character field filename: {simid}")
    destination = output / name
    shutil.copytree(base, destination)
    inp = destination / "inpsd.dat"
    text = inp.read_text()
    records = {
        "simid": simid,
        "ncell": f"{shape[0]} {shape[1]} {shape[2]}",
        "BC": " ".join(bc),
        "do_dip": str(do_dip),
        "Nstep": "1",
        "Mensemble": str(ensembles),
        "temp": "0.0",
        "plotenergy": "1",
        "ene_step": "1",
        "do_avrg": "N",
        # prn_totalbfields writes befftot.<simid>.out in Tesla.  This is the
        # observable checked against analytic_reference.json, not CPU output.
        "do_prn_beff": "Y",
        "beff_step": "1",
        "beff_buff": "1",
    }
    if macro_block:
        records.update({
            "do_macro_cells": "Y",
            "block_size_x": str(macro_block[0]),
            "block_size_y": str(macro_block[1]),
            "block_size_z": str(macro_block[2]),
        })
    for key, value in records.items():
        text = set_record(text, key, value)
    inp.write_text(text)
    (destination / "momfile").write_text(
        f"1 1 1.0 {moment[0]:.17g} {moment[1]:.17g} {moment[2]:.17g}\n"
    )
    if analytic:
        reference = dipole_reference(shape, moment, 1.0e-9)
        (destination / "analytic_reference.json").write_text(json.dumps(reference, indent=2) + "\n")
    else:
        metadata = {
            "reference_model": "macrocell aggregation/layout only",
            "note": "This case is not a CPU benchmark or an atomistic analytic oracle. "
                    "Use it to validate macro-cell populations and edge-block handling.",
            "shape": list(shape),
            "macro_block": list(macro_block or ()),
            "ensembles": ensembles,
        }
        (destination / "analytic_reference.json").write_text(json.dumps(metadata, indent=2) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tests/dipole_validation/generated")
    parser.add_argument("--force", action="store_true", help="replace an existing output directory")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output.exists():
        if not args.force:
            raise SystemExit(f"{args.output} exists; pass --force to replace it")
        shutil.rmtree(args.output)
    args.output.mkdir(parents=True)
    sc3d = ROOT / "tests/SC3d_dipole"
    sc2d = ROOT / "tests/SC2d_dipole"
    for base in (sc3d, sc2d):
        if not (base / "inpsd.dat").is_file():
            raise SystemExit(f"missing base fixture: {base}")

    # All analytic cases explicitly use open boundaries. The brute-force CPU
    # dipole implementation is a finite sum, so leaving PBC in the input here
    # would label a finite reference as a periodic one.
    write_case("sc3d_single_x", sc3d, args.output, (1, 1, 1), (1.0, 0.0, 0.0),
               do_dip=1, bc=("0", "0", "0"), simid="d3single")
    write_case("sc3d_pair_axial_x", sc3d, args.output, (2, 1, 1), (1.0, 0.0, 0.0),
               do_dip=1, bc=("0", "0", "0"), simid="d3pairx")
    write_case("sc3d_pair_transverse_z", sc3d, args.output, (2, 1, 1), (0.0, 0.0, 1.0),
               do_dip=1, bc=("0", "0", "0"), simid="d3pairz")
    write_case("sc2d_finite_slab_z", sc2d, args.output, (4, 4, 1), (0.0, 0.0, 1.0),
               do_dip=1, bc=("0", "0", "0"), simid="d2slabz")
    write_case("sc2d_finite_slab_x", sc2d, args.output, (4, 4, 1), (1.0, 0.0, 0.0),
               do_dip=1, bc=("0", "0", "0"), simid="d2slabx")

    # Uniform macro grid: the first CPU/GPU capture baseline.  It uses four
    # full 2x2 cells, so the FFT kernel is translationally invariant.
    write_case("sc2d_macro_uniform", sc2d, args.output, (4, 4, 1), (0.0, 0.0, 1.0),
               do_dip=2, bc=("0", "0", "0"), macro_block=(2, 2, 1), analytic=False,
               simid="d2unifrm")
    write_case("sc2d_macro_uniform_m4", sc2d, args.output, (4, 4, 1), (0.0, 0.0, 1.0),
               do_dip=2, bc=("0", "0", "0"), macro_block=(2, 2, 1), analytic=False,
               simid="d2unifm4", ensembles=4)

    # Partial-edge layout test.  This stays generated even when a local CPU
    # build cannot yet run it: an edge capture is a CV6 gate, not a reason to
    # silently validate only divisible grids.
    write_case("sc2d_macro_edges", sc2d, args.output, (65, 63, 1), (0.0, 0.0, 1.0),
               do_dip=2, bc=("0", "0", "0"), macro_block=(8, 8, 1), analytic=False,
               simid="d2macro")
    print(f"generated dipole validation cases in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
