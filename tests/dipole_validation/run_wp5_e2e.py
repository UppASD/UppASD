#!/usr/bin/env python3
"""Run narrow WP5 GPU-SD acceptance cases against an independent oracle.

The golden data is versioned rather than generated during CTest.  This runner
only launches the production executable and reads its Tesla/mRy writers.
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
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).with_name("wp5_e2e")
GOLDENS = json.loads(Path(__file__).with_name("wp5_production_goldens_v1.json").read_text())["cases"]
MU_B = 9.274009994e-24
MRY = 2.179872325e-21
PREFAC = 1.0e-7 * MU_B  # alat=1 in wp5_e2e/inpsd.dat
RELATIVE_TOLERANCE = 2.0e-8


def numeric_rows(path: Path) -> list[list[float]]:
    rows: list[list[float]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        rows.append([float(token.replace("D", "E").replace("d", "e")) for token in line.split()])
    return rows


def require_close(label: str, actual: float, expected: float) -> float:
    # Output writers use eight significant figures.  The absolute floor covers
    # the deliberately tiny Tesla/mRy values in these fixtures.
    tolerance = max(RELATIVE_TOLERANCE * abs(expected), 1.0e-42)
    error = abs(actual - expected)
    if not math.isclose(actual, expected, rel_tol=RELATIVE_TOLERANCE, abs_tol=tolerance):
        raise AssertionError(f"{label}: got {actual:.17g}, expected {expected:.17g}, tolerance {tolerance:.3g}")
    return error


def set_record(text: str, key: str, value: str) -> str:
    expression = re.compile(rf"(?mi)^(\s*{re.escape(key)})(?:\s+.*)?$")
    replacement = f"{key} {value}"
    text, count = expression.subn(replacement, text, count=1)
    return text if count else text + f"\n{replacement}\n"


def configure(run_dir: Path, updates: dict[str, str], jfile: str | None = None) -> None:
    input_file = run_dir / "inpsd.dat"
    text = input_file.read_text()
    for key, value in updates.items():
        # mcnstep must be parsed before Nstep in an MC input.  It is absent
        # from the SD fixture on purpose, so insert it at that point for the
        # MC rejection case rather than appending it after Nstep.
        if key == "mcnstep" and not re.search(r"(?mi)^\s*mcnstep\b", text):
            text, count = re.subn(r"(?mi)^(\s*nstep\b.*)$", f"mcnstep {value}\n\\1", text, count=1)
            if count:
                continue
        text = set_record(text, key, value)
    input_file.write_text(text)
    if jfile is not None:
        (run_dir / "jfile").write_text(jfile)


def launch(binary: Path, parent: Path, name: str, updates: dict[str, str], *,
           jfile: str | None = None, environment: dict[str, str] | None = None,
           expect_success: bool = True, posfile: str | None = None,
           momfile: str | None = None) -> tuple[Path, str]:
    run_dir = parent / name
    shutil.copytree(FIXTURE, run_dir)
    configure(run_dir, updates, jfile)
    if posfile is not None:
        (run_dir / "posfile").write_text(posfile)
    if momfile is not None:
        (run_dir / "momfile").write_text(momfile)
    env = None if environment is None else {**os.environ, **environment}
    completed = subprocess.run([str(binary)], cwd=run_dir, text=True, env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (run_dir / "run.log").write_text(completed.stdout)
    if expect_success != (completed.returncode == 0):
        state = "succeeded" if completed.returncode == 0 else f"failed ({completed.returncode})"
        raise RuntimeError(f"WP5 case {name} unexpectedly {state}; log follows:\n{completed.stdout}")
    return run_dir, completed.stdout


def final_fields(run_dir: Path, simid: str, count: int) -> list[tuple[float, float, float]]:
    rows = numeric_rows(run_dir / f"befftot.{simid}.out")
    result = [tuple(row[3:6]) for row in rows if int(round(row[0])) == 1]
    if len(result) != count:
        raise AssertionError(f"{simid}: expected {count} final field rows, found {len(result)}")
    return result


def final_energy(run_dir: Path, simid: str) -> list[float]:
    rows = numeric_rows(run_dir / f"totenergy.{simid}.out")
    if not rows:
        raise AssertionError(f"{simid}: no energy rows were written")
    # SD measures the state at the start of a time step.  Depending on the
    # configured buffer cadence, its final flush may therefore leave iteration
    # zero as the latest physical energy record for a one-step fixture.  Check
    # the latest writer record, rather than inventing a second measurement.
    return rows[-1]


def expected_energy_per_atom(case: dict) -> float:
    return case["energy_total"] / math.prod(case["grid"]) * PREFAC * MU_B / MRY


def check_dipole_case(run_dir: Path, simid: str, case: dict, ensembles: int = 1) -> tuple[float, float]:
    expected_fields = [tuple(PREFAC * value for value in field) for field in case["fields"]]
    measured_fields = final_fields(run_dir, simid, len(expected_fields) * ensembles)
    maximum_field = 0.0
    for ensemble in range(ensembles):
        for site, expected in enumerate(expected_fields):
            for component, (actual, wanted) in enumerate(zip(measured_fields[ensemble * len(expected_fields) + site], expected)):
                maximum_field = max(maximum_field, require_close(
                    f"{simid}: B[e={ensemble},site={site},component={component}] (T)", actual, wanted))
    energy = final_energy(run_dir, simid)
    expected = expected_energy_per_atom(case)
    try:
        maximum_energy = max(require_close(f"{simid}: total energy (mRy/atom)", energy[1], expected),
                             require_close(f"{simid}: Dip energy (mRy/atom)", energy[8], expected),
                             require_close(f"{simid}: total-minus-Dip (mRy/atom)", energy[1] - energy[8], 0.0))
    except AssertionError as error:
        log = (run_dir / "run.log").read_text()
        trace = "\n".join(line for line in log.splitlines() if "GPU_TEST_" in line)
        raise AssertionError(f"{error}; complete final output row={energy}; raw dipole trace:\n{trace}") from error
    return maximum_field, maximum_energy


def check_rejection(binary: Path, parent: Path, name: str, updates: dict[str, str], message: str,
                    *, jfile: str | None = None, posfile: str | None = None, momfile: str | None = None,
                    environment: dict[str, str] | None = None) -> None:
    run_dir = parent / name
    shutil.copytree(FIXTURE, run_dir)
    configure(run_dir, updates, jfile)
    if posfile is not None:
        (run_dir / "posfile").write_text(posfile)
    if momfile is not None:
        (run_dir / "momfile").write_text(momfile)
    env = None if environment is None else {**os.environ, **environment}
    completed = subprocess.run([str(binary)], cwd=run_dir, text=True, env=env,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (run_dir / "run.log").write_text(completed.stdout)
    if completed.returncode == 0 or message not in completed.stdout:
        raise AssertionError(f"{name}: expected rejection containing {message!r}; log follows:\n{completed.stdout}")


def main() -> None:
    global RELATIVE_TOLERANCE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--fp32", action="store_true", help="apply the WP9 fp32 physical-output budget")
    args = parser.parse_args()
    if args.fp32:
        RELATIVE_TOLERANCE = 5.0e-5
    binary = args.binary.resolve()
    if not binary.is_file():
        raise FileNotFoundError(binary)

    with tempfile.TemporaryDirectory(prefix="uppasd-wp5-e2e-") as temporary:
        parent = Path(temporary)
        maximum_field = maximum_energy = maximum_exchange_delta = 0.0

        # SDEalgh=5 executes the predictor/corrector path.  The M=1 fixture
        # checks the production Tesla scatter and measurement energy boundary.
        run_dir, _ = launch(binary, parent, "dipole_111", {"simid": "wp5_111"})
        field, energy = check_dipole_case(run_dir, "wp5_111", GOLDENS["one_cell_uniform"])
        maximum_field, maximum_energy = max(maximum_field, field), max(maximum_energy, energy)

        # A real M=4 GPU-SD run, with four identical physical replicas.  The
        # standalone GPU acceptance test separately verifies nonidentical
        # ensemble isolation; this catches SD initialization/measurement shape.
        run_dir, _ = launch(binary, parent, "dipole_m4", {"simid": "wp5_m4", "Mensemble": "4"})
        field, energy = check_dipole_case(run_dir, "wp5_m4", GOLDENS["one_cell_uniform"], ensembles=4)
        maximum_field, maximum_energy = max(maximum_field, field), max(maximum_energy, energy)

        # N=3 is deliberately not warp-aligned.  Its output is compared to a
        # separately frozen 3-cell independent Ewald calculation.
        n3_updates = {"simid": "wp5_n3", "ncell": "3 1 1"}
        run_dir, _ = launch(binary, parent, "dipole_n3", n3_updates)
        field, energy = check_dipole_case(run_dir, "wp5_n3", GOLDENS["three_cell_uniform"])
        maximum_field, maximum_energy = max(maximum_field, field), max(maximum_energy, energy)

        # WP8 integration seam: with a state uniform inside each 2x1x1
        # block, the projected Hamiltonian must reproduce the block-one
        # run exactly.  This exercises the Fortran map, device reduction,
        # coarse kernel construction, scatter to both member atoms, and the
        # per-atom measurement normalization in the real SD executable.
        fine_dir, _ = launch(binary, parent, "wp8_fine", {"simid": "wp8fin", "ncell": "2 1 1"})
        coarse_dir, _ = launch(binary, parent, "wp8_coarse", {
            "simid": "wp8crs", "ncell": "2 1 1", "block_size_x": "2"})
        fine_fields = final_fields(fine_dir, "wp8fin", 2)
        coarse_fields = final_fields(coarse_dir, "wp8crs", 2)
        for site, (fine_field, coarse_field) in enumerate(zip(fine_fields, coarse_fields)):
            for component, (fine_value, coarse_value) in enumerate(zip(fine_field, coarse_field)):
                maximum_field = max(maximum_field, require_close(
                    f"WP8 uniform block field [site={site},component={component}] (T)", coarse_value, fine_value))
        fine_energy, coarse_energy = final_energy(fine_dir, "wp8fin"), final_energy(coarse_dir, "wp8crs")
        maximum_energy = max(maximum_energy,
            require_close("WP8 uniform block total energy (mRy/atom)", coarse_energy[1], fine_energy[1]),
            require_close("WP8 uniform block Dip energy (mRy/atom)", coarse_energy[8], fine_energy[8]))

        # The same restricted-state identity must retain basis channels.  The
        # basis-resolved map has two atoms per primitive cell and two cells per
        # coarse block, so each coarse channel contains exactly two atoms.
        na2_posfile = "1 1 0.0 0.0 0.0\n2 1 0.5 0.0 0.0\n"
        na2_momfile = "1 1 1.0 1.0 -0.2 0.4\n2 1 1.0 -0.3 0.5 0.1\n"
        na2_fine, _ = launch(binary, parent, "wp8_na2_fine", {"simid": "wp8n2f", "ncell": "2 1 1"},
                              posfile=na2_posfile, momfile=na2_momfile)
        na2_coarse, _ = launch(binary, parent, "wp8_na2_coarse", {
            "simid": "wp8n2c", "ncell": "2 1 1", "block_size_x": "2"},
            posfile=na2_posfile, momfile=na2_momfile)
        for site, (fine_field, coarse_field) in enumerate(zip(final_fields(na2_fine, "wp8n2f", 4),
                                                               final_fields(na2_coarse, "wp8n2c", 4))):
            for component, (fine_value, coarse_value) in enumerate(zip(fine_field, coarse_field)):
                maximum_field = max(maximum_field, require_close(
                    f"WP8 NA=2 uniform block field [site={site},component={component}] (T)", coarse_value, fine_value))
        na2_fine_energy, na2_coarse_energy = final_energy(na2_fine, "wp8n2f"), final_energy(na2_coarse, "wp8n2c")
        maximum_energy = max(maximum_energy,
            require_close("WP8 NA=2 uniform block total energy (mRy/atom)", na2_coarse_energy[1], na2_fine_energy[1]),
            require_close("WP8 NA=2 uniform block Dip energy (mRy/atom)", na2_coarse_energy[8], na2_fine_energy[8]))

        # The exchange-plus-dipole total is checked by an OFF/on pair with an
        # identical exchange Hamiltonian.  Their field and total-energy delta
        # must be the independent dipole oracle; the exchange column itself is
        # invariant.  This proves additive scatter and total-energy accounting.
        # The fixture's primitive lattice vector is 3 Å, so the exchange
        # neighbour must be one primitive translation away, not one Å.
        # Keep exchange active but below the text writer's resolution relative
        # to the deliberately tiny Tesla-scale dipole fixture.  A conventional
        # J=1 produces ~10^2 T, which rounds the additive ~10^-31 T dipole
        # delta away in the 8-significant-digit production output.
        exchange_jfile = "1 1 3.0 0.0 0.0 1.0e-40\n"
        off_dir, off_log = launch(binary, parent, "exchange_off", {
            "simid": "wp5_xof", "gpu_dipole_mode": "OFF"}, jfile=exchange_jfile)
        on_dir, on_log = launch(binary, parent, "exchange_on", {
            "simid": "wp5_xon"}, jfile=exchange_jfile)
        off_fields = final_fields(off_dir, "wp5_xof", 1)
        on_fields = final_fields(on_dir, "wp5_xon", 1)
        expected_fields = GOLDENS["one_cell_uniform"]["fields"]
        try:
            for site, expected in enumerate(expected_fields):
                for component, expected_component in enumerate(expected):
                    maximum_exchange_delta = max(maximum_exchange_delta, require_close(
                        f"exchange-plus-dipole field delta [{site},{component}] (T)",
                        on_fields[site][component] - off_fields[site][component], PREFAC * expected_component))
        except AssertionError as error:
            on_dipole_log = "\n".join(line for line in on_log.splitlines()
                                        if "EWALD3D_FFT" in line or "GPU_TEST_" in line)
            raise AssertionError(f"{error}; OFF fields={off_fields}; ON fields={on_fields}; "
                                 f"ON dipole log:\n{on_dipole_log}") from error
        off_energy, on_energy = final_energy(off_dir, "wp5_xof"), final_energy(on_dir, "wp5_xon")
        expected = expected_energy_per_atom(GOLDENS["one_cell_uniform"])
        maximum_exchange_delta = max(maximum_exchange_delta,
            require_close("exchange-plus-dipole total delta (mRy/atom)", on_energy[1] - off_energy[1], expected),
            require_close("exchange-plus-dipole Dip (mRy/atom)", on_energy[8], expected),
            require_close("exchange energy unchanged (mRy/atom)", on_energy[2] - off_energy[2], 0.0))

        # Input gates execute before device allocations.  Regular NA>1 and
        # divisible coarse blocks are covered by their dedicated WP6/WP8
        # operator suites; this production harness retains only unsupported
        # execution modes and the required partial-edge rejection.
        check_rejection(binary, parent, "reject_mc", {"mode": "M", "mcnstep": "1"}, "GPU spin dynamics only")
        check_rejection(binary, parent, "reject_partial_coarse", {"ncell": "3 1 1", "block_size_x": "2"},
                        "requires positive macrocell blocks that divide N1, N2, and N3")
        check_rejection(binary, parent, "reject_open", {"BC": "P P 0"}, "requires periodic boundary conditions P P P")
        check_rejection(binary, parent, "reject_slab", {"BC": "P P F"}, "requires periodic boundary conditions P P P")
        check_rejection(binary, parent, "reject_surface", {"gpu_dipole_surface": "VACUUM"}, "requires gpu_dipole_surface TINFOIL")
        check_rejection(binary, parent, "reject_legacy", {"do_dip": "1"}, "requires do_dip=0")

        # Deterministic test seam for the real preflight: force one byte free
        # and require its diagnostic before normal device allocations proceed.
        check_rejection(binary, parent, "reject_preflight", {}, "GPU MEMORY BUDGET",
                        environment={"UPPASD_GPU_TEST_FREE_BYTES": "1", "UPPASD_GPU_SKIP_BUDGET": ""})

        print("wp5-e2e max_field_tesla_error={:.17g} max_energy_mry_error={:.17g} "
              "max_exchange_delta_error={:.17g}".format(maximum_field, maximum_energy, maximum_exchange_delta))


if __name__ == "__main__":
    main()
