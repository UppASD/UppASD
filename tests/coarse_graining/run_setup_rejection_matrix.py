#!/usr/bin/env python3
"""Require input-reachable adaptive-CG capability errors at setup time."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from fixture_dependencies import STATIC_MIXED_CASE


def replace_line(text: str, keyword: str, value: str) -> str:
    pattern = re.compile(rf"(?im)^{re.escape(keyword)}\s+.*$")
    replacement = f"{keyword} {value}"
    if not pattern.search(text):
        raise AssertionError(f"base fixture has no {keyword} line")
    return pattern.sub(replacement, text, count=1)


def remove_line(text: str, keyword: str) -> str:
    pattern = re.compile(rf"(?im)^{re.escape(keyword)}\s+.*\n?")
    if not pattern.search(text):
        raise AssertionError(f"base fixture has no {keyword} line")
    return pattern.sub("", text, count=1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).with_name("e2e")
    base = root / STATIC_MIXED_CASE
    source = (base / "inpsd.dat").read_text()

    cases = {
        "mode": (replace_line(source, "mode", "M"), "mode:"),
        "initial-phase": (replace_line(source, "ip_mode", "X"), "ip_mode:"),
        "restart": (replace_line(source, "initmag", "4"), "initmag=4"),
        "temperature": (replace_line(source, "temp", "10.0"), "Temp/do_qhb/do_3tm"),
        "stochastic-integrator": (replace_line(source, "SDEalgh", "2"), "SDEalgh/llg"),
        "boundary": (replace_line(source, "BC", "O P P"), "BC:"),
        "missing-characteristic-length": (remove_line(source, "alat"), "alat:"),
        "partial-block": (replace_line(source, "block_size_y", "3"), "block_size_x/y/z"),
        "legacy-dipole": (replace_line(source, "do_dip", "1") if re.search(r"(?im)^do_dip\s+", source) else source + "\ndo_dip 1\n", "do_dip:"),
        "tensor-exchange": (source + "\ndo_jtensor 1\n", "Hamiltonian capability"),
        "symmetric-anisotropic-exchange": (source + "\nsa ../jfile\n", "Hamiltonian capability"),
        "pseudo-dipolar": (source + "\npd ../jfile\n", "Hamiltonian capability"),
        "biquadratic": (source + "\nbq ../jfile\n", "Hamiltonian capability"),
        "ring-exchange": (source + "\nring ../jfile\n", "Hamiltonian capability"),
        "chiral": (source + "\nchir ../jfile\n", "Hamiltonian capability"),
        "four-site": (source + "\nfourx ../jfile\n", "Hamiltonian capability"),
        "biquadratic-dmi": (source + "\nbiqdm ../jfile\n", "Hamiltonian capability"),
        "random-anisotropy": (source + "\nrandom_anisotropy Y\n", "random_anisotropy"),
        "local-spin-fluctuation": (source + "\ndo_lsf Y\n", "do_ralloy/Nchmax/do_lsf/ind_mom_flag"),
        "induced-moment": (source + "\nind_mom_flag Y\n", "do_ralloy/Nchmax/do_lsf/ind_mom_flag"),
        "sparse": (source + "\ndo_sparse Y\n", "do_sparse/do_reduced"),
        "reduced": (source + "\ndo_reduced Y\n", "do_sparse/do_reduced"),
        "fixed-moment": (source + "\ndo_fixed_mom Y\n", "do_sparse/do_reduced"),
        "legacy-sampler": (source + "\ndo_autocorr Y\n", "do_autocorr/do_spintemp/do_kmc"),
        "spin-temperature": (source + "\ndo_spintemp Y\n", "do_autocorr/do_spintemp/do_kmc"),
        "kmc": (source + "\ndo_kmc Y\n", "do_autocorr/do_spintemp/do_kmc"),
        "demag": (source + "\ndemag Y\n", "hfield/do_bpulse/demag"),
        "pulse": (source + "\ndo_bpulse 1\n", "hfield/do_bpulse/demag"),
        "external-field": (source + "\nhfield 1.0 0.0 0.0\n", "hfield/do_bpulse/demag"),
        "legacy-energy-output": (replace_line(source, "plotenergy", "1"), "plotenergy:"),
    }

    binary = args.binary.resolve()
    failures: list[str] = []
    with tempfile.TemporaryDirectory(prefix="cg13-rejections-") as temporary:
        temporary_root = Path(temporary)
        for name, (input_text, diagnostic) in cases.items():
            case = temporary_root / name
            shutil.copytree(base, case)
            (case / "inpsd.dat").write_text(input_text)
            result = subprocess.run(
                [str(binary)], cwd=case, text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, timeout=60, check=False,
            )
            setup_error = any(marker in result.stdout for marker in (
                "AdaptiveCG setup rejected:", "Input data is inconsistent",
                "**** ERROR ****", "Fortran runtime error", "ERROR STOP",
            ))
            if result.returncode == 0 or not setup_error or "AdaptiveCG: capability accepted" in result.stdout:
                failures.append(f"{name}: expected setup-time rejection ({diagnostic})\n{result.stdout}")

    if failures:
        raise AssertionError("\n".join(failures))
    print(f"CG-13 setup-rejection matrix passed ({len(cases)} cases)")


if __name__ == "__main__":
    main()
