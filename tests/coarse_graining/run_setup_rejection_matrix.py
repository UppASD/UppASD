#!/usr/bin/env python3
"""Require input-reachable adaptive-CG capability errors at setup time."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from fixture_dependencies import CLUSTER_ANISOTROPY_REJECTION_CASE, STATIC_MIXED_CASE


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
        # Switching to mode M without mcnstep trips an unrelated legacy
        # input-consistency fallback in read_parameters (nstep is read
        # before mcnstep is known, so it silently seeds mcnstep from nstep
        # and marks sane_input=.false.) before AdaptiveCG's own mode
        # rejection is ever reached. Setting mcnstep in the same position as
        # "mode" keeps it ahead of the fixture's existing Nstep line, which
        # is required for the fallback not to fire.
        "mode": (replace_line(source, "mode", "M\nmcnstep 2"), "mode:"),
        "initial-phase": (replace_line(source, "ip_mode", "X"), "ip_mode:"),
        # initmag=4 (restart) used to be rejected here unconditionally; it
        # no longer is (RCG-04B follow-up -- see
        # docs/RCG-04_MOVING_E2E_EVIDENCE.md section 10.4/10.7 and
        # tests/coarse_graining/e2e/initmag_restart_atomistic, the positive
        # regression proving the capability now works end to end). There is
        # therefore no longer a setup-time rejection to assert for this
        # case; it has been removed rather than left asserting stale
        # behaviour.
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
        "canonical-sparse": (source + "\ncpu_ham_backend sparse\n", "cpu_ham_backend/do_sparse/do_reduced"),
        "canonical-convolution": (source + "\ncpu_ham_backend convolution\n", "cpu_ham_backend/do_sparse/do_reduced"),
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
        # The base fixture's inpsd.dat reads "posfile ../posfile",
        # "momfile ../momfile", and "exchange ../jfile" (and every
        # Hamiltonian-capability case appends its own "<kw> ../jfile" line),
        # all resolved relative to each case directory one level below
        # temporary_root. Without copies of these tracked siblings here,
        # posfile/momfile silently fall back to a single-basis-atom default
        # that misclassifies the geometry as an explicit-device rejection,
        # and jfile-dependent cases fail on a missing-file stop, in both
        # cases before the case can reach its intended setup-time rejection.
        for shared_fixture in ("posfile", "momfile", "jfile", "kfile_cg_x"):
            shutil.copy(root / shared_fixture, temporary_root / shared_fixture)
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

        # RCG-03 ANI-NONUNIFORM-REJECT: unlike the cases above, this fixture
        # is not a text mutation of static_mixed -- it is its own tracked
        # do_cluster input (posfile_clus/momfile_clus/exchange_clus/
        # anisotropy_clus), copied verbatim, that reaches
        # build_production_anisotropy's cell-periodicity check through the
        # ordinary do_cluster embedding path rather than fault injection.
        cluster_case = temporary_root / "anisotropy-cluster-nonuniform"
        shutil.copytree(root / CLUSTER_ANISOTROPY_REJECTION_CASE, cluster_case)
        result = subprocess.run(
            [str(binary)], cwd=cluster_case, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=60, check=False,
        )
        if result.returncode == 0 or "anisotropy: basis" not in result.stdout or \
                "AdaptiveCG: capability accepted" in result.stdout:
            failures.append(
                "anisotropy-cluster-nonuniform: expected setup-time rejection "
                f"(anisotropy: basis)\n{result.stdout}"
            )

    if failures:
        raise AssertionError("\n".join(failures))
    print(f"CG-13 setup-rejection matrix passed ({len(cases) + 1} cases)")


if __name__ == "__main__":
    main()
