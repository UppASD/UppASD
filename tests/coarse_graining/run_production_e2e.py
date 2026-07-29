#!/usr/bin/env python3
"""Launch the normal UppASD executable for CG-10.5 production cases."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


def run_case(binary: Path, case: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary)],
        cwd=case,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
        check=False,
    )


def metric(output: str, name: str) -> int:
    match = re.search(rf"{re.escape(name)}=(\d+)", output)
    if not match:
        raise AssertionError(f"missing {name} in output\n{output}")
    return int(match.group(1))

def restart_state(case: Path) -> list[list[float]]:
    restart = next(case.glob("restart.*.out"))
    state: list[list[float]] = []
    for line in restart.read_text().splitlines():
        fields = line.split()
        if len(fields) == 7 and fields[0].lstrip("+-").isdigit():
            state.append([float(value) for value in fields[3:7]])
    return state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--gpu-binary", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = Path(__file__).with_name("e2e")

    off = run_case(binary, root / "feature_off")
    assert off.returncode == 0, off.stdout
    assert "AdaptiveCG:" not in off.stdout, off.stdout

    outputs: dict[str, str] = {}
    for name in ("static_all_fine", "static_all_coarse", "static_mixed", "adaptive_mixed"):
        result = run_case(binary, root / name)
        assert result.returncode == 0, f"{name}\n{result.stdout}"
        assert "AdaptiveCG: capability accepted" in result.stdout
        outputs[name] = result.stdout

    fine_updates = metric(outputs["static_all_fine"], "active_atom_updates")
    coarse_updates = metric(outputs["static_all_coarse"], "active_atom_updates")
    mixed_updates = metric(outputs["static_mixed"], "active_atom_updates")
    baseline = metric(outputs["static_all_coarse"], "baseline_atom_updates")
    reduced = metric(outputs["static_all_coarse"], "reduced_atom_updates")
    assert fine_updates == baseline
    assert coarse_updates < mixed_updates < fine_updates
    assert reduced > 0
    assert metric(outputs["adaptive_mixed"], "accepted_transitions") > 0
    off_state = restart_state(root / "feature_off")
    fine_state = restart_state(root / "static_all_fine")
    assert len(off_state) == len(fine_state) == 48
    assert max(
        abs(reference - adaptive)
        for reference_row, adaptive_row in zip(off_state, fine_state)
        for reference, adaptive in zip(reference_row, adaptive_row)
    ) <= 1.0e-10

    bad = run_case(binary, root / "invalid_partial_block")
    assert bad.returncode != 0
    assert "block_size_x/y/z" in bad.stdout
    bad_mask = run_case(binary, root / "invalid_mask")
    assert bad_mask.returncode != 0
    assert "duplicate block id" in bad_mask.stdout
    bad_temp = run_case(binary, root / "unsupported_temperature")
    assert bad_temp.returncode != 0
    assert "Temp/do_qhb/do_3tm" in bad_temp.stdout

    if args.gpu_binary:
        gpu_binary = args.gpu_binary.resolve()
        gpu_outputs: dict[str, str] = {}
        for name in ("gpu_static_mixed", "gpu_adaptive_mixed"):
            result = run_case(gpu_binary, root / name)
            if "no CUDA-capable device" in result.stdout or "no HIP-capable device" in result.stdout:
                print("CG-10.5 GPU executable cases skipped: no backend device")
                break
            assert result.returncode == 0, f"{name}\n{result.stdout}"
            assert "Gpu: AdaptiveCG initial active_atoms=" in result.stdout
            gpu_outputs[name] = result.stdout
        if len(gpu_outputs) == 2:
            static_active = metric(gpu_outputs["gpu_static_mixed"], "active_atoms")
            adaptive_counts = [
                int(value)
                for value in re.findall(
                    r"AdaptiveCG (?:initial|final) active_atoms=(\d+)",
                    gpu_outputs["gpu_adaptive_mixed"],
                )
            ]
            assert 0 < static_active < 48
            assert len(adaptive_counts) == 2 and adaptive_counts[1] < adaptive_counts[0]
            print("CG-10.5 GPU production executable tests passed")
    print("CG-10.5 production executable tests passed")


if __name__ == "__main__":
    main()
