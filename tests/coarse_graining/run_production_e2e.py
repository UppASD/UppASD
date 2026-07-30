#!/usr/bin/env python3
"""Launch the normal UppASD executable for CG-10.5 production cases."""

from __future__ import annotations

import argparse
import math
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


def float_metric(output: str, name: str) -> float:
    matches = re.findall(
        rf"{re.escape(name)}=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)",
        output,
    )
    if not matches:
        raise AssertionError(f"missing {name} in output\n{output}")
    return float(matches[-1])


def final_state(output: str) -> list[int]:
    matches = re.findall(
        r"(?:resolution_state label=final step=\d+|final_state) values=([0-9,\-]+)",
        output,
    )
    if not matches:
        raise AssertionError(f"missing final adaptive state\n{output}")
    return [int(value) for value in matches[-1].split(",")]


def close(reference: float, candidate: float, relative: float, absolute: float) -> bool:
    return math.isclose(reference, candidate, rel_tol=relative, abs_tol=absolute)


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
    examples = Path(__file__).resolve().parents[2] / "examples" / "AdaptiveCoarseGraining"

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
    bad_initial_phase = run_case(binary, root / "unsupported_initial_phase_x")
    assert bad_initial_phase.returncode != 0
    assert "ip_mode" in bad_initial_phase.stdout
    assert "Calls parallel tempering" not in bad_initial_phase.stdout

    parity_cpu: dict[str, str] = {}
    for mode in ("static", "adaptive"):
        name = f"parity_{mode}_cpu"
        result = run_case(binary, root / name)
        assert result.returncode == 0, f"{name}\n{result.stdout}"
        assert "last_energy_j atomistic_bilinear=" in result.stdout
        assert "last_field_checksums_t atom_sum=" in result.stdout
        assert "trajectory_checksums direction_sum=" in result.stdout
        parity_cpu[mode] = result.stdout
    assert abs(float_metric(parity_cpu["static"], "direction_sum")) < 40.0
    assert metric(parity_cpu["adaptive"], "accepted_transitions") > 0

    initial_phase_cpu = run_case(binary, root / "initial_phase_sd_cpu")
    assert initial_phase_cpu.returncode == 0, initial_phase_cpu.stdout
    assert "Performing initial phase:" in initial_phase_cpu.stdout
    assert "initial_state_source=completed_atomistic_sd" in initial_phase_cpu.stdout
    assert initial_phase_cpu.stdout.index("Performing initial phase:") < (
        initial_phase_cpu.stdout.index("AdaptiveCG: capability accepted")
    )
    assert metric(initial_phase_cpu.stdout, "accepted_transitions") > 0
    assert close(
        float_metric(initial_phase_cpu.stdout, "direction_norm2"),
        48.0,
        2.0e-12,
        2.0e-12,
    )
    print("CG-10.5 Initmag=2 atomistic SD initial-phase handoff passed")

    atomistic_initial_phases = {
        "initial_phase_mc": ("M", "Performing MC initial phase:"),
        "initial_phase_heat_bath": ("H", "Performing MC initial phase:"),
        "initial_phase_q": ("Q", "Line search minimization done"),
        "initial_phase_y": ("Y", "Line search minimization done"),
        "initial_phase_z": ("Z", "Line search minimization done"),
        "initial_phase_g": ("G", "Performing initial phase: energy minimization"),
    }
    atomistic_outputs: dict[str, str] = {}
    for name, (ip_mode, runner_banner) in atomistic_initial_phases.items():
        result = run_case(binary, root / name)
        assert result.returncode == 0, f"{name}\n{result.stdout}"
        assert runner_banner.lower() in result.stdout.lower()
        assert f"handoff_state_validated mode={ip_mode}" in result.stdout
        assert f"initial_state_source=completed_atomistic_{ip_mode}" in result.stdout
        assert result.stdout.index("Enter initial phase:") < result.stdout.index(
            "AdaptiveCG: capability accepted"
        )
        assert close(
            float_metric(result.stdout, "direction_norm2"),
            48.0,
            2.0e-12,
            2.0e-12,
        )
        atomistic_outputs[ip_mode] = result.stdout
    for ip_mode in ("Q", "Y", "Z"):
        assert abs(float_metric(atomistic_outputs[ip_mode], "direction_sum")) < 40.0
    print("CG-10.5 M/H/Q/Y/Z/G validated atomistic handoffs passed")

    spiral = run_case(binary, root / "initmag_spin_spiral")
    assert spiral.returncode == 0, spiral.stdout
    assert "spin spiral modulation" in spiral.stdout
    assert "initial_state_source=initmag mode=8" in spiral.stdout
    assert abs(float_metric(spiral.stdout, "direction_sum")) < 20.0
    print("CG-10.5 Initmag=8 inhomogeneous spin-spiral seed passed")

    for name in ("static_mixed", "adaptive", "initial_phase_texture"):
        result = run_case(binary, examples / name)
        assert result.returncode == 0, f"example {name}\n{result.stdout}"
        assert "AdaptiveCG: capability accepted" in result.stdout
        assert "AdaptiveCG: last_energy_j" in result.stdout
        if name == "initial_phase_texture":
            assert "initial_state_source=completed_atomistic_sd" in result.stdout
    print("CG-10.5 user-facing input examples passed")

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
            initial_phase_gpu = run_case(gpu_binary, root / "initial_phase_sd_gpu")
            assert initial_phase_gpu.returncode == 0, initial_phase_gpu.stdout
            assert "GpuSDSimulation: SD initial phase starting" in initial_phase_gpu.stdout
            assert "initial_state_source=completed_atomistic_sd" in initial_phase_gpu.stdout
            assert initial_phase_gpu.stdout.index(
                "GpuSDSimulation: SD initial phase starting"
            ) < initial_phase_gpu.stdout.index("AdaptiveCG: capability accepted")
            assert metric(initial_phase_gpu.stdout, "accepted_transitions") > 0
            mc_initial_phase_gpu = run_case(
                gpu_binary, root / "initial_phase_mc_gpu"
            )
            assert mc_initial_phase_gpu.returncode == 0, mc_initial_phase_gpu.stdout
            assert "GpuMCSimulation: MC initial phase starting" in (
                mc_initial_phase_gpu.stdout
            )
            assert "handoff_state_validated mode=M" in mc_initial_phase_gpu.stdout
            assert "initial_state_source=completed_atomistic_M" in (
                mc_initial_phase_gpu.stdout
            )
            assert mc_initial_phase_gpu.stdout.index(
                "GpuMCSimulation: MC initial phase starting"
            ) < mc_initial_phase_gpu.stdout.index("AdaptiveCG: capability accepted")
            q_initial_phase_gpu = run_case(
                gpu_binary, root / "initial_phase_q_gpu"
            )
            assert q_initial_phase_gpu.returncode == 0, q_initial_phase_gpu.stdout
            assert "line search minimization done" in q_initial_phase_gpu.stdout.lower()
            assert "handoff_state_validated mode=Q" in q_initial_phase_gpu.stdout
            assert "initial_state_source=completed_atomistic_Q" in (
                q_initial_phase_gpu.stdout
            )
            assert q_initial_phase_gpu.stdout.index(
                "Enter initial phase:"
            ) < q_initial_phase_gpu.stdout.index("AdaptiveCG: capability accepted")
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
            for mode in ("static", "adaptive"):
                name = f"parity_{mode}_gpu"
                result = run_case(gpu_binary, root / name)
                assert result.returncode == 0, f"{name}\n{result.stdout}"
                assert "Gpu: AdaptiveCG last_energy_j" in result.stdout
                assert "Gpu: AdaptiveCG last_field_checksums_t" in result.stdout
                cpu = parity_cpu[mode]
                gpu = result.stdout
                assert final_state(cpu) == final_state(gpu), (
                    f"{mode} CPU/GPU state decisions differ\n"
                    f"CPU {final_state(cpu)}\nGPU {final_state(gpu)}"
                )
                for key in (
                    "atomistic_bilinear",
                    "coarse_exchange",
                    "coarse_spiralization",
                    "coarse_anisotropy",
                    "coarse_external",
                    "coarse_dipole",
                    "total",
                ):
                    cpu_key = "last_total_energy_j" if key == "total" else key
                    reference = float_metric(cpu, cpu_key)
                    candidate = float_metric(gpu, key)
                    assert close(reference, candidate, 5.0e-4, 2.0e-24), (
                        f"{mode} CPU/GPU {key} energy differs: "
                        f"{reference} vs {candidate}"
                    )
                for key in ("atom_sum", "atom_norm2", "coarse_sum", "coarse_norm2"):
                    reference = float_metric(cpu, key)
                    candidate = float_metric(gpu, key)
                    assert close(reference, candidate, 8.0e-4, 2.0e-6), (
                        f"{mode} CPU/GPU {key} field checksum differs: "
                        f"{reference} vs {candidate}"
                    )
                cpu_state = restart_state(root / f"parity_{mode}_cpu")
                gpu_state = restart_state(root / f"parity_{mode}_gpu")
                assert len(cpu_state) == len(gpu_state) == 48
                assert max(
                    abs(reference - candidate)
                    for reference_row, candidate_row in zip(cpu_state, gpu_state)
                    for reference, candidate in zip(reference_row, candidate_row)
                ) <= 8.0e-5
                if mode == "adaptive":
                    assert metric(cpu, "accepted_transitions") == metric(
                        gpu, "accepted_transitions"
                    )
            fft = run_case(gpu_binary, root / "gpu_fft_static_mixed")
            assert fft.returncode == 0, f"gpu_fft_static_mixed\n{fft.stdout}"
            assert "EWALD3D_FFT production field/energy operator ready" in fft.stdout
            assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0
            assert float_metric(fft.stdout, "fft") > 0.0
            assert final_state(fft.stdout) == final_state(parity_cpu["static"])
            print("CG-10.5 adaptive uniform FFT dipole production coupling passed")
            print("CG-10.5 CPU/GPU state-field-energy-transition-trajectory parity passed")
            print("CG-10.5 GPU production executable tests passed")
    print("CG-10.5 production executable tests passed")


if __name__ == "__main__":
    main()
