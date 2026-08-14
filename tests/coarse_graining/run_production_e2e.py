#!/usr/bin/env python3
"""Launch the normal UppASD executable for CG-10.5 production cases."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
from pathlib import Path

from fixture_dependencies import (
    ADAPTIVE_MIXED_CASE,
    CPU_ANISOTROPY_CASES,
    CPU_INITIAL_PHASE_CASES,
    CPU_INITIAL_PHASE_SD_CASE,
    CPU_PARITY_CASES,
    CPU_REJECTION_CASES,
    EXAMPLE_CASES,
    FEATURE_OFF_CASE,
    FINITE_TEMPERATURE_CASE,
    GPU_DMI_CASE,
    GPU_EXECUTABLE_CASES,
    GPU_FFT_CASE,
    GPU_PARITY_CASES,
    GPU_INITIAL_PHASE_CASES,
    POLARIZATION_GATE_CASE,
    POLARIZATION_GATE_GPU_CASE,
    RESTART_INITMAG_CASE,
    SPIN_SPIRAL_CASE,
    STATIC_ALL_COARSE_CASE,
    STATIC_ALL_FINE_CASE,
    STATIC_MIXED_CASE,
)


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


def static_stage_trace(output: str) -> list[tuple[int, str, dict[str, float]]]:
    """Return the diagnostic RCG-09B handoff snapshots in emission order."""
    traces: list[tuple[int, str, dict[str, float]]] = []
    pattern = re.compile(
        r"(?:Gpu:\s+)?AdaptiveCG: stage_trace step=(\d+) stage=(\w+) (.*)"
    )
    for match in pattern.finditer(output):
        values = {
            name: float(value)
            for name, value in re.findall(
                r"([a-z0-9_]+)=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)",
                match.group(3),
            )
        }
        traces.append((int(match.group(1)), match.group(2), values))
    return traces


def assert_static_stage_parity(cpu: str, gpu: str) -> None:
    """Fail at the first static Depondt/Heun handoff that diverges."""
    reference = static_stage_trace(cpu)
    candidate = static_stage_trace(gpu)
    assert reference, f"missing CPU static stage trace\\n{cpu}"
    assert len(reference) == len(candidate), (
        f"static stage trace count differs: CPU={len(reference)} GPU={len(candidate)}\\n"
        f"GPU output\\n{gpu}"
    )
    tolerances = {
        "atom_direction_sum": (8.0e-5, 8.0e-5),
        "atom_direction_norm2": (8.0e-5, 8.0e-5),
        "atomistic_direction_sum": (8.0e-5, 8.0e-5),
        "atomistic_direction_norm2": (8.0e-5, 8.0e-5),
        "coarse_direction_sum": (8.0e-5, 8.0e-5),
        "coarse_direction_norm2": (8.0e-5, 8.0e-5),
        "atom_field_sum": (8.0e-4, 2.0e-6),
        "atom_field_norm2": (8.0e-4, 2.0e-6),
        "coarse_field_sum": (8.0e-4, 2.0e-6),
        "coarse_field_norm2": (8.0e-4, 2.0e-6),
        "atomistic_bilinear": (5.0e-4, 2.0e-24),
    }
    for (step, stage, expected), (actual_step, actual_stage, actual) in zip(
        reference, candidate
    ):
        assert (step, stage) == (actual_step, actual_stage), (
            "static stage order differs: "
            f"CPU=({step}, {stage}) GPU=({actual_step}, {actual_stage})"
        )
        for metric_name, (relative, absolute) in tolerances.items():
            if stage == "corrector" and metric_name.startswith("atom_direction"):
                # The subset Depondt contract leaves non-atomistic entries in
                # emom2 untouched until the following reconstruction. Compare
                # its accepted active atoms here; full directions are checked
                # after the reconstructed stage.
                continue
            assert metric_name in expected and metric_name in actual, (
                f"static stage {step}/{stage} is missing {metric_name}"
            )
            assert close(expected[metric_name], actual[metric_name], relative, absolute), (
                f"static CPU/GPU first divergence at step={step} stage={stage} "
                f"metric={metric_name}: {expected[metric_name]} vs {actual[metric_name]}"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--gpu-binary", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = Path(__file__).with_name("e2e")
    examples = Path(__file__).resolve().parents[2] / "examples" / "AdaptiveCoarseGraining"

    off = run_case(binary, root / FEATURE_OFF_CASE)
    assert off.returncode == 0, off.stdout
    assert "AdaptiveCG:" not in off.stdout, off.stdout

    outputs: dict[str, str] = {}
    for name in (
        STATIC_ALL_FINE_CASE,
        STATIC_ALL_COARSE_CASE,
        STATIC_MIXED_CASE,
        ADAPTIVE_MIXED_CASE,
    ):
        result = run_case(binary, root / name)
        assert result.returncode == 0, f"{name}\n{result.stdout}"
        assert "AdaptiveCG: capability accepted" in result.stdout
        outputs[name] = result.stdout

    finite_temperature = run_case(binary, root / FINITE_TEMPERATURE_CASE)
    assert finite_temperature.returncode == 0, finite_temperature.stdout
    assert "AdaptiveCG: capability accepted" in finite_temperature.stdout
    print("CG-10.5 finite-temperature adaptive fine-region execution passed")

    fine_updates = metric(outputs[STATIC_ALL_FINE_CASE], "active_atom_updates")
    coarse_updates = metric(outputs[STATIC_ALL_COARSE_CASE], "active_atom_updates")
    mixed_updates = metric(outputs[STATIC_MIXED_CASE], "active_atom_updates")
    baseline = metric(outputs[STATIC_ALL_COARSE_CASE], "baseline_atom_updates")
    reduced = metric(outputs[STATIC_ALL_COARSE_CASE], "reduced_atom_updates")
    assert fine_updates == baseline
    assert coarse_updates < mixed_updates < fine_updates
    assert reduced > 0
    assert metric(outputs[ADAPTIVE_MIXED_CASE], "accepted_transitions") > 0
    off_state = restart_state(root / FEATURE_OFF_CASE)
    fine_state = restart_state(root / STATIC_ALL_FINE_CASE)
    assert len(off_state) == len(fine_state) == 48
    assert max(
        abs(reference - adaptive)
        for reference_row, adaptive_row in zip(off_state, fine_state)
        for reference, adaptive in zip(reference_row, adaptive_row)
    ) <= 1.0e-10

    bad = run_case(binary, root / CPU_REJECTION_CASES[0])
    assert bad.returncode != 0
    assert "block_size_x/y/z" in bad.stdout
    bad_mask = run_case(binary, root / CPU_REJECTION_CASES[1])
    assert bad_mask.returncode != 0
    assert "duplicate block id" in bad_mask.stdout
    bad_initial_phase = run_case(binary, root / CPU_REJECTION_CASES[2])
    assert bad_initial_phase.returncode != 0
    assert "ip_mode" in bad_initial_phase.stdout
    assert "Calls parallel tempering" not in bad_initial_phase.stdout
    missing_alat = run_case(binary, root / CPU_REJECTION_CASES[3])
    assert missing_alat.returncode != 0
    assert "explicit positive SI lattice parameter in metres" in missing_alat.stdout

    chiral = run_case(binary, root / CPU_ANISOTROPY_CASES[0])
    assert chiral.returncode == 0, chiral.stdout
    assert "AdaptiveCG: capability accepted" in chiral.stdout
    for term in (
        "atomistic_onsite",
        "coarse_spiralization",
        "coarse_anisotropy",
    ):
        assert abs(float_metric(chiral.stdout, term)) > 1.0e-30, (
            f"adaptive DMI/anisotropy term {term} was not exercised\n{chiral.stdout}"
        )
    print("CG-10.5 mixed DMI/anisotropy production path and SI alat guard passed")

    anisotropy_fine = run_case(binary, root / CPU_ANISOTROPY_CASES[1])
    anisotropy_coarse = run_case(binary, root / CPU_ANISOTROPY_CASES[2])
    assert anisotropy_fine.returncode == 0, anisotropy_fine.stdout
    assert anisotropy_coarse.returncode == 0, anisotropy_coarse.stdout
    expected_anisotropy_j = 24.0 * (-0.002 - 0.003) * 2.179872325e-21
    assert close(
        float_metric(anisotropy_fine.stdout, "atomistic_onsite"),
        expected_anisotropy_j,
        2.0e-12,
        1.0e-32,
    )
    assert close(
        float_metric(anisotropy_coarse.stdout, "coarse_anisotropy"),
        expected_anisotropy_j,
        2.0e-12,
        1.0e-32,
    )
    assert abs(float_metric(anisotropy_fine.stdout, "coarse_anisotropy")) < 1.0e-32
    assert abs(float_metric(anisotropy_coarse.stdout, "atomistic_onsite")) < 1.0e-32
    print("CG-10.5 uniaxial mRy-to-J and atomistic/coarse ownership passed")

    parity_cpu: dict[str, str] = {}
    for mode, name in CPU_PARITY_CASES.items():
        result = run_case(binary, root / name)
        assert result.returncode == 0, f"{name}\n{result.stdout}"
        assert "last_energy_j atomistic_bilinear=" in result.stdout
        assert "last_field_checksums_t atom_sum=" in result.stdout
        assert "trajectory_checksums direction_sum=" in result.stdout
        parity_cpu[mode] = result.stdout
    assert abs(float_metric(parity_cpu["static"], "direction_sum")) < 40.0
    assert metric(parity_cpu["adaptive"], "accepted_transitions") > 0

    # RCG-03 POL-CANCELLATION production regression: identical geometry and a
    # selector configuration deliberately neutralized to always request
    # coarsening (cg_refine_threshold/cg_coarsen_threshold saturated near the
    # misalignment score's 2.0 ceiling) as parity_adaptive_cpu above, but
    # starting from a genuinely incoherent (random) per-atom state instead of
    # an aligned one. With the misalignment criterion neutralized, the only
    # thing that can still block coarsening is the polarization gate itself;
    # parity_adaptive_cpu (aligned, same neutralized thresholds) coarsening
    # while this case does not is the negative control proving the gate
    # (rather than the ordinary selector) is what holds the line here.
    polarization_gate = run_case(binary, root / POLARIZATION_GATE_CASE)
    assert polarization_gate.returncode == 0, polarization_gate.stdout
    assert metric(polarization_gate.stdout, "accepted_transitions") == 0
    final_coarse_counts = re.findall(
        r"resolution_counts fine=\d+ interface=\d+ coarse=(\d+)", polarization_gate.stdout
    )
    assert final_coarse_counts and int(final_coarse_counts[-1]) == 0, polarization_gate.stdout
    assert all(state != 0 for state in final_state(polarization_gate.stdout))
    print("RCG-03 polarization gate blocks coarsening of a low-polarization block")

    initial_phase_cpu = run_case(binary, root / CPU_INITIAL_PHASE_SD_CASE)
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
        case: (mode, banner)
        for mode, (case, banner) in CPU_INITIAL_PHASE_CASES.items()
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

    spiral = run_case(binary, root / SPIN_SPIRAL_CASE)
    assert spiral.returncode == 0, spiral.stdout
    assert "spin spiral modulation" in spiral.stdout
    assert "initial_state_source=initmag mode=8" in spiral.stdout
    assert abs(float_metric(spiral.stdout, "direction_sum")) < 20.0
    print("CG-10.5 Initmag=8 inhomogeneous spin-spiral seed passed")

    # RCG-04B follow-up: validate_configuration used to reject initmag=4
    # (restart) outright. It no longer does; this is a setup/capability
    # smoke test only (Nstep 1) proving AdaptiveCG accepts a restart-format
    # initial condition end to end, not a moving-dynamics or parity claim.
    restart_initmag = run_case(binary, root / RESTART_INITMAG_CASE)
    assert restart_initmag.returncode == 0, restart_initmag.stdout
    assert "AdaptiveCG: capability accepted" in restart_initmag.stdout
    assert "initial_state_source=initmag mode=4" in restart_initmag.stdout
    print("CG-10.5 Initmag=4 restart capability passed")

    for name in EXAMPLE_CASES:
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
        for name in GPU_EXECUTABLE_CASES:
            result = run_case(gpu_binary, root / name)
            if "no CUDA-capable device" in result.stdout or "no HIP-capable device" in result.stdout:
                print("CG-10.5 GPU executable cases skipped: no backend device")
                break
            assert result.returncode == 0, f"{name}\n{result.stdout}"
            assert "Gpu: AdaptiveCG initial active_atoms=" in result.stdout
            gpu_outputs[name] = result.stdout
        if len(gpu_outputs) == 2:
            initial_phase_gpu = run_case(gpu_binary, root / GPU_INITIAL_PHASE_CASES[0])
            assert initial_phase_gpu.returncode == 0, initial_phase_gpu.stdout
            assert "GpuSDSimulation: SD initial phase starting" in initial_phase_gpu.stdout
            assert "initial_state_source=completed_atomistic_sd" in initial_phase_gpu.stdout
            assert initial_phase_gpu.stdout.index(
                "GpuSDSimulation: SD initial phase starting"
            ) < initial_phase_gpu.stdout.index("AdaptiveCG: capability accepted")
            assert metric(initial_phase_gpu.stdout, "accepted_transitions") > 0
            mc_initial_phase_gpu = run_case(
                gpu_binary, root / GPU_INITIAL_PHASE_CASES[1]
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
                gpu_binary, root / GPU_INITIAL_PHASE_CASES[2]
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
            for mode, name in GPU_PARITY_CASES.items():
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
                if mode == "static":
                    assert_static_stage_parity(cpu, gpu)
                for key in (
                    "atomistic_bilinear",
                    "atomistic_onsite",
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
                restart_max_error = max(
                    abs(reference - candidate)
                    for reference_row, candidate_row in zip(cpu_state, gpu_state)
                    for reference, candidate in zip(reference_row, candidate_row)
                )
                assert restart_max_error <= 8.0e-5, (
                    f"{mode} CPU/GPU restart state differs: "
                    f"max_abs={restart_max_error:.16e} tolerance=8.0e-5"
                )
                if mode == "adaptive":
                    assert metric(cpu, "accepted_transitions") == metric(
                        gpu, "accepted_transitions"
                    )

            # RCG-03 POL-CANCELLATION production regression, GPU backend:
            # same neutralized-selector / random-Initmag construction as the
            # CPU polarization_gate_cpu case above, run through the real GPU
            # dispatch path. Without the GPU-side gate (evaluatePolarizationGate
            # wired into proposeSelectorState's hardAtomisticBlockMask), this
            # would coarsen every step since the misalignment criterion is
            # saturated to always request it.
            polarization_gate_gpu = run_case(gpu_binary, root / POLARIZATION_GATE_GPU_CASE)
            assert polarization_gate_gpu.returncode == 0, polarization_gate_gpu.stdout
            assert metric(polarization_gate_gpu.stdout, "accepted_transitions") == 0
            gpu_final_coarse_counts = re.findall(
                r"resolution_counts fine=\d+ interface=\d+ coarse=(\d+)",
                polarization_gate_gpu.stdout,
            )
            assert gpu_final_coarse_counts and int(gpu_final_coarse_counts[-1]) == 0, (
                polarization_gate_gpu.stdout
            )
            print("RCG-03 GPU polarization gate blocks coarsening of a low-polarization block")

            chiral_gpu = run_case(gpu_binary, root / GPU_DMI_CASE)
            assert chiral_gpu.returncode == 0, chiral_gpu.stdout
            assert final_state(chiral.stdout) == final_state(chiral_gpu.stdout)
            for key in (
                "atomistic_bilinear",
                "atomistic_onsite",
                "coarse_exchange",
                "coarse_spiralization",
                "coarse_anisotropy",
                "coarse_external",
                "coarse_dipole",
                "total",
            ):
                cpu_key = "last_total_energy_j" if key == "total" else key
                reference = float_metric(chiral.stdout, cpu_key)
                candidate = float_metric(chiral_gpu.stdout, key)
                assert close(reference, candidate, 5.0e-4, 2.0e-24), (
                    f"DMI/anisotropy CPU/GPU {key} differs: "
                    f"{reference} vs {candidate}"
                )
            for key in ("atom_sum", "atom_norm2", "coarse_sum", "coarse_norm2"):
                reference = float_metric(chiral.stdout, key)
                candidate = float_metric(chiral_gpu.stdout, key)
                assert close(reference, candidate, 8.0e-4, 2.0e-6), (
                    f"DMI/anisotropy CPU/GPU {key} differs: "
                    f"{reference} vs {candidate}"
                )
            print("CG-10.5 DMI/anisotropy CPU/GPU production parity passed")
            fft = run_case(gpu_binary, root / GPU_FFT_CASE)
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
