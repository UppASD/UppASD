#!/usr/bin/env python3
"""RCG-06C: CPU/GPU adaptive-CG phase-timing reconciliation (F-18).

Reuses ``run_moving_backend_parity.py``'s fixture list, workspace
preparation, and single-fixture runner rather than inventing a new fixture
or copy mechanism. Runs ``moving_wall_adaptive`` -- the only tracked
fixture that engages every CPU/GPU adaptive-CG phase, including the
selector/transition path (ADAPTIVE mask mode; every other tracked fixture
is either feature-off or STATIC-mask, so their selector phase is
legitimately always zero) -- and parses the new phase-reconciliation
diagnostic lines this slice's fix added:

- CPU: ``AdaptiveCG: phase_wall_seconds wall=... phase_sum=... unaccounted=...``
  (``print_adaptive_cg_summary``, ``source/CoarseGraining/adaptivecgproduction.f90``).
- GPU: ``Gpu: AdaptiveCG final ... wall_ms=... phase_sum_ms=... unaccounted_ms=...``
  (``GpuSimulation::release``, ``source/gpu_files/gpuSimulation.cpp``).

What this proves, and what it does not:

- It proves the reconciliation identity holds on a real run
  (``unaccounted == wall - phase_sum`` as printed), that ``phase_sum`` is
  never more than a small clock-domain-noise margin above ``wall`` (a
  genuine correctness invariant -- named phases cannot legitimately take
  longer than the wall-clock interval that contains them), and that
  ``phase_sum`` is itself nonzero (proving phases are actually being
  measured, not silently wired to a stub).
- It does NOT construct a negative control that distinguishes
  ``system_clock`` (wall clock, what this slice's fix uses) from
  ``cpu_time`` (process CPU time, what it replaced) at the source level,
  because nothing in this codebase yet introduces OpenMP threading or
  blocking I/O inside the adaptive-CG CPU path that would make the two
  clocks diverge (RCG-07 is where CPU parallelism is planned) -- there is
  currently no runtime behavior difference to inject a fault against. The
  clock-semantics checklist item is verified by source inspection/diff
  instead (recorded in the evidence document), not by this fixture.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from run_moving_backend_parity import FIXTURES, prepare_workspace, run_fixture

FIXTURE_NAME = "moving_wall_adaptive"

_FLOAT_RE = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
_CPU_LINE_RE = re.compile(
    rf"AdaptiveCG: phase_wall_seconds wall=\s*({_FLOAT_RE}) phase_sum=\s*({_FLOAT_RE}) "
    rf"unaccounted=\s*({_FLOAT_RE})"
)
_GPU_LINE_RE = re.compile(
    rf"Gpu: AdaptiveCG final .*? wall_ms=({_FLOAT_RE}) phase_sum_ms=({_FLOAT_RE}) "
    rf"unaccounted_ms=({_FLOAT_RE})"
)


class ReconciliationError(AssertionError):
    pass


def check_cpu(stdout: str) -> None:
    match = _CPU_LINE_RE.search(stdout)
    if not match:
        raise ReconciliationError("CPU stdout has no 'phase_wall_seconds' reconciliation line")
    wall, phase_sum, unaccounted = (float(x) for x in match.groups())
    if wall <= 0.0:
        raise ReconciliationError(f"CPU wall time is not positive: {wall}")
    if phase_sum <= 0.0:
        raise ReconciliationError(f"CPU phase_sum is not positive: {phase_sum} -- phases are not being measured")
    if abs(unaccounted - (wall - phase_sum)) > 1.0e-9:
        raise ReconciliationError(
            f"CPU unaccounted={unaccounted} does not equal wall-phase_sum={wall - phase_sum} as printed"
        )
    if phase_sum > wall * 1.001:
        raise ReconciliationError(f"CPU phase_sum {phase_sum} exceeds wall {wall} by more than rounding noise")
    if unaccounted < -wall * 0.01:
        raise ReconciliationError(f"CPU unaccounted time {unaccounted} is materially negative")
    print(f"CPU reconciliation: wall={wall:.6f}s phase_sum={phase_sum:.6f}s "
          f"unaccounted={unaccounted:.6f}s ({100.0 * unaccounted / wall:.2f}% of wall time)")


def check_gpu(stdout: str) -> None:
    match = _GPU_LINE_RE.search(stdout)
    if not match:
        raise ReconciliationError("GPU stdout has no 'wall_ms=.../phase_sum_ms=.../unaccounted_ms=' reconciliation line")
    wall, phase_sum, unaccounted = (float(x) for x in match.groups())
    if wall <= 0.0:
        raise ReconciliationError(f"GPU wall time is not positive: {wall}")
    if phase_sum <= 0.0:
        raise ReconciliationError(f"GPU phase_sum is not positive: {phase_sum} -- phases are not being measured")
    # gpuSimulation.cpp prints this line with "%.3f" (3 decimals), and wall,
    # phase_sum, and unaccounted are each rounded independently from the
    # underlying doubles -- unlike the CPU line, which prints all three with
    # es24.16 (full double precision, no rounding to worry about, hence its
    # 1.0e-9 tolerance below). Three independent +/-0.0005ms roundings can
    # combine to a worst case of 0.0015ms; found by this line failing
    # intermittently (tolerance was 1.0e-6, tighter than the print format
    # itself supports) while collecting RCG-06D's regression evidence --
    # pre-existing since RCG-06C, unrelated to RCG-06D's RNG work.
    if abs(unaccounted - (wall - phase_sum)) > 2.0e-3:
        raise ReconciliationError(
            f"GPU unaccounted={unaccounted} does not equal wall_ms-phase_sum_ms={wall - phase_sum} as printed"
        )
    # Device-event timers (CUDA/HIP GPU_EVENT_ELAPSED_TIME) and the host wall
    # clock (std::chrono::steady_clock) are different clock domains with
    # independent measurement noise, so this allows more slack than the
    # CPU check's single-clock reconciliation.
    if phase_sum > wall * 1.05:
        raise ReconciliationError(
            f"GPU phase_sum {phase_sum} exceeds wall {wall} by more than expected cross-clock-domain noise"
        )
    if unaccounted < -wall * 0.05:
        raise ReconciliationError(f"GPU unaccounted time {unaccounted} is materially negative")
    print(f"GPU reconciliation: wall={wall:.3f}ms phase_sum={phase_sum:.3f}ms "
          f"unaccounted={unaccounted:.3f}ms ({100.0 * unaccounted / wall:.2f}% of wall time)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--gpu-binary", type=Path)
    parser.add_argument("--workspace-root", type=Path, required=True)
    args = parser.parse_args()

    args.workspace_root.mkdir(parents=True, exist_ok=True)
    fixture = next(f for f in FIXTURES if f.name == FIXTURE_NAME)

    cpu_workspace = prepare_workspace(args.workspace_root, "cpu", gpu=False)
    cpu_result = run_fixture(args.cpu_binary.resolve(), cpu_workspace, fixture)
    check_cpu(cpu_result.stdout)

    if args.gpu_binary:
        gpu_workspace = prepare_workspace(args.workspace_root, "gpu", gpu=True)
        gpu_result = run_fixture(args.gpu_binary.resolve(), gpu_workspace, fixture)
        check_gpu(gpu_result.stdout)
    else:
        print("GPU binary not supplied; GPU reconciliation not run (deferred, not blocking)")

    print("RCG-06C timing reconciliation: PASS")


if __name__ == "__main__":
    main()
