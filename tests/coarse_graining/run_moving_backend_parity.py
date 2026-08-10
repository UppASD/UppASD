#!/usr/bin/env python3
"""RCG-04I: CPU/GPU backend precision parity across the accepted moving fixtures.

Scope (restated from ``docs/RCG-04_MOVING_E2E_PROMPT_PACK.md`` section 13):
this slice proves that the RCG-04D-H moving fixtures -- already accepted as
CPU-only evidence, each explicitly annotated in ``CMakeLists.txt`` with
"CPU-only in this slice; CPU/GPU backend parity is RCG-04I's responsibility,
not claimed here" -- are backend-equivalent for the *geometries already
accepted there*. It does not claim the unequal-width/skew-cell ownership
equivalence reserved for RCG-05.

**Mechanism, verified once by hand before writing this harness (see the
RCG-04I evidence document for the exact smoke-test transcript):** a single
CUDA-enabled build (``UPPASD_GPU_BACKEND=CUDA``) produces one executable
that runs the *ordinary* Fortran/CPU LLG path unless the input file also
requests ``do_gpu Y``/``do_gpu_llg Y`` (``source/uppasd.f90``, existing
``gpu_static_mixed``/``gpu_adaptive_mixed`` fixtures already use exactly
this flag block). None of the RCG-04D-H moving fixtures set ``do_gpu`` at
all, which is *why* they were CPU-only even when built with CUDA support.
This harness therefore does not need a second, backend-specific fixture: it
copies each accepted fixture's byte-identical inputs into an isolated
workspace per backend and appends the same four-line GPU-enable block
(never modifying momfile, jfile, posfile, mask.dat, kfile, or dmfile), then
runs the appropriate binary. "CPU fp64", "CUDA fp64", and "CUDA fp32" are
still three genuinely separate, freshly built executables (``UPPASD_PRECISION``
is a compile-time storage-type switch, not a runtime flag) -- see the
evidence document for each build's exact configure command and commit.

**Known, documented backend capability gap (not a defect introduced or
hidden here; restated from ``trajectory_evidence.py``'s module docstring,
first found by RCG-04C):** the GPU AdaptiveCG diagnostic path
(``source/gpu_files/gpuSimulation.cpp``) prints only ``initial``/``final``
resolution-state snapshots, never a per-event transition log, and its one
aggregate summary line hardcodes ``rejected_transitions=0``. This harness
therefore compares CPU and GPU ``accepted_transitions`` aggregate counts
(both backends track that one correctly) but does not claim per-event
transition-history parity on GPU, and explicitly records the
``rejected_transitions`` hardcode as an open capability gap rather than
silently comparing 0 to 0.

Every backend named on the command line is run against every fixture; a
backend simply not passed on the command line is never claimed as passing
-- see ``main()``'s final open-items report.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import trajectory_evidence as te
from run_moving_adaptive_wall import (
    ANISOTROPY_AXIS as WALL_AXIS,
    AXIS_LENGTH_CELLS as WALL_AXIS_LENGTH_CELLS,
    BLOCK_LENGTH_CELLS as WALL_BLOCK_LENGTH_CELLS,
    GEOMETRY as WALL_GEOMETRY,
)
from run_moving_dmi_chiral import CHIRALITY_AXIS, GEOMETRY as DMI_GEOMETRY

SOURCE_E2E = Path(__file__).with_name("e2e")

GPU_ENABLE_BLOCK = (
    "\ngpu_mode 1\ndo_gpu Y\ndo_gpu_llg Y\ndo_gpu_measurements N\ndo_gpu_correlations N\n"
)

_ACCEPTED_RE = re.compile(r"\baccepted_transitions=(\d+)")
_REJECTED_RE = re.compile(r"\brejected_transitions=(\d+)")


@dataclasses.dataclass(frozen=True)
class Fixture:
    name: str
    simid: str
    natom: int
    nstep: int
    adaptive_cg: bool
    chirality: bool = False
    wall: bool = False


# Every fixture RCG-04D-H accepted as CPU-only moving-dynamics evidence
# (source of truth: fixture_dependencies.py's MOVING_* constants). Two
# natural size/length scales already exist in the accepted set and are
# reused rather than inventing new ones: 48 atoms (moving_feature_off /
# moving_all_fine, ncell 6 2 2) vs 192 atoms (every other fixture,
# ncell 24 2 2), and 50 steps vs 900 steps (moving_wall_* ).
FIXTURES: tuple[Fixture, ...] = (
    Fixture("moving_feature_off", "cg104mvo", 48, 50, adaptive_cg=False),
    Fixture("moving_all_fine", "cg104mva", 48, 50, adaptive_cg=True),
    Fixture("moving_feature_off_wide", "cg104ewo", 192, 50, adaptive_cg=False),
    Fixture("moving_all_fine_wide", "cg104ewf", 192, 50, adaptive_cg=True),
    Fixture("moving_all_coarse_bs1", "cg104ec1", 192, 50, adaptive_cg=True),
    Fixture("moving_all_coarse_bs2", "cg104ec2", 192, 50, adaptive_cg=True),
    Fixture("moving_all_coarse_bs4", "cg104ec4", 192, 50, adaptive_cg=True),
    Fixture("moving_all_coarse_bs8", "cg104ec8", 192, 50, adaptive_cg=True),
    Fixture("moving_static_mixed_bs1", "cg104fb1", 192, 50, adaptive_cg=True),
    Fixture("moving_static_mixed_bs2", "cg104fb2", 192, 50, adaptive_cg=True),
    Fixture("moving_static_mixed_bs1_shifted", "cg104fs1", 192, 50, adaptive_cg=True),
    Fixture("moving_wall_feature_off", "cg105gof", 192, 900, adaptive_cg=False, wall=True),
    Fixture("moving_wall_adaptive", "cg105gad", 192, 900, adaptive_cg=True, wall=True),
    Fixture("moving_dmi_chiral_all_fine_plus", "cg104h1a", 192, 50, adaptive_cg=True, chirality=True),
    Fixture("moving_dmi_chiral_bs1_plus", "cg104h2a", 192, 50, adaptive_cg=True, chirality=True),
    Fixture("moving_dmi_chiral_bs1_minus", "cg104h2b", 192, 50, adaptive_cg=True, chirality=True),
    Fixture("moving_dmi_chiral_bs2_plus", "cg104h3a", 192, 50, adaptive_cg=True, chirality=True),
    Fixture("moving_dmi_chiral_bs1_plus_reversed", "cg104h4a", 192, 50, adaptive_cg=True, chirality=True),
    Fixture("moving_dmi_chiral_bs1_minus_reversed", "cg104h4b", 192, 50, adaptive_cg=True, chirality=True),
)

# --- Frozen fp64/fp32 backend-parity budgets -------------------------------
# Derived from this slice's own exploratory run (``--mode explore``) across
# every fixture above, i.e. two system sizes (48/192 atoms) and two
# trajectory lengths (50/900 steps); see the RCG-04I evidence document for
# the complete raw per-fixture table these were derived from. Frozen BEFORE
# the acceptance run below (governing rule: freeze before the final
# acceptance run, then re-verify negative controls under the frozen value).
#
# Observed scaling (raw table in the evidence document), by fixture class:
#
# - Ordinary atomistic GPU LLG (do_adaptive_cg=N: moving_feature_off[_wide],
#   moving_wall_feature_off): CPU/GPU differ by ~1e-6 to ~1.3e-4 rad, growing
#   with trajectory length (50 -> 900 steps), and this floor is essentially
#   IDENTICAL between the fp64 and fp32 GPU builds (e.g. moving_feature_off:
#   2.113e-5 rad fp64 vs 2.137e-5 rad fp32) -- consistent with a genuine
#   CPU/GPU floating-point-associativity (summation-order) difference at
#   double precision dominating, with fp32's additional truncation being
#   negligible on top of it.
# - STATIC-mask AdaptiveCG (every moving_all_fine/all_coarse/static_mixed/
#   dmi_chiral fixture): fp64 CPU vs GPU trajectories are exactly bit-
#   identical (0.0 error) over all 12 such fixtures; fp32 shows a small
#   nonzero ~1e-7-6e-7 rad error consistent with pure fp32 storage
#   truncation. See this file's ``load_run_result``/``run_backend``
#   docstring and the evidence document for why exact fp64 identity is
#   plausible here (small per-block reductions execute in the same order on
#   both backends at this system size) and is flagged as an open item for a
#   larger-scale re-check, not treated as a stronger or weaker result than
#   the other classes.
# - ADAPTIVE-mask AdaptiveCG (moving_wall_adaptive, the only fixture that
#   engages the genuine GPU adaptive-transition runtime,
#   ``gpuAdaptiveRuntime``): the largest observed divergence, ~0.062 rad
#   max angular error -- but, notably, *this is essentially identical
#   between the fp64 and fp32 GPU builds* (0.061900 rad fp64 vs 0.061900 rad
#   fp32), even though ``accepted_transitions`` and every wall
#   block-crossing event match exactly (same step, same from/to block) on
#   both backends and both precisions. This rules out a floating-point-
#   precision explanation for this fixture's error floor: it is dominated
#   by a genuine CPU/GPU algorithmic-order difference in the adaptive
#   selector/reconstruction path (not a missing transition, not a
#   precision-scaling effect), which is why the two budgets below are
#   deliberately close rather than differing by fp32/fp64's ~1e8x unit-
#   roundoff ratio -- fitting a large fp32/fp64 gap here would misrepresent
#   what was actually observed.
#
# A single flat budget across all 19 fixtures was tried first and rejected
# by this slice's own negative-control run (governing rule: re-verify
# defect sensitivity under the frozen budget): a budget loose enough to
# admit ``moving_wall_adaptive``'s genuine ~0.062 rad backend divergence
# (5x headroom -> 0.3 rad) is *larger than several fixtures' own total
# physical displacement* (e.g. ``moving_feature_off_wide``/
# ``moving_all_fine_wide`` move only ~0.0074 rad end to end -- see the
# displacement table in the evidence document), so a defect that produced
# total nonsense on one of those fixtures could still hide under a
# wall-adaptive-sized budget. The budget is therefore split by the fixture
# class that actually produces the observed error floor, each still
# checked to stay well under its class's own smallest observed physical
# displacement (so a freeze/sign/missing-transition defect on any covered
# fixture is still caught with margin):
#
# - "ordinary" (every fixture except ``moving_wall_adaptive``): observed
#   worst-case backend divergence 2.11e-5 rad fp64 / 2.14e-5 rad fp32
#   (``moving_feature_off``); observed smallest own physical displacement
#   in this class is 0.00742 rad (``moving_feature_off_wide``/
#   ``moving_all_fine_wide``). Budget below sits ~50-70x above the observed
#   noise floor and ~7x below the smallest displacement in the class.
# - "wall_adaptive" (``moving_wall_adaptive`` only, the sole fixture that
#   engages the genuine GPU adaptive-transition runtime): observed ~0.062
#   rad divergence, essentially identical at fp64 and fp32 (see above --
#   dominated by CPU/GPU algorithmic-order difference in the adaptive
#   selector/reconstruction path, not by storage precision); its own
#   physical displacement is 0.725 rad. Budget below gives ~2.4x headroom
#   over the observed divergence and stays ~4.8x below its own
#   displacement.
FROZEN_BUDGET_FP64_ORDINARY = dict(max_angular_rad=1.0e-3, max_component=1.0e-3, max_restart=1.0e-3)
FROZEN_BUDGET_FP64_WALL_ADAPTIVE = dict(max_angular_rad=1.5e-1, max_component=1.5e-1, max_restart=1.5e-1)
FROZEN_BUDGET_FP32_ORDINARY = dict(max_angular_rad=1.5e-3, max_component=1.5e-3, max_restart=1.5e-3)
FROZEN_BUDGET_FP32_WALL_ADAPTIVE = dict(max_angular_rad=1.8e-1, max_component=1.8e-1, max_restart=1.8e-1)

WALL_ADAPTIVE_FIXTURE_NAME = "moving_wall_adaptive"


def budget_for(precision: str, fixture: Fixture) -> dict:
    if fixture.name == WALL_ADAPTIVE_FIXTURE_NAME:
        return FROZEN_BUDGET_FP64_WALL_ADAPTIVE if precision == "fp64" else FROZEN_BUDGET_FP32_WALL_ADAPTIVE
    return FROZEN_BUDGET_FP64_ORDINARY if precision == "fp64" else FROZEN_BUDGET_FP32_ORDINARY


class BackendParityError(AssertionError):
    pass


class NegativeControlDidNotFailError(AssertionError):
    """Raised when a deliberately broken/mismatched trajectory pair passed the budget anyway."""


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_workspace(root: Path, backend: str, *, gpu: bool) -> Path:
    """Copy the whole ``e2e`` tree (fixtures + shared jfile/posfile/kfile/dmfile)
    into an isolated per-backend workspace, then (for a GPU backend) append the
    GPU-enable block to every listed fixture's ``inpsd.dat`` copy only.

    Copying the full tree (not just one fixture directory) is required
    because every fixture's ``inpsd.dat`` references shared inputs by a
    parent-relative path (``posfile ../posfile``, ``exchange ../jfile``,
    DMI fixtures additionally ``../dmfile_chiral``/``../kfile_cg_x``) --
    confirmed by grep across every ``moving_*/inpsd.dat`` before writing
    this function.
    """
    workspace = root / backend / "e2e"
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE_E2E, workspace)

    for fixture in FIXTURES:
        src_inpsd = SOURCE_E2E / fixture.name / "inpsd.dat"
        dst_inpsd = workspace / fixture.name / "inpsd.dat"
        src_bytes = src_inpsd.read_bytes()
        assert dst_inpsd.read_bytes() == src_bytes, (
            f"{backend}/{fixture.name}: inpsd.dat copy diverged from source before any "
            "GPU-enable block was appended"
        )
        src_mom = SOURCE_E2E / fixture.name / "momfile"
        dst_mom = workspace / fixture.name / "momfile"
        assert sha256_file(src_mom) == sha256_file(dst_mom), (
            f"{backend}/{fixture.name}: momfile is not byte-identical to the source fixture"
        )
        if gpu:
            with dst_inpsd.open("a") as handle:
                handle.write(GPU_ENABLE_BLOCK)
    return workspace


def run_fixture(binary: Path, workspace: Path, fixture: Fixture) -> subprocess.CompletedProcess[str]:
    case_dir = workspace / fixture.name
    result = subprocess.run(
        [str(binary)], cwd=case_dir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180, check=False,
    )
    if result.returncode != 0:
        raise BackendParityError(
            f"{fixture.name}: binary {binary} exited {result.returncode}:\n{result.stdout}"
        )
    return result


@dataclasses.dataclass
class FixtureRunResult:
    fixture: Fixture
    stdout: str
    workspace: Path
    trajectory: te.Trajectory
    restart: te.Trajectory


def load_run_result(
    fixture: Fixture, workspace: Path, stdout: str, *, normalization_tolerance: float
) -> FixtureRunResult:
    case_dir = workspace / fixture.name
    trajectory = te.load_moment_trajectory(
        case_dir, simid=fixture.simid, normalization_tolerance=normalization_tolerance
    )
    restart = te.load_restart_state(case_dir, normalization_tolerance=normalization_tolerance)
    return FixtureRunResult(
        fixture=fixture, stdout=stdout, workspace=workspace, trajectory=trajectory, restart=restart
    )


def run_backend(binary: Path, workspace: Path, *, normalization_tolerance: float) -> dict[str, FixtureRunResult]:
    """``normalization_tolerance`` is the strict parser's accepted deviation of a
    parsed direction vector's norm from 1.0 (``trajectory_evidence.py``'s own
    parser-level sanity check, unrelated to this harness's own CPU/GPU
    component/angular budgets below). fp64 keeps the parser's tight default
    (1e-6); fp32 needs a looser one -- fp32 has roughly 7 significant decimal
    digits, so a normalization drift of a few 1e-6 to 1e-5 after many
    accumulated LLG steps is the expected fp32 noise floor, not a defect
    (confirmed here: this harness first hit ``NonUnitDirectionError`` at
    ``norm=1.0000011...`` on a 900-step fp32 fixture, i.e. ~1.1e-6 drift --
    well inside fp32's expected precision, not a corrupted state).
    """
    results: dict[str, FixtureRunResult] = {}
    for fixture in FIXTURES:
        proc = run_fixture(binary, workspace, fixture)
        results[fixture.name] = load_run_result(
            fixture, workspace, proc.stdout, normalization_tolerance=normalization_tolerance
        )
    return results


def aggregate_counts(stdout: str) -> tuple[int | None, int | None]:
    accepted = _ACCEPTED_RE.search(stdout)
    rejected = _REJECTED_RE.search(stdout)
    return (
        int(accepted.group(1)) if accepted else None,
        int(rejected.group(1)) if rejected else None,
    )


@dataclasses.dataclass
class FixtureComparison:
    fixture: Fixture
    component_max: float
    component_rms: float
    angular_max_rad: float
    angular_rms_rad: float
    restart_max: float
    energy_terms: dict | None
    transitions_note: str
    chirality_cpu: float | None = None
    chirality_gpu: float | None = None
    wall_note: str | None = None

    def as_dict(self) -> dict:
        return {
            "fixture": self.fixture.name, "natom": self.fixture.natom, "nstep": self.fixture.nstep,
            "component_max": self.component_max, "component_rms": self.component_rms,
            "angular_max_rad": self.angular_max_rad, "angular_rms_rad": self.angular_rms_rad,
            "restart_max": self.restart_max, "energy_terms": self.energy_terms,
            "transitions_note": self.transitions_note,
            "chirality_cpu": self.chirality_cpu, "chirality_gpu": self.chirality_gpu,
            "wall_note": self.wall_note,
        }


def _restart_content_only(trajectory: te.Trajectory) -> te.Trajectory:
    """Normalize a restart trajectory's step label to 0 before comparison.

    **Finding, reported here per governing rule 3.1, not fixed (out of
    RCG-04I's validation-only scope):** ``source/uppasd.f90:363`` computes
    the restart file's ``mstep`` label as ``mstep = mstep - 1`` from the
    host Fortran step counter. Under ``do_gpu_llg=Y`` the measurement loop
    runs on the GPU and this host counter is never incremented, so every
    GPU-backend restart file (independent of fixture or ``Nstep``) is
    written with the sentinel label ``mstep=-1`` instead of the true final
    iteration -- confirmed by direct inspection of
    ``moving_feature_off``'s CPU (``mstep=50``) vs GPU (``mstep=-1``)
    restart file in this slice's evidence document. The physical
    ``|Mom|/Mx/My/Mz`` columns are unaffected (verified: their values are
    within the same fp64 numerical-parity range as the moment trajectory,
    not corrupted or zeroed) -- this is a diagnostic-label gap, not a
    physics defect, so it is normalized away for the content comparison
    below and reported as an open item, not silently masked.
    """
    normalized_steps = tuple(dataclasses.replace(step, step=0) for step in trajectory.steps)
    return dataclasses.replace(trajectory, steps=normalized_steps)


def compare_fixture(cpu: FixtureRunResult, gpu: FixtureRunResult) -> FixtureComparison:
    fixture = cpu.fixture
    component = te.component_trajectory_error(cpu.trajectory, gpu.trajectory)
    angular = te.angular_trajectory_error(cpu.trajectory, gpu.trajectory)
    restart_component = te.component_trajectory_error(
        _restart_content_only(cpu.restart), _restart_content_only(gpu.restart)
    )

    energy_terms = None
    if fixture.adaptive_cg:
        cpu_series = te.parse_energy_field_series(cpu.stdout)
        gpu_series = te.parse_energy_field_series(gpu.stdout)
        if cpu_series and gpu_series:
            energy_terms = te.compare_energy_field_series(cpu_series, gpu_series)

    transitions_note = "n/a (feature-off, no AdaptiveCG diagnostic)"
    if fixture.adaptive_cg:
        cpu_acc, cpu_rej = aggregate_counts(cpu.stdout)
        gpu_acc, gpu_rej = aggregate_counts(gpu.stdout)
        cpu_events = te.parse_transition_events(cpu.stdout)
        gpu_events = te.parse_transition_events(gpu.stdout)
        transitions_note = (
            f"accepted_transitions cpu={cpu_acc} gpu={gpu_acc}; "
            f"rejected_transitions cpu={cpu_rej} gpu={gpu_rej} "
            f"(GPU hardcodes 0 -- documented capability gap, not compared as parity); "
            f"per-event transition log: cpu={len(cpu_events)} events, gpu={len(gpu_events)} events "
            "(GPU has no per-event log -- documented capability gap, not claimed as parity)"
        )
        if cpu_acc is not None and gpu_acc is not None and cpu_acc != gpu_acc:
            transitions_note += " *** accepted_transitions MISMATCH ***"

    chirality_cpu = chirality_gpu = None
    if fixture.chirality:
        bonds = te.axis_chain_bonds(DMI_GEOMETRY, axis_cell_index=0)
        chirality_cpu = te.signed_chirality(cpu.trajectory.steps[-1], directed_bonds=bonds, axis=CHIRALITY_AXIS)
        chirality_gpu = te.signed_chirality(gpu.trajectory.steps[-1], directed_bonds=bonds, axis=CHIRALITY_AXIS)

    wall_note = None
    if fixture.wall:
        def wall_summary(run: FixtureRunResult) -> te.WallCrossingSummary:
            series = [
                (step.step, te.domain_wall_centers(
                    step, WALL_GEOMETRY, axis_cell_index=0, easy_axis=WALL_AXIS))
                for step in run.trajectory.steps
            ]
            return te.track_wall_crossings(
                series, axis_length_cells=WALL_AXIS_LENGTH_CELLS,
                block_length_cells=WALL_BLOCK_LENGTH_CELLS,
            )
        cpu_wall = wall_summary(cpu)
        gpu_wall = wall_summary(gpu)
        wall_note = (
            f"crossings cpu={cpu_wall.crossing_events} gpu={gpu_wall.crossing_events}; "
            f"displacement cpu={cpu_wall.displacement} gpu={gpu_wall.displacement}"
        )
        if cpu_wall.crossing_events != gpu_wall.crossing_events:
            wall_note += " *** CROSSING-EVENT MISMATCH ***"

    return FixtureComparison(
        fixture=fixture, component_max=component.max_abs_error_overall,
        component_rms=component.rms_error_overall, angular_max_rad=angular.max_radians,
        angular_rms_rad=angular.rms_radians, restart_max=restart_component.max_abs_error_overall,
        energy_terms=energy_terms, transitions_note=transitions_note,
        chirality_cpu=chirality_cpu, chirality_gpu=chirality_gpu, wall_note=wall_note,
    )


def print_comparison_table(precision: str, comparisons: list[FixtureComparison]) -> None:
    print(f"\n--- {precision} backend-parity results (CPU fp64 reference) ---")
    for comparison in comparisons:
        f = comparison.fixture
        print(
            f"  {f.name:38s} natom={f.natom:3d} nstep={f.nstep:4d} "
            f"component_max={comparison.component_max:.6g} "
            f"angular_max_rad={comparison.angular_max_rad:.6g} "
            f"({math.degrees(comparison.angular_max_rad):.4g} deg) "
            f"restart_max={comparison.restart_max:.6g}"
        )
        print(f"      transitions: {comparison.transitions_note}")
        if comparison.wall_note:
            print(f"      wall: {comparison.wall_note}")
        if comparison.chirality_cpu is not None:
            print(
                f"      chirality: cpu={comparison.chirality_cpu:.6g} gpu={comparison.chirality_gpu:.6g}"
            )


def assert_budget(precision: str, comparisons: list[FixtureComparison]) -> None:
    failures = []
    for comparison in comparisons:
        budget = budget_for(precision, comparison.fixture)
        if comparison.component_max > budget["max_component"]:
            failures.append(f"{comparison.fixture.name}: component_max {comparison.component_max} > {budget['max_component']}")
        if comparison.angular_max_rad > budget["max_angular_rad"]:
            failures.append(f"{comparison.fixture.name}: angular_max_rad {comparison.angular_max_rad} > {budget['max_angular_rad']}")
        if comparison.restart_max > budget["max_restart"]:
            failures.append(f"{comparison.fixture.name}: restart_max {comparison.restart_max} > {budget['max_restart']}")
        if comparison.fixture.chirality:
            same_sign = (comparison.chirality_cpu >= 0) == (comparison.chirality_gpu >= 0)
            if not same_sign:
                failures.append(f"{comparison.fixture.name}: chirality sign mismatch cpu={comparison.chirality_cpu} gpu={comparison.chirality_gpu}")
        if comparison.fixture.wall and "MISMATCH" in (comparison.wall_note or ""):
            failures.append(f"{comparison.fixture.name}: {comparison.wall_note}")
        if "MISMATCH" in comparison.transitions_note:
            failures.append(f"{comparison.fixture.name}: {comparison.transitions_note}")
    if failures:
        raise BackendParityError(f"{precision}: {len(failures)} fixture(s) exceeded the frozen budget:\n" + "\n".join(failures))
    ordinary = FROZEN_BUDGET_FP64_ORDINARY if precision == "fp64" else FROZEN_BUDGET_FP32_ORDINARY
    wall = FROZEN_BUDGET_FP64_WALL_ADAPTIVE if precision == "fp64" else FROZEN_BUDGET_FP32_WALL_ADAPTIVE
    print(f"{precision}: all {len(comparisons)} fixtures within their frozen class budget "
          f"(ordinary={ordinary}, wall_adaptive={wall})")


def run_negative_controls(precision: str, cpu_results: dict[str, FixtureRunResult],
                           gpu_results: dict[str, FixtureRunResult]) -> None:
    """Prove the frozen budget above is defect-sensitive, not vacuous.

    Two independent negative controls, both against already-computed
    in-memory trajectories (no tracked file is ever written or mutated; the
    next run re-parses the real, unmodified production output from disk):

    1. The same freeze/drop-step/perturb-one-component mutation battery
       RCG-04D introduced, applied to one representative fixture's GPU
       trajectory before comparing it to the CPU reference.
    2. A genuine physics-defect stand-in that costs no source mutation at
       all: the accepted ``+q`` and ``-q`` chiral partner fixtures (already
       tracked, already used by RCG-04H to assert opposite chirality) are
       byte-for-byte different physical states, not numerical noise around
       the same state. Comparing GPU ``bs1_plus`` against CPU ``bs1_minus``
       under the *same* frozen backend-parity budget must fail hugely, or
       the budget would be loose enough to wave through a real
       sign/handedness defect.
    """
    # ``moving_wall_feature_off`` is the representative fixture for the
    # mutation battery: it has the largest own physical displacement of any
    # "ordinary"-class fixture (0.733 rad, see the evidence document's
    # displacement table), so its tight ordinary-class budget (~1e-3 rad)
    # is guaranteed to be far smaller than a freeze/perturb defect's effect
    # regardless of which fixture the budget was calibrated against.
    # (``moving_all_fine_wide`` was tried first and rejected here: its own
    # displacement is only 0.0074 rad, so a freeze mutation on it produced
    # an error smaller than even the tight ordinary budget -- exactly the
    # vacuous-negative-control failure mode this rule exists to catch.)
    rep_name = "moving_wall_feature_off"
    rep = cpu_results[rep_name]
    gpu_rep = gpu_results[rep_name]
    budget = budget_for(precision, rep.fixture)

    frozen = dataclasses.replace(gpu_rep.trajectory, steps=tuple(
        dataclasses.replace(step, records=gpu_rep.trajectory.steps[0].records)
        for step in gpu_rep.trajectory.steps
    ))
    frozen_error = te.angular_trajectory_error(rep.trajectory, frozen)
    if frozen_error.max_radians <= budget["max_angular_rad"]:
        raise NegativeControlDidNotFailError(
            f"{precision} freeze control did not fail: {frozen_error.max_radians} <= {budget['max_angular_rad']}"
        )
    print(f"{precision} negative control [freeze]: max_radians={frozen_error.max_radians:.6g} "
          f"> budget {budget['max_angular_rad']:.3g} -- failed as expected")

    dropped = dataclasses.replace(
        gpu_rep.trajectory, steps=gpu_rep.trajectory.steps[:1] + gpu_rep.trajectory.steps[2:]
    )
    try:
        te.angular_trajectory_error(rep.trajectory, dropped)
    except te.TrajectoryKeyMismatchError:
        print(f"{precision} negative control [drop-step]: TrajectoryKeyMismatchError raised as expected")
    else:
        raise NegativeControlDidNotFailError(f"{precision} drop-step control did not fail")

    final = gpu_rep.trajectory.steps[-1]
    key = next(iter(final.records))
    record = final.records[key]
    corrupted = dataclasses.replace(record, direction=(record.direction[0] + 0.5, record.direction[1], record.direction[2]))
    perturbed_records = dict(final.records)
    perturbed_records[key] = corrupted
    perturbed = dataclasses.replace(
        gpu_rep.trajectory,
        steps=gpu_rep.trajectory.steps[:-1] + (dataclasses.replace(final, records=perturbed_records),),
    )
    perturbed_error = te.angular_trajectory_error(rep.trajectory, perturbed)
    if perturbed_error.max_radians <= budget["max_angular_rad"]:
        raise NegativeControlDidNotFailError(
            f"{precision} perturb-one-component control did not fail: "
            f"{perturbed_error.max_radians} <= {budget['max_angular_rad']}"
        )
    print(f"{precision} negative control [perturb-one-component]: max_radians={perturbed_error.max_radians:.6g} "
          f"> budget {budget['max_angular_rad']:.3g} -- failed as expected")

    plus_gpu = gpu_results["moving_dmi_chiral_bs1_plus"].trajectory
    minus_cpu = cpu_results["moving_dmi_chiral_bs1_minus"].trajectory
    sign_error = te.angular_trajectory_error(minus_cpu, plus_gpu)
    if sign_error.max_radians <= budget["max_angular_rad"]:
        raise NegativeControlDidNotFailError(
            f"{precision} +q-vs--q chirality sign control did not fail: "
            f"{sign_error.max_radians} <= {budget['max_angular_rad']} -- the frozen budget would "
            "silently accept a genuine DMI handedness defect"
        )
    print(f"{precision} negative control [+q vs -q chirality]: max_radians={sign_error.max_radians:.6g} "
          f"> budget {budget['max_angular_rad']:.3g} -- failed as expected "
          "(the frozen budget correctly rejects a real physics-sign defect)")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cuda-fp64-binary", type=Path)
    parser.add_argument("--cuda-fp32-binary", type=Path)
    parser.add_argument("--hip-fp64-binary", type=Path)
    parser.add_argument("--hip-fp32-binary", type=Path)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--mode", choices=("explore", "accept"), default="accept")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    args.workspace_root.mkdir(parents=True, exist_ok=True)

    cpu_workspace = prepare_workspace(args.workspace_root, "cpu", gpu=False)
    print(f"CPU fp64 workspace: {cpu_workspace}")
    cpu_results = run_backend(args.cpu_binary.resolve(), cpu_workspace, normalization_tolerance=1.0e-6)
    print(f"CPU fp64: {len(cpu_results)}/{len(FIXTURES)} fixtures ran successfully")

    all_json: dict = {"cpu_binary": str(args.cpu_binary), "backends": {}}
    backends_run: list[str] = []
    backends_deferred: list[str] = []

    def do_backend(label: str, binary: Path | None) -> list[FixtureComparison] | None:
        if binary is None:
            backends_deferred.append(label)
            return None
        workspace = prepare_workspace(args.workspace_root, label, gpu=True)
        print(f"\n{label} workspace: {workspace}")
        precision = "fp32" if label.endswith("fp32") else "fp64"
        normalization_tolerance = 1.0e-4 if precision == "fp32" else 1.0e-6
        results = run_backend(binary.resolve(), workspace, normalization_tolerance=normalization_tolerance)
        print(f"{label}: {len(results)}/{len(FIXTURES)} fixtures ran successfully")
        comparisons = [compare_fixture(cpu_results[f.name], results[f.name]) for f in FIXTURES]
        print_comparison_table(label, comparisons)
        all_json["backends"][label] = {"comparisons": [c.as_dict() for c in comparisons]}
        if args.mode == "accept":
            assert_budget(precision, comparisons)
            run_negative_controls(precision, cpu_results, results)
        backends_run.append(label)
        return comparisons

    do_backend("cuda_fp64", args.cuda_fp64_binary)
    do_backend("cuda_fp32", args.cuda_fp32_binary)
    do_backend("hip_fp64", args.hip_fp64_binary)
    do_backend("hip_fp32", args.hip_fp32_binary)

    print(f"\nBackends run: {backends_run}")
    print(f"Backends deferred (no binary supplied): {backends_deferred}")

    if args.json_out:
        args.json_out.write_text(json.dumps(all_json, indent=2))
        print(f"Machine-readable results written to {args.json_out}")

    if args.mode == "accept" and not backends_run:
        raise BackendParityError("no GPU backend binary was supplied; nothing to accept")

    print("RCG-04I backend-parity harness completed" + (" (explore mode, no assertions)" if args.mode == "explore" else ""))


if __name__ == "__main__":
    main()
