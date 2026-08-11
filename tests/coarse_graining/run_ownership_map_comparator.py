#!/usr/bin/env python3
"""RCG-05C: run the CPU/GPU ownership-map comparator through the real
executable.

Two things are demonstrated, both through the real `sd.f95`/`sd.f95.cuda`(/
`sd.f95.hip`) binaries, never by reading production source alone:

1. **Regression** -- every existing RCG-04 `cg_mask_mode ADAPTIVE`/`STATIC`
   moving fixture (`run_moving_backend_parity.FIXTURES`, filtered to
   `adaptive_cg=True`; those cubic-ish geometries, `block_size_y=
   block_size_z=2` everywhere) must produce *identical* CPU/GPU final
   ownership maps, block for block -- reusing
   `run_moving_backend_parity.prepare_workspace`'s isolated-workspace/
   hash-verification pattern (RCG-04I), not reinventing it.
2. **Defect demonstration** -- the new, tracked, git-committed
   `e2e/ownership_aniso_buffer` fixture (RCG-05C; genuinely anisotropic
   `buffer_width_blocks=(2,1,1)`, `cg_mask_mode ADAPTIVE` so the GPU dilation
   kernel actually runs) must show CPU matching the correct (per-axis)
   dilation oracle and GPU matching the isotropic (scalar-collapsed) one --
   the same self-consistency mechanism `ownership_map_comparator.py`
   documents, run here against the real executable's output rather than a
   captured excerpt.

Per the RCG-05 prompt pack's evidence policy (SS3.3): every comparison below
prints the *actual* compared maps and, on any mismatch, the differing block
ids and their two states -- never only a match/mismatch count.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import shutil
import subprocess
from pathlib import Path

import ownership_map_comparator as omc
import static_topology_oracle as sto
from moving_state_generator import Geometry
from run_moving_backend_parity import FIXTURES as RCG04I_FIXTURES
from run_moving_backend_parity import GPU_ENABLE_BLOCK, sha256_file
from torque_oracle import parse_jfile_shells

SOURCE_E2E = Path(__file__).with_name("e2e")

ANISO_FIXTURE_NAME = "ownership_aniso_buffer"
REGRESSION_FIXTURE_NAMES: tuple[str, ...] = tuple(
    sorted({f.name for f in RCG04I_FIXTURES if f.adaptive_cg})
)
# The only RCG-04 fixture that engages the genuine GPU adaptive-transition
# runtime (`GpuAdaptiveRuntime::proposeSelectorState`/`dilateAdaptiveState`);
# every other regression fixture is `cg_mask_mode STATIC`, where GPU never
# re-dilates at all (it only copies the CPU-computed initial state once), so
# an exact match there is a necessary but much weaker check. See
# `docs/RCG-04_MOVING_E2E_EVIDENCE.md` SS13 / `run_moving_backend_parity.py`'s
# own module docstring for how this was established.
DILATION_ENGAGING_REGRESSION_FIXTURE = "moving_wall_adaptive"

# ownership_aniso_buffer's own geometry/jfile parameters (README.md /
# GENERATOR_MANIFEST.json), needed for the self-consistency oracle
# cross-check -- not re-derived from the fixture files at runtime, to keep
# this an independent check against what the fixture *claims* to be.
ANISO_GEOMETRY = Geometry(na=2, n1=6, n2=10, n3=9, basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
ANISO_BLOCK_SIZE = dict(block_size_x=1, block_size_y=2, block_size_z=3)


class OwnershipComparatorError(AssertionError):
    pass


def prepare_workspace(root: Path, backend: str, *, gpu: bool, fixtures: tuple[str, ...]) -> Path:
    """``run_moving_backend_parity.prepare_workspace``'s exact pattern
    (RCG-04I), generalized to an arbitrary fixture-name list so it also
    covers RCG-05C's own new fixture, which is not in RCG-04I's own
    hardcoded ``FIXTURES`` tuple."""
    workspace = root / backend / "e2e"
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE_E2E, workspace)

    for name in fixtures:
        src_inpsd = SOURCE_E2E / name / "inpsd.dat"
        dst_inpsd = workspace / name / "inpsd.dat"
        assert dst_inpsd.read_bytes() == src_inpsd.read_bytes(), (
            f"{backend}/{name}: inpsd.dat copy diverged from source before any GPU-enable block was appended"
        )
        if gpu:
            with dst_inpsd.open("a") as handle:
                handle.write(GPU_ENABLE_BLOCK)
    return workspace


def run_fixture(binary: Path, workspace: Path, name: str, *, timeout: int = 300) -> str:
    case_dir = workspace / name
    result = subprocess.run(
        [str(binary)], cwd=case_dir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, check=False,
    )
    if result.returncode != 0:
        raise OwnershipComparatorError(
            f"{name}: binary {binary} exited {result.returncode}:\n{result.stdout}"
        )
    return result.stdout


@dataclasses.dataclass
class RegressionOutcome:
    name: str
    comparison: omc.MapComparison


def run_regression_set(cpu_binary: Path, gpu_binary: Path, workspace_root: Path, label: str) -> list[RegressionOutcome]:
    fixtures = REGRESSION_FIXTURE_NAMES
    cpu_ws = prepare_workspace(workspace_root, f"{label}-cpu", gpu=False, fixtures=fixtures)
    gpu_ws = prepare_workspace(workspace_root, f"{label}-gpu", gpu=True, fixtures=fixtures)
    outcomes = []
    for name in fixtures:
        cpu_stdout = run_fixture(cpu_binary, cpu_ws, name)
        gpu_stdout = run_fixture(gpu_binary, gpu_ws, name)
        cpu_map = omc.parse_final_ownership_map(cpu_stdout, "CPU")
        gpu_map = omc.parse_final_ownership_map(gpu_stdout, label)
        comparison = omc.compare_ownership_maps(cpu_map, gpu_map, label=f"regression:{name}")
        outcomes.append(RegressionOutcome(name=name, comparison=comparison))
    return outcomes


def print_regression_table(label: str, outcomes: list[RegressionOutcome]) -> None:
    print(f"\n--- RCG-05C regression set ({label}) ---")
    for outcome in outcomes:
        c = outcome.comparison
        status = "MATCH" if c.matches else f"MISMATCH ({len(c.mismatched_block_ids)} blocks)"
        print(f"  {outcome.name:38s} nblocks={c.a.nblocks:4d} {status}")
        if not c.matches:
            print(f"      mismatched block ids: {list(c.mismatched_block_ids)}")
            print(f"      mismatch detail (block: cpu,{label}): {c.mismatch_detail}")


def assert_regression_matches_exactly(outcomes: list[RegressionOutcome]) -> None:
    failures = [o.name for o in outcomes if not o.comparison.matches]
    if failures:
        raise OwnershipComparatorError(
            f"regression fixtures disagreed between CPU and GPU (expected exact match on "
            f"cubic-ish RCG-04 geometries): {failures}"
        )


@dataclasses.dataclass
class AnisoOutcome:
    cpu_map: omc.OwnershipMap
    gpu_map: omc.OwnershipMap
    direct: omc.MapComparison
    cpu_self_consistency: omc.SelfConsistencyResult
    gpu_self_consistency: omc.SelfConsistencyResult
    wrap_axes: tuple[bool, bool, bool]
    coverage: omc.BondCoverageResult


def run_aniso_fixture(cpu_binary: Path, gpu_binary: Path, workspace_root: Path, label: str) -> AnisoOutcome:
    fixtures = (ANISO_FIXTURE_NAME,)
    cpu_ws = prepare_workspace(workspace_root, f"{label}-cpu", gpu=False, fixtures=fixtures)
    gpu_ws = prepare_workspace(workspace_root, f"{label}-gpu", gpu=True, fixtures=fixtures)
    cpu_stdout = run_fixture(cpu_binary, cpu_ws, ANISO_FIXTURE_NAME)
    gpu_stdout = run_fixture(gpu_binary, gpu_ws, ANISO_FIXTURE_NAME)

    cpu_map = omc.parse_final_ownership_map(cpu_stdout, "CPU")
    gpu_map = omc.parse_final_ownership_map(gpu_stdout, label)
    direct = omc.compare_ownership_maps(cpu_map, gpu_map, label=f"aniso:{label}")

    jfile_text = (SOURCE_E2E / ANISO_FIXTURE_NAME / "jfile").read_text()
    shells = parse_jfile_shells(jfile_text)
    cpu_sc = omc.self_consistency_check(cpu_map, ANISO_GEOMETRY, shells, **ANISO_BLOCK_SIZE)
    gpu_sc = omc.self_consistency_check(gpu_map, ANISO_GEOMETRY, shells, **ANISO_BLOCK_SIZE)

    topology = sto.compute_expected_topology(
        ANISO_GEOMETRY, shells, fine_block_ids=cpu_map.fine_block_ids(), **ANISO_BLOCK_SIZE,
    )
    wrap_axes = omc.periodic_wrap_axes(topology)
    coverage = omc.bond_coverage(ANISO_GEOMETRY, shells, topology.atom_block_by_atom, cpu_map, gpu_map)

    return AnisoOutcome(
        cpu_map=cpu_map, gpu_map=gpu_map, direct=direct,
        cpu_self_consistency=cpu_sc, gpu_self_consistency=gpu_sc,
        wrap_axes=wrap_axes, coverage=coverage,
    )


def print_aniso_report(label: str, outcome: AnisoOutcome) -> None:
    print(f"\n--- RCG-05C defect demonstration: {ANISO_FIXTURE_NAME} (CPU vs {label}) ---")
    print(f"  CPU block counts:        {outcome.cpu_map.block_counts()}")
    print(f"  {label} block counts:{' ' * max(0, 8 - len(label))}{outcome.gpu_map.block_counts()}")
    print(f"  CPU fine seed blocks:    {sorted(outcome.cpu_map.fine_block_ids())}")
    print(f"  {label} fine seed blocks:{' ' * max(0, 1 - len(label))}{sorted(outcome.gpu_map.fine_block_ids())}")
    print(f"  direct CPU-vs-{label} match: {outcome.direct.matches} "
          f"({len(outcome.direct.mismatched_block_ids)} of {outcome.direct.a.nblocks} blocks differ)")
    print(f"  mismatched block ids: {list(outcome.direct.mismatched_block_ids)}")
    print(f"  mismatch detail (block: cpu,{label}): {outcome.direct.mismatch_detail}")
    print(f"  CPU  correct_buffer_width={outcome.cpu_self_consistency.correct_buffer_width} "
          f"matches_correct_oracle={outcome.cpu_self_consistency.matches_correct_oracle} "
          f"matches_isotropic_oracle={outcome.cpu_self_consistency.matches_isotropic_oracle}")
    print(f"  {label} isotropic_buffer_width={outcome.gpu_self_consistency.isotropic_buffer_width} "
          f"matches_correct_oracle={outcome.gpu_self_consistency.matches_correct_oracle} "
          f"matches_isotropic_oracle={outcome.gpu_self_consistency.matches_isotropic_oracle}")
    print(f"  periodic wrap axes exercised (x,y,z): {outcome.wrap_axes}")
    print(f"  cross-interface bond coverage: total={outcome.coverage.total_bonds} "
          f"cpu_interface_bonds={outcome.coverage.interface_bonds_a} "
          f"{label}_interface_bonds={outcome.coverage.interface_bonds_b} "
          f"disagreeing_endpoints={len(outcome.coverage.disagreeing_bond_endpoints)}")


def assert_aniso_outcome(outcome: AnisoOutcome, *, expect_match: bool) -> None:
    """``expect_match=False`` (the default, and this slice's own registered
    CTest expectation): the buffer-width scalarization defect is not yet
    fixed, so CPU and GPU must currently *disagree*, and the disagreement
    must be attributable specifically to that defect (CPU matches the
    correct per-axis oracle, GPU matches the isotropic one) -- RCG-05C's own
    negative-control evidence.

    ``expect_match=True`` (``--expect-aniso-match``, for RCG-05D to pass once
    its fix lands, re-running this same script and CTest registration
    unchanged rather than duplicating it): CPU and GPU must now agree
    exactly, restoring the ordinary regression-test meaning of "the
    ownership-map comparator passes" the prompt pack's RCG-05D section
    describes ("after RCG-05D's fix, it must pass").
    """
    if not outcome.cpu_self_consistency.matches_correct_oracle:
        raise OwnershipComparatorError(
            "CPU's own reported ownership map on the anisotropic fixture does not match the "
            "correct per-axis dilation oracle -- either the fixture or the CPU path itself is "
            "broken; this is not a GPU-only defect demonstration if CPU disagrees with its own "
            "documented algorithm"
        )
    if expect_match:
        if not outcome.direct.matches:
            raise OwnershipComparatorError(
                "CPU and GPU ownership maps still disagree on the anisotropic fixture under "
                "--expect-aniso-match: the buffer-width scalarization fix has not restored "
                f"exact CPU/GPU agreement ({len(outcome.direct.mismatched_block_ids)} blocks "
                f"differ): {outcome.direct.mismatch_detail}"
            )
        return
    if outcome.direct.matches:
        raise OwnershipComparatorError(
            "CPU and GPU ownership maps matched exactly on the anisotropic fixture: the "
            "buffer-width scalarization defect was NOT reproduced (either already fixed -- pass "
            "--expect-aniso-match once that is confirmed -- or this fixture no longer exercises "
            "it) -- this comparator's own pre-fix negative-control purpose requires this fixture "
            "to currently fail"
        )
    if not outcome.gpu_self_consistency.matches_isotropic_oracle:
        raise OwnershipComparatorError(
            "GPU's ownership map disagreed with CPU, but does not match the isotropic-dilation "
            "oracle either -- the disagreement is real but not attributable to the specific "
            "buffer-width scalarization this slice is chartered to demonstrate; investigate "
            "before treating this as RCG-05C's defect evidence"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cuda-fp64-binary", type=Path)
    parser.add_argument("--cuda-fp32-binary", type=Path)
    parser.add_argument("--hip-fp64-binary", type=Path)
    parser.add_argument("--hip-fp32-binary", type=Path)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument(
        "--mode", choices=("explore", "accept"), default="accept",
        help="accept (default): assert regression matches and check the aniso fixture against "
             "--expect-aniso-match. explore: print only, assert nothing.",
    )
    parser.add_argument(
        "--expect-aniso-match", action="store_true",
        help="RCG-05C default (omitted): expect CPU and GPU to still DISAGREE on "
             "ownership_aniso_buffer (the buffer-width scalarization defect is not yet fixed) "
             "and require the disagreement to match the isotropic-dilation oracle specifically. "
             "Pass this flag (RCG-05D, post-fix): expect CPU and GPU to now AGREE exactly.",
    )
    args = parser.parse_args()
    args.workspace_root.mkdir(parents=True, exist_ok=True)

    backends: list[tuple[str, Path]] = []
    if args.cuda_fp64_binary:
        backends.append(("CUDA-fp64", args.cuda_fp64_binary))
    if args.cuda_fp32_binary:
        backends.append(("CUDA-fp32", args.cuda_fp32_binary))
    if args.hip_fp64_binary:
        backends.append(("HIP-fp64", args.hip_fp64_binary))
    if args.hip_fp32_binary:
        backends.append(("HIP-fp32", args.hip_fp32_binary))
    if not backends:
        raise OwnershipComparatorError("no GPU backend binary was supplied; nothing to compare")

    all_json: dict = {"cpu_binary": str(args.cpu_binary), "backends": {}}
    for label, gpu_binary in backends:
        workspace_root = args.workspace_root / label
        regression = run_regression_set(args.cpu_binary, gpu_binary, workspace_root, label)
        print_regression_table(label, regression)

        aniso = run_aniso_fixture(args.cpu_binary, gpu_binary, workspace_root, label)
        print_aniso_report(label, aniso)

        aniso_record = omc.evidence_record(
            fixture_name=ANISO_FIXTURE_NAME, geometry=ANISO_GEOMETRY,
            block_size=(1, 2, 3), cell_vectors=sto._IDENTITY_CELL, cg_buffer_blocks=0,
            maps={"CPU": aniso.cpu_map, label: aniso.gpu_map},
            comparisons=[aniso.direct],
            self_consistency=[aniso.cpu_self_consistency, aniso.gpu_self_consistency],
            wrap_axes=aniso.wrap_axes, coverage=aniso.coverage,
        )
        all_json["backends"][label] = {
            "regression": [
                {"fixture": o.name, "comparison": o.comparison.as_dict()} for o in regression
            ],
            "aniso_fixture": aniso_record,
        }

        if args.mode == "accept":
            assert_regression_matches_exactly(regression)
            assert_aniso_outcome(aniso, expect_match=args.expect_aniso_match)
            dilation_engaging = next(
                (o for o in regression if o.name == DILATION_ENGAGING_REGRESSION_FIXTURE), None
            )
            if dilation_engaging is None:
                raise OwnershipComparatorError(
                    f"{DILATION_ENGAGING_REGRESSION_FIXTURE!r} (the only regression fixture that "
                    "engages the real GPU dilation kernel) is missing from the regression set"
                )
            if args.expect_aniso_match:
                print(f"\n{label}: regression set ({len(regression)} fixtures) matched exactly; "
                      f"{ANISO_FIXTURE_NAME} now also matches exactly under --expect-aniso-match "
                      "(RCG-05D's fix restored CPU/GPU agreement).")
            else:
                print(f"\n{label}: regression set ({len(regression)} fixtures, including the "
                      f"dilation-engaging {DILATION_ENGAGING_REGRESSION_FIXTURE!r}) matched exactly; "
                      f"{ANISO_FIXTURE_NAME}'s buffer-width scalarization defect was reproduced "
                      "(CPU matches the correct oracle, GPU matches the isotropic one).")

    if args.json_out:
        args.json_out.write_text(json.dumps(all_json, indent=2))
        print(f"\nMachine-readable evidence written to {args.json_out}")

    print("\nRCG-05C ownership-map comparator completed"
          + (" (explore mode, no assertions)" if args.mode == "explore" else ""))


if __name__ == "__main__":
    main()
