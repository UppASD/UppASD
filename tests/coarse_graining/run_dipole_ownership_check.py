#!/usr/bin/env python3
"""RCG-05F: CPU-vs-GPU cross-check that the regular-grid FFT dipole does not
perturb fine/buffer/coarse ownership, at RCG-05C's full per-block identity
resolution (never aggregate counts alone), across a genuinely anisotropic
geometry (RCG-05B's own generator), plus the existing RCG-04-era
`gpu_fft_static_mixed`/`parity_static_cpu` dipole fixture pair. A skew-cell
sibling was attempted and dropped -- see the `PAIRS` comment below for the
exact, unrelated `neighbourmap.f90` blocker; skew coverage for this claim
comes from the Fortran unit tests instead (also noted below).

Why this is a CPU-vs-GPU cross-check and not a CPU dipole computation: CPU
never computes a nonzero dipole field in production at all --
`setup_static_hybrid_operator` rejects `use_uniform_coarse_dipole`
outright (`source/CoarseGraining/statichybridoperator.f90:115-119`,
"FFT dipole is an independent all-grid owner and cannot be embedded in
short-range hybrid dispatch"), and `evaluate_static_hybrid_operator`'s own
argument list has no dipole parameter at all
(`source/CoarseGraining/adaptivecgproduction.f90:1002-1006` never passes
one). The regular-grid FFT dipole (`gpu_dipole_mode EWALD3D_FFT`) is a
GPU-exclusive production capability requiring `do_gpu=Y`/`do_gpu_llg=Y`
(`adaptivecgproduction.f90:766-773`). So "CPU-vs-GPU cross-check" here means:
does turning GPU's dipole on change GPU's own reported ownership map,
relative to the same geometry's CPU (necessarily dipole-off) run? -- not
"do CPU and GPU compute the same dipole energy" (CPU cannot compute one at
all).

Why `cg_mask_mode STATIC`, not ADAPTIVE: under STATIC, GPU never re-dilates
the mask at runtime -- it only ever copies CPU's setup-time `block_state`
once (`tests/coarse_graining/e2e/ownership_aniso_buffer/README.md`
documents this same fact for why *that* fixture had to use ADAPTIVE to
expose the buffer-width dilation defect at all). That makes STATIC mode the
one case where an *exact* CPU/GPU ownership-map match is the correct,
unconfounded expectation for a dipole cross-check: a nonzero dipole field
genuinely changes the LLG trajectory, which could legitimately (not as a
bug) change which blocks an ADAPTIVE selector later refines/coarsens, which
would make a literal map-equality assertion meaningless under ADAPTIVE.
STATIC's ownership is fixed at setup, before any field (dipole or
otherwise) is even evaluated, so dipole cannot reach it.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import shutil
import subprocess
from pathlib import Path

import ownership_map_comparator as omc
from run_moving_backend_parity import GPU_ENABLE_BLOCK

SOURCE_E2E = Path(__file__).with_name("e2e")

# Matches gpu_fft_static_mixed/inpsd.dat's already-validated, already-
# production-tested dipole settings exactly (CG-10.5's own GPU_FFT_CASE),
# reused rather than re-chosen, so this check is not silently exercising a
# different, unvalidated dipole configuration.
DIPOLE_ENABLE_BLOCK = (
    "\ngpu_dipole_mode EWALD3D_FFT\ngpu_dipole_surface TINFOIL\n"
    "gpu_dipole_tol 1.0d-10\ngpu_dipole_mesh 0 0 0\n"
)


class DipoleOwnershipCheckError(AssertionError):
    pass


@dataclasses.dataclass(frozen=True)
class DipolePair:
    label: str
    cpu_fixture: str  # dipole-free, run unmodified as the CPU reference
    gpu_fixture: str  # GPU comparand; may be the same directory as cpu_fixture
    gpu_needs_dipole_block: bool  # False when gpu_fixture already bakes dipole in (gpu_fft_static_mixed)
    gpu_needs_gpu_block: bool  # False when gpu_fixture already bakes do_gpu in


# The new RCG-05B-generator-built, genuinely anisotropic fixture (RCG-05F),
# plus the pre-existing RCG-04-era CG-10.5 dipole pair, upgraded here from an
# aggregate `final_state` count comparison (`run_production_e2e.py` line
# ~437) to RCG-05C's full per-block identity comparator.
#
# A genuinely skew-cell sibling (`ownership_dipole_skew`, same construction,
# non-orthogonal `cell_vectors`) was also generated and attempted here, but
# hit a pre-existing `neighbourmap.f90` limitation unrelated to buffer-width
# or dipole ownership -- `setup_nm: no basis match` for the same single
# same-sublattice bond that maps cleanly on every orthogonal fixture in this
# repository -- and was removed rather than committed broken. RCG-05F's
# skew coverage for the dipole-exactly-once/non-overlap claims instead comes
# from `test_coarse_tensor_operator.f90`'s
# `test_dipole_unmasked_and_exactly_once` and
# `test_static_hybrid_operator.f90`'s
# `test_anisotropic_skew_ownership_non_overlap`, which exercise a skew cell
# directly against the operator (bypassing `neighbourmap.f90` entirely, since
# unit tests build bonds/topology by hand). The full-executable skew
# cross-check remains an open item; see the RCG-05F evidence doc.
PAIRS: tuple[DipolePair, ...] = (
    DipolePair(
        label="ownership_dipole_unequal_width",
        cpu_fixture="ownership_dipole_unequal_width",
        gpu_fixture="ownership_dipole_unequal_width",
        gpu_needs_dipole_block=True, gpu_needs_gpu_block=True,
    ),
    DipolePair(
        label="gpu_fft_static_mixed (RCG-04-era, block_size 1/2/2)",
        cpu_fixture="parity_static_cpu",
        gpu_fixture="gpu_fft_static_mixed",
        gpu_needs_dipole_block=False, gpu_needs_gpu_block=False,
    ),
)


def prepare_workspace(root: Path, tag: str, *, gpu_block: bool, dipole_block: bool) -> Path:
    workspace = root / tag
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE_E2E, workspace)
    if gpu_block or dipole_block:
        for fixture_dir in workspace.iterdir():
            inpsd = fixture_dir / "inpsd.dat"
            if not inpsd.is_file():
                continue
            with inpsd.open("a") as handle:
                if gpu_block:
                    handle.write(GPU_ENABLE_BLOCK)
                if dipole_block:
                    handle.write(DIPOLE_ENABLE_BLOCK)
    return workspace


def run_fixture(binary: Path, workspace: Path, name: str, *, timeout: int = 300) -> str:
    case_dir = workspace / name
    result = subprocess.run(
        [str(binary)], cwd=case_dir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, check=False,
    )
    if result.returncode != 0:
        raise DipoleOwnershipCheckError(
            f"{name}: binary {binary} exited {result.returncode}:\n{result.stdout}"
        )
    return result.stdout


@dataclasses.dataclass
class PairOutcome:
    pair: DipolePair
    cpu_map: omc.OwnershipMap
    gpu_map: omc.OwnershipMap
    comparison: omc.MapComparison


def run_pair(cpu_binary: Path, gpu_binary: Path, workspace_root: Path, gpu_label: str, pair: DipolePair) -> PairOutcome:
    cpu_ws = prepare_workspace(workspace_root, f"{pair.label}-cpu", gpu_block=False, dipole_block=False)
    gpu_ws = prepare_workspace(
        workspace_root, f"{pair.label}-gpu",
        gpu_block=pair.gpu_needs_gpu_block, dipole_block=pair.gpu_needs_dipole_block,
    )
    cpu_stdout = run_fixture(cpu_binary, cpu_ws, pair.cpu_fixture)
    gpu_stdout = run_fixture(gpu_binary, gpu_ws, pair.gpu_fixture)
    cpu_map = omc.parse_final_ownership_map(cpu_stdout, "CPU")
    gpu_map = omc.parse_final_ownership_map(gpu_stdout, gpu_label)
    comparison = omc.compare_ownership_maps(cpu_map, gpu_map, label=f"dipole:{pair.label}")
    return PairOutcome(pair=pair, cpu_map=cpu_map, gpu_map=gpu_map, comparison=comparison)


def print_outcome(label: str, outcome: PairOutcome) -> None:
    c = outcome.comparison
    status = "MATCH" if c.matches else f"MISMATCH ({len(c.mismatched_block_ids)} blocks)"
    print(f"  {outcome.pair.label:55s} nblocks={c.a.nblocks:4d} {status}")
    print(f"      CPU (no dipole) block counts: {outcome.cpu_map.block_counts()}")
    print(f"      {label} (dipole on) block counts: {outcome.gpu_map.block_counts()}")
    if not c.matches:
        print(f"      mismatched block ids: {list(c.mismatched_block_ids)}")
        print(f"      mismatch detail (block: cpu,{label}): {c.mismatch_detail}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cuda-fp64-binary", type=Path)
    parser.add_argument("--cuda-fp32-binary", type=Path)
    parser.add_argument("--hip-fp64-binary", type=Path)
    parser.add_argument("--hip-fp32-binary", type=Path)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--mode", choices=("explore", "accept"), default="accept")
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
        raise DipoleOwnershipCheckError("no GPU backend binary was supplied; nothing to compare")

    all_json: dict = {"cpu_binary": str(args.cpu_binary), "backends": {}}
    any_mismatch = False
    for gpu_label, gpu_binary in backends:
        print(f"\n--- RCG-05F dipole ownership cross-check ({gpu_label}) ---")
        outcomes = []
        for pair in PAIRS:
            outcome = run_pair(args.cpu_binary, gpu_binary, args.workspace_root / gpu_label, gpu_label, pair)
            outcomes.append(outcome)
            print_outcome(gpu_label, outcome)
            any_mismatch = any_mismatch or not outcome.comparison.matches
        all_json["backends"][gpu_label] = [
            {
                "label": o.pair.label,
                "cpu_map": o.cpu_map.as_dict(),
                "gpu_map": o.gpu_map.as_dict(),
                "comparison": o.comparison.as_dict(),
            }
            for o in outcomes
        ]

    if args.json_out:
        args.json_out.write_text(json.dumps(all_json, indent=2))
        print(f"\nMachine-readable evidence written to {args.json_out}")

    if args.mode == "accept" and any_mismatch:
        raise DipoleOwnershipCheckError(
            "at least one CPU(no dipole)/GPU(dipole) ownership map pair disagreed under "
            "cg_mask_mode STATIC, where GPU never re-dilates and only copies CPU's setup-time "
            "block_state -- the FFT dipole leaked into ownership resolution somewhere"
        )
    print("\nRCG-05F dipole ownership cross-check completed"
          + (" (explore mode, no assertions)" if args.mode == "explore" else ""))


if __name__ == "__main__":
    main()
