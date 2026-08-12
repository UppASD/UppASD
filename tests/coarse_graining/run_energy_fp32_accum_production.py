#!/usr/bin/env python3
"""RCG-06B ENERGY-FP32-ACCUM, production-scale layer (F-11).

The standalone ``adaptive-cg-energy-fp32-accum`` microbenchmark
(``test_energy_fp32_accum.cpp``) measures the accumulator-precision defect
class in isolation, over a wide sweep of accumulated-term counts. This
script complements it with the real production binary, reusing RCG-04I's
already-accepted backend-parity infrastructure
(``run_moving_backend_parity.py``) rather than inventing a new fixture set
or workspace/parsing mechanism -- the same 19 tracked ``moving_*`` fixtures,
the same CPU/CUDA workspace preparation, and the same
``trajectory_evidence.parse_energy_field_series``/``compare_energy_field_series``
term-by-term energy parser RCG-04I already fixed to understand both
backends' ``total=``/``last_total_energy_j=`` spellings (RCG-04I section
17.5).

Purpose (checklist items this fixture is direct evidence for):

- "Energy error scaling is measured over increasing system size": this
  tracked set already spans two system sizes (48 atoms:
  ``moving_all_fine``; 192 atoms: every other AdaptiveCG fixture), the
  widest span available without inventing new fixtures. Energy relative
  error is reported per fixture, grouped by ``natom``.
- "FP32 field parity and FP64 energy budgets are distinct": measured here,
  not assumed -- and the measurement itself is the finding worth recording.
  At these two production scales, CUDA fp64 vs CPU fp64 and CUDA fp32 vs
  CPU fp64 total-energy relative error are *nearly identical*
  (observed worst case: fp64 8.588e-05, fp32 8.592e-05), because both are
  dominated by the same already-documented CPU/GPU floating-point-
  associativity (summation-order) floor RCG-04I's own trajectory budgets
  are built on, not by accumulator precision -- reproducing RCG-04I section
  17.8's "fp32 error tracks fp64 error almost exactly" signature for this
  fixture class. This script therefore freezes one production-scale energy
  budget from the observed data (RCG-04I section 17.9 methodology) rather
  than two artificially separated ones. The genuinely distinct,
  order-of-magnitude-separated FP64-accumulator vs FP32-accumulator budgets
  are demonstrated by the accumulator-isolated standalone microbenchmark
  instead, where summation order is not a confound.

This script only exercises the ``total`` energy term (the one every
AdaptiveCG fixture's diagnostic emits and the one the fixed accumulator
sums into ``energyTerms_[7]``); per-term breakdowns are available in the
underlying ``FixtureComparison.energy_terms`` dict for anyone who wants
them but are not separately budgeted here.
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

from run_moving_backend_parity import (
    FIXTURES,
    FixtureComparison,
    compare_fixture,
    prepare_workspace,
    run_backend,
)


class EnergyBudgetError(AssertionError):
    pass


def energy_rel_errors(comparisons: list[FixtureComparison]) -> dict[int, list[float]]:
    """Map natom -> list of per-fixture 'total' energy relative errors."""
    by_size: dict[int, list[float]] = {}
    for comparison in comparisons:
        if comparison.energy_terms is None:
            continue  # feature-off fixtures have no AdaptiveCG energy diagnostic
        for sample in comparison.energy_terms["samples"]:
            rel_error = sample["energies_j"]["total"]["rel_error"]
            by_size.setdefault(comparison.fixture.natom, []).append(rel_error)
    return by_size


def print_size_scaling(label: str, by_size: dict[int, list[float]]) -> None:
    print(f"\n--- {label}: 'total' energy relative error by system size ---")
    for natom in sorted(by_size):
        values = by_size[natom]
        print(
            f"  natom={natom:4d}  n_fixtures={len(values):2d}  "
            f"mean={statistics.mean(values):.3e}  max={max(values):.3e}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--cuda-fp64-binary", type=Path, required=True)
    parser.add_argument("--cuda-fp32-binary", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    args.workspace_root.mkdir(parents=True, exist_ok=True)

    cpu_workspace = prepare_workspace(args.workspace_root, "cpu", gpu=False)
    cpu_results = run_backend(args.cpu_binary.resolve(), cpu_workspace, normalization_tolerance=1.0e-6)
    print(f"CPU fp64: {len(cpu_results)}/{len(FIXTURES)} fixtures ran successfully")

    fp64_workspace = prepare_workspace(args.workspace_root, "cuda_fp64", gpu=True)
    fp64_results = run_backend(args.cuda_fp64_binary.resolve(), fp64_workspace, normalization_tolerance=1.0e-6)
    fp64_comparisons = [compare_fixture(cpu_results[f.name], fp64_results[f.name]) for f in FIXTURES]

    fp32_workspace = prepare_workspace(args.workspace_root, "cuda_fp32", gpu=True)
    fp32_results = run_backend(args.cuda_fp32_binary.resolve(), fp32_workspace, normalization_tolerance=1.0e-4)
    fp32_comparisons = [compare_fixture(cpu_results[f.name], fp32_results[f.name]) for f in FIXTURES]

    fp64_by_size = energy_rel_errors(fp64_comparisons)
    fp32_by_size = energy_rel_errors(fp32_comparisons)
    print_size_scaling("CUDA fp64", fp64_by_size)
    print_size_scaling("CUDA fp32", fp32_by_size)

    fp64_all = [v for values in fp64_by_size.values() for v in values]
    fp32_all = [v for values in fp32_by_size.values() for v in values]
    fp64_worst = max(fp64_all)
    fp32_worst = max(fp32_all)

    # Finding, not an assumption this script started with: at the two
    # production scales this tracked fixture set offers (48/192 atoms), the
    # 'total' energy comparison against the CPU reference is dominated by
    # the same already-documented CPU/GPU floating-point-associativity
    # (summation-order) floor RCG-04I's own trajectory-error budgets are
    # built on (docs/RCG-04_MOVING_E2E_EVIDENCE.md section 17.7b/17.8) --
    # NOT by F-11's accumulator-precision defect, which this script's first
    # run (see the raw evidence in docs/RCG-06_MEMORY_TIMING_PRECISION_EVIDENCE.md)
    # showed producing an fp32 worst case (8.592e-05) essentially identical
    # to the fp64 worst case (8.588e-05) -- the same "fp32 error tracks fp64
    # error almost exactly" signature RCG-04I section 17.8 already
    # documented for this class of fixture, now re-confirmed after this
    # slice's fix. This is expected, not a gap: it means production-scale
    # 'total' energy comparison at these two sizes cannot, by itself, tell
    # accumulator precision apart from algorithmic summation order. The
    # accumulator-isolated ENERGY-FP32-ACCUM microbenchmark
    # (test_energy_fp32_accum.cpp) is what demonstrates the genuinely
    # distinct, order-of-magnitude-separated FP64-accumulator vs
    # FP32-accumulator budgets (double stays within 1e-9 up to N=3e6; float
    # grows to ~3e-2) -- this script's role is the complementary one of
    # confirming the fix does not regress production-scale accuracy, using
    # one budget derived from what is actually observed here (RCG-04I
    # section 17.9 methodology: derive from data, freeze with headroom, do
    # not assume near-machine-epsilon).
    PRODUCTION_ENERGY_BUDGET = 1.0e-3  # ~12x headroom over the observed 8.6e-5 worst case;
    # deliberately matches RCG-04I's own frozen "ordinary" fp64/fp32
    # trajectory budget (FROZEN_BUDGET_FP64_ORDINARY/FROZEN_BUDGET_FP32_ORDINARY
    # in run_moving_backend_parity.py), since both are budgeting the same
    # underlying CPU/GPU summation-order phenomenon.

    print(f"\nProduction 'total' energy budget: {PRODUCTION_ENERGY_BUDGET:.3e} "
          f"(observed worst fp64={fp64_worst:.3e}, fp32={fp32_worst:.3e})")

    failures = []
    if fp64_worst > PRODUCTION_ENERGY_BUDGET:
        failures.append(f"CUDA fp64 worst total-energy rel_error {fp64_worst:.3e} exceeds budget {PRODUCTION_ENERGY_BUDGET:.3e}")
    if fp32_worst > PRODUCTION_ENERGY_BUDGET:
        failures.append(f"CUDA fp32 worst total-energy rel_error {fp32_worst:.3e} exceeds budget {PRODUCTION_ENERGY_BUDGET:.3e}")

    if args.json_out:
        args.json_out.write_text(json.dumps({
            "fp64_by_size": fp64_by_size, "fp32_by_size": fp32_by_size,
            "production_energy_budget": PRODUCTION_ENERGY_BUDGET,
        }, indent=2))

    if failures:
        raise EnergyBudgetError("; ".join(failures))

    print("\nENERGY-FP32-ACCUM production layer: PASS -- production-scale energy "
          "comparison holds at both precisions; the distinct FP64-accumulator "
          "vs FP32-accumulator budgets are demonstrated by the accumulator-"
          "isolated microbenchmark, not by this whole-run comparison, which "
          "is dominated by CPU/GPU summation order at these sizes")


if __name__ == "__main__":
    main()
