#!/usr/bin/env python3
"""CGP-00 Part D: transition-decision stability evidence.

docs/CGP_work.md CGP-00 Part D asks whether changing the adaptive-CG energy
reduction's precision can flip a resolution-transition accept/reject
decision, and asks for the sweep to be centered on real observed jumps
rather than only tested far from the boundary.

Two things had to be established by reading the code (not assumed) before
this script could be written, and both are worth recording plainly:

1. The energy-jump accept/reject gate (``abs(energy_jump) <=
   energy_jump_limit_j``) is implemented exactly once in this codebase, in
   ``source/CoarseGraining/adaptivehybridsolver.f90``'s
   ``apply_adaptive_transitions``, which is the **CPU-only** adaptive-CG
   path. Its ``energy_before``/``energy_after``/``energy_jump`` are computed
   entirely in ``dblprec`` (Fortran double) -- there is no FP32 anywhere in
   that call chain, so this task's Part A precision sweep (SINGLE_PREC local
   contributions cast into an FP64 GPU reduction) does not apply to it
   directly.
2. The **GPU** adaptive-CG runtime (``source/gpu_files/gpuAdaptiveRuntime.cpp``,
   this phase's actual primary scope) has no energy-jump gate at all in its
   transition path: ``proposeSelectorState``/``publishProposedState`` decide
   and publish transitions from selector/polarization thresholds only, never
   evaluate or compare an energy jump. ``adaptive_energy_jump_limit_j`` is
   read on the GPU side (``gpuSimulation.cpp``) only to print it in a
   diagnostic line; nothing compares against it. This is independently
   confirmed by this repository's own RCG-04I docstring
   (``run_moving_backend_parity.py``): "the GPU AdaptiveCG diagnostic path
   ... prints only initial/final resolution-state snapshots, never a
   per-event transition log, and its one aggregate summary line hardcodes
   ``rejected_transitions=0``". A hardcoded zero-rejection count is exactly
   what "no gate exists" looks like from the outside.

So there is currently no live GPU decision for a precision change to flip.
This script instead does the next most useful thing: extract the real
``energy_jump_j`` magnitudes the CPU-only path actually produces for a
tracked, transition-heavy fixture (``moving_wall_adaptive``, a domain wall
swept through an adaptive block grid over 900 steps -- the only tracked
fixture whose adaptive transitions are driven by real physics rather than a
synthetic construction), and reports how they sit relative to the
configured ``cg_energy_jump_limit``. Combined with
``test_energy_hierarchical_precision.cpp``'s synthetic
``energy_difference(delta=...)`` case (CGP-00 Part B/A), which sweeps delta
down through 1e-9 and shows exactly how much each reduction-precision mode
would corrupt a jump of that size, this brackets the question: is the real
jump magnitude distribution inside a regime where the swept synthetic deltas
already show a precision-dependent risk? That comparison is drawn in the
CGP-00 evidence document, not computed here.
"""

from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

from run_moving_backend_parity import FIXTURES, prepare_workspace, run_fixture

TRANSITION_RE = re.compile(
    r"AdaptiveCG: transition step=(?P<step>\d+) block=(?P<block>\d+) "
    r"old_state=(?P<old>\d+) new_state=(?P<new>\d+) accepted=(?P<accepted>[TF]) "
    r"reason=(?P<reason>\S*) outcome=(?P<outcome>\S*) "
    r"energies_j=\s*(?P<before>[-+0-9.eE]+)\s+(?P<after>[-+0-9.eE]+)\s+(?P<jump>[-+0-9.eE]+)\s+"
    r"polarization_ratio=\s*(?P<ratio>[-+0-9.eE]+)"
)

REJECTED_COUNT_RE = re.compile(r"rejected_transitions=(\d+)")


def parse_energy_jump_limit(inpsd_path: Path) -> float:
    for line in inpsd_path.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == "cg_energy_jump_limit":
            return float(parts[1])
    raise AssertionError(f"{inpsd_path}: no cg_energy_jump_limit line found")


def parse_transition_events(stdout: str) -> list[dict]:
    events = []
    for match in TRANSITION_RE.finditer(stdout):
        events.append({
            "step": int(match.group("step")),
            "block": int(match.group("block")),
            "accepted": match.group("accepted") == "T",
            "reason": match.group("reason"),
            "outcome": match.group("outcome"),
            "energy_before_j": float(match.group("before")),
            "energy_after_j": float(match.group("after")),
            "energy_jump_j": float(match.group("jump")),
        })
    return events


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--fixture", default="moving_wall_adaptive")
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    fixture = next((f for f in FIXTURES if f.name == args.fixture), None)
    if fixture is None:
        raise SystemExit(f"unknown tracked fixture {args.fixture!r}; see FIXTURES in "
                         "run_moving_backend_parity.py")

    args.workspace_root.mkdir(parents=True, exist_ok=True)
    workspace = prepare_workspace(args.workspace_root, "cpu_transition_sweep", gpu=False)
    result = run_fixture(args.cpu_binary.resolve(), workspace, fixture)

    limit = parse_energy_jump_limit(workspace / fixture.name / "inpsd.dat")
    events = parse_transition_events(result.stdout)
    reject_matches = REJECTED_COUNT_RE.findall(result.stdout)

    print(f"CGP-00 Part D: {fixture.name} ({fixture.natom} atoms, {fixture.nstep} steps), "
          f"cg_energy_jump_limit={limit:.6e} J")
    print(f"  transition events parsed from stdout: {len(events)}")
    if reject_matches:
        print(f"  print_adaptive_cg_summary rejected_transitions field(s) seen: {reject_matches}")

    if not events:
        print("  No transition events were printed. Either the fixture produced no "
              "block-state transitions in this run, or cg_diagnostics is below 2 in "
              "its inpsd.dat -- both are evidence outcomes worth recording, not a "
              "script defect.")
        payload = {"fixture": fixture.name, "energy_jump_limit_j": limit, "events": []}
        if args.json_out:
            args.json_out.write_text(json.dumps(payload, indent=2))
        return

    jumps = [abs(event["energy_jump_j"]) for event in events]
    accepted = [event for event in events if event["accepted"]]
    rejected = [event for event in events if not event["accepted"]]
    energy_rejected = [event for event in rejected if event["outcome"] == "energy-jump-rejected"]

    print(f"  accepted={len(accepted)} rejected={len(rejected)} "
          f"(energy-jump-rejected={len(energy_rejected)})")
    print(f"  |energy_jump_j| over all events: min={min(jumps):.6e} "
          f"median={statistics.median(jumps):.6e} max={max(jumps):.6e}")

    # The quantity CGP-00 Part D actually cares about: how close each real
    # jump sits to the configured limit. A ratio near 1.0 is exactly the
    # boundary-adjacent regime the task asks to sweep, not avoid.
    ratios = [j / limit for j in jumps if limit > 0]
    if ratios:
        near_boundary = sum(1 for r in ratios if 0.1 <= r <= 10.0)
        print(f"  jump / limit ratio: min={min(ratios):.3e} median={statistics.median(ratios):.3e} "
              f"max={max(ratios):.3e}; {near_boundary}/{len(ratios)} events within 10x of the "
              "boundary in either direction")

    payload = {
        "fixture": fixture.name,
        "energy_jump_limit_j": limit,
        "n_events": len(events),
        "n_accepted": len(accepted),
        "n_rejected": len(rejected),
        "n_energy_jump_rejected": len(energy_rejected),
        "jump_abs_j": {"min": min(jumps), "median": statistics.median(jumps), "max": max(jumps)},
        "events": events,
    }
    if args.json_out:
        args.json_out.write_text(json.dumps(payload, indent=2))

    print("\nCGP-00 Part D data collection: DONE (this script only gathers evidence; "
          "the CGP-00 evidence document draws the precision-risk conclusion by "
          "comparing these real magnitudes against "
          "test_energy_hierarchical_precision.cpp's synthetic energy_difference sweep)")


if __name__ == "__main__":
    main()
