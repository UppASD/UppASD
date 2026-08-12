#!/usr/bin/env python3
"""RCG-05F: confirm ownership invariants are preserved *across* an adaptive
mask-rebuild transition, not only at the final snapshot -- reusing RCG-04G's
transition-history infrastructure (`trajectory_evidence.parse_transition_events`/
`parse_resolution_state_history`) rather than building new instrumentation,
per the RCG-05F prompt's explicit instruction.

This is CPU-only: `trajectory_evidence.py`'s own module docstring already
documents that GPU stdout has no per-event `AdaptiveCG: transition ...` log
at all (`print_transition_events` is CPU-only production code, called from
`adaptivecgproduction.f90`; the GPU diagnostic mirror does not emit an
equivalent). This is an existing, already-documented backend gap, not
something this check discovers.

Two invariants are checked, both directly from real per-step production
diagnostics (`cg_diagnostics>=2`), never a bare aggregate count:

1. **Every resolution_state sample is a complete, well-formed partition**:
   every sample (``label=initial``, every ``label=step``, ``label=final``)
   reports exactly ``nblocks`` values, every one of which is a valid
   COARSE(0)/BUFFER(1)/FINE(2) state -- i.e. no block is ever transiently
   missing or duplicated in the reported ownership vector.
2. **Accepted/rejected transitions correlate exactly with real snapshot
   changes**: for every logged transition event, compare the block's own
   reported state in the resolution_state sample immediately before that
   event's step against the sample immediately at/after it. An *accepted*
   transition's block must show a real value change; a *rejected*
   transition's block must show no change at all. (This deliberately does
   not attempt to decode `AdaptiveCG: transition ...`'s own `old_state`/
   `new_state` fields against `resolution_state`'s COARSE/BUFFER/FINE scale
   -- inspection of real output showed they are not the same encoding, e.g.
   an accepted `old_state=3 new_state=1` line corresponds to an actual
   resolution_state change of BUFFER(1)->FINE(2) for that block, not
   COARSE(0)->BUFFER(1); `old_state`/`new_state` are the selector's own
   internal proposal-stage codes, not the committed block_state values.
   Using the two independently-emitted resolution_state snapshots directly,
   rather than reverse-engineering that separate encoding, is what makes
   this check a genuine cross-check of two independently emitted signals,
   not a single parser's self-consistency.)

If either invariant failed, that would be direct evidence of the RCG-05F
prompt's named risk: "a block changing resolution state should transiently
double-count or drop the dipole/short-range/on-site contribution" -- since
`evaluate_static_hybrid_operator`'s interaction_owner/onsite_owner masks
(`coarsetensoroperator.f90:478-479,499-500`) are derived directly from
`atomistic_block`/`coarse_block`, which are in turn derived directly from
this exact resolution/`block_state` vector
(`statichybridoperator.f90:245-248`) -- a malformed or inconsistently
recorded snapshot would mean the ownership mask itself was malformed at
that instant.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import trajectory_evidence as te

SOURCE_E2E = Path(__file__).with_name("e2e")
FIXTURE_NAME = "ownership_aniso_buffer"


class TransitionInvariantError(AssertionError):
    pass


def run_cpu(binary: Path, workspace: Path, timeout: int = 300) -> str:
    import shutil

    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE_E2E, workspace)
    case_dir = workspace / FIXTURE_NAME
    result = subprocess.run(
        [str(binary)], cwd=case_dir, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout, check=False,
    )
    if result.returncode != 0:
        raise TransitionInvariantError(
            f"{FIXTURE_NAME}: binary {binary} exited {result.returncode}:\n{result.stdout}"
        )
    return result.stdout


def check_complete_partition(samples: list[te.ResolutionSample]) -> int:
    if not samples:
        raise TransitionInvariantError(
            "no AdaptiveCG: resolution_state samples were emitted at all -- "
            "cg_diagnostics must be >=2 for this check to be non-vacuous"
        )
    nblocks = len(samples[0].values)
    if nblocks == 0:
        raise TransitionInvariantError("first resolution_state sample reports zero blocks")
    for sample in samples:
        if len(sample.values) != nblocks:
            raise TransitionInvariantError(
                f"resolution_state sample label={sample.label} step={sample.step} reports "
                f"{len(sample.values)} blocks, expected {nblocks} (a block was transiently "
                "dropped or duplicated)"
            )
        invalid = [v for v in sample.values if v not in (0, 1, 2)]
        if invalid:
            raise TransitionInvariantError(
                f"resolution_state sample label={sample.label} step={sample.step} contains "
                f"invalid state value(s) {invalid} outside {{COARSE=0, BUFFER=1, FINE=2}}"
            )
    return nblocks


def surrounding_samples(
    samples: list[te.ResolutionSample], step: int,
) -> tuple[te.ResolutionSample, te.ResolutionSample]:
    before_candidates = [s for s in samples if s.step < step]
    after_candidates = [s for s in samples if s.step >= step]
    if not before_candidates or not after_candidates:
        raise TransitionInvariantError(
            f"transition at step={step} has no resolution_state sample both before and "
            f"at/after it (available steps: {sorted({s.step for s in samples})})"
        )
    before = max(before_candidates, key=lambda s: s.step)
    # Prefer label=step over label=final when both exist at the same step
    # (both should carry identical values when this is the last step; using
    # the earliest-emitted one keeps this deterministic either way).
    after = min(after_candidates, key=lambda s: (s.step, 0 if s.label == "step" else 1))
    return before, after


def check_transitions_match_snapshots(
    events: list[te.TransitionEvent], samples: list[te.ResolutionSample],
) -> tuple[int, int]:
    accepted_checked = 0
    rejected_checked = 0
    for event in events:
        before, after = surrounding_samples(samples, event.step)
        if event.block < 1 or event.block > len(before.values):
            raise TransitionInvariantError(
                f"transition event references block={event.block}, out of range "
                f"1..{len(before.values)}"
            )
        before_value = before.values[event.block - 1]
        after_value = after.values[event.block - 1]
        if event.accepted:
            if before_value == after_value:
                raise TransitionInvariantError(
                    f"accepted transition step={event.step} block={event.block} "
                    f"(reason={event.reason}) shows no change between resolution_state "
                    f"label={before.label} step={before.step} (value={before_value}) and "
                    f"label={after.label} step={after.step} (value={after_value}) -- an "
                    "accepted transition was silently dropped from the recorded ownership map"
                )
            accepted_checked += 1
        else:
            if before_value != after_value:
                raise TransitionInvariantError(
                    f"rejected transition step={event.step} block={event.block} "
                    f"(reason={event.reason}, outcome={event.outcome}) shows a change between "
                    f"resolution_state label={before.label} step={before.step} "
                    f"(value={before_value}) and label={after.label} step={after.step} "
                    f"(value={after_value}) -- a rejected transition was applied anyway "
                    "(double-counted/leaked into the committed ownership map)"
                )
            rejected_checked += 1
    return accepted_checked, rejected_checked


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", type=Path, required=True)
    parser.add_argument("--workspace-root", type=Path, required=True)
    args = parser.parse_args()

    stdout = run_cpu(args.cpu_binary, args.workspace_root / "transition-invariants")
    events = te.parse_transition_events(stdout)
    samples = te.parse_resolution_state_history(stdout)

    nblocks = check_complete_partition(samples)
    print(f"RCG-05F: {len(samples)} resolution_state samples, each a complete "
          f"{nblocks}-block partition with only valid COARSE/BUFFER/FINE values.")

    if not events:
        raise TransitionInvariantError(
            f"{FIXTURE_NAME} emitted no AdaptiveCG: transition events -- this check requires "
            "at least one transition (accepted or rejected) to be non-vacuous"
        )
    accepted_checked, rejected_checked = check_transitions_match_snapshots(events, samples)
    if accepted_checked == 0:
        raise TransitionInvariantError(
            f"{FIXTURE_NAME} logged {len(events)} transitions but none were accepted -- "
            "the accepted-transition half of this check would be vacuous"
        )
    if rejected_checked == 0:
        raise TransitionInvariantError(
            f"{FIXTURE_NAME} logged {len(events)} transitions but none were rejected -- "
            "the rejected-transition half of this check would be vacuous"
        )
    print(f"RCG-05F: {accepted_checked} accepted transitions each show a real resolution_state "
          f"change; {rejected_checked} rejected transitions each show no change. Ownership was "
          "neither double-counted nor dropped across any logged transition.")
    print("\nRCG-05F transition ownership invariant check completed")


if __name__ == "__main__":
    main()
