#!/usr/bin/env python3
"""RCG-06A MEM-LARGE-HOST: run the ordinary UppASD executable on a large
(Natom=8000) adaptive-CG host under a constrained process stack limit, and
confirm a minimal fault-injection reproduction of the eliminated
automatic-array pattern (F-13) is killed under the identical limit.

See tests/coarse_graining/e2e/mem_large_host/README.md for the full
provenance, stack-limit arithmetic, and negative-control rationale.
"""

from __future__ import annotations

import argparse
import resource
import subprocess
import sys
from pathlib import Path

STACK_LIMIT_BYTES = 512 * 1024  # 512 KiB; see README.md for the arithmetic
CASE = "mem_large_host"


def _limit_stack() -> None:
    resource.setrlimit(resource.RLIMIT_STACK, (STACK_LIMIT_BYTES, STACK_LIMIT_BYTES))


def run_production(binary: Path, case_dir: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary)],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
        check=False,
        preexec_fn=_limit_stack,
    )


def run_fault_injection(binary: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(binary)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
        check=False,
        preexec_fn=_limit_stack,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--fault-injection-binary", required=True, type=Path)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent
    case_dir = repo_root / "e2e" / CASE
    if not case_dir.is_dir():
        raise AssertionError(f"missing fixture directory {case_dir}")

    # Negative control first: the eliminated automatic-array pattern, at the
    # exact same Natom and stack limit, must fail.
    injected = run_fault_injection(args.fault_injection_binary)
    if injected.returncode == 0:
        raise AssertionError(
            "fault-injection reproduction of the pre-RCG-06A automatic "
            "arrays did NOT crash under ulimit -s "
            f"{STACK_LIMIT_BYTES // 1024} KiB -- the negative control is "
            f"not discriminating at Natom=8000. stdout:\n{injected.stdout}"
        )
    print(
        "fault-injection negative control: killed as expected "
        f"(returncode={injected.returncode})"
    )

    # Positive evidence: the real, fixed production binary must complete
    # cleanly at the same Natom under the same stack limit.
    result = run_production(args.binary, case_dir)
    if result.returncode != 0:
        raise AssertionError(
            "adaptive-CG production run at Natom=8000 failed under ulimit -s "
            f"{STACK_LIMIT_BYTES // 1024} KiB (returncode={result.returncode}). "
            f"stdout:\n{result.stdout}"
        )
    if "AdaptiveCG: capability accepted" not in result.stdout:
        raise AssertionError(
            f"adaptive CG setup was not accepted\n{result.stdout}"
        )
    if "AdaptiveCG: resolution_state label=final" not in result.stdout:
        raise AssertionError(
            f"missing final resolution_state diagnostic\n{result.stdout}"
        )
    print(
        "MEM-LARGE-HOST: production run completed under ulimit -s "
        f"{STACK_LIMIT_BYTES // 1024} KiB at Natom=8000"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
