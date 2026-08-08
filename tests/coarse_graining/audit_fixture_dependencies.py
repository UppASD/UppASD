#!/usr/bin/env python3
"""Check executable adaptive-CG fixture dependencies against ``git ls-files``."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

from fixture_dependencies import EXAMPLE_CASES, all_e2e_cases

REPO_ROOT = Path(__file__).resolve().parents[2]
E2E_ROOT = REPO_ROOT / "tests/coarse_graining/e2e"
EXAMPLES_ROOT = REPO_ROOT / "examples/AdaptiveCoarseGraining"
INPUT_FILE_KEYS = {
    "anisotropy", "cg_static_mask_file", "dm", "exchange", "jfile",
    "kfile", "momfile", "pdfile", "posfile", "qfile", "restartfile",
}


def tracked_files() -> set[str]:
    result = subprocess.run(
        ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
        check=True, stdout=subprocess.PIPE,
    )
    return {path for path in result.stdout.decode().split("\0") if path}


def tracked_path(path: Path, tracked: set[str]) -> bool:
    relative = path.as_posix()
    if relative in tracked:
        return relative in tracked
    prefix = relative.rstrip("/") + "/"
    return any(item.startswith(prefix) for item in tracked)


def input_references(case: Path) -> list[Path]:
    input_path = case / "inpsd.dat"
    references = [input_path]
    if not input_path.is_file():
        return references
    for line in input_path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0].lower() in INPUT_FILE_KEYS:
            references.append((case / fields[1]).resolve())
    return references


def main() -> int:
    tracked = tracked_files()
    missing: set[tuple[str, str]] = set()
    checked_files: set[str] = set()
    cases = [(E2E_ROOT / name, "either executable harness") for name in all_e2e_cases()]
    cases.extend((EXAMPLES_ROOT / name, "run_production_e2e.py") for name in EXAMPLE_CASES)

    for case, owner in cases:
        relative_case = case.relative_to(REPO_ROOT)
        if not tracked_path(relative_case, tracked):
            missing.add((relative_case.as_posix(), f"fixture directory ({owner})"))
        for path in input_references(case):
            relative = path.relative_to(REPO_ROOT)
            checked_files.add(relative.as_posix())
            if not tracked_path(relative, tracked):
                missing.add((relative.as_posix(), f"input reference from {relative_case}/inpsd.dat"))

    if missing:
        print("adaptive-CG fixture dependency audit: FAIL")
        for path, reason in sorted(missing):
            print(f"MISSING tracked path: {path} ({reason})")
        print(f"checked {len(cases)} fixture directories and {len(checked_files)} input paths")
        return 1
    print(f"adaptive-CG fixture dependency audit: PASS ({len(cases)} fixture directories, {len(checked_files)} input paths)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
