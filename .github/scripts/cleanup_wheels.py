#!/usr/bin/env python3
"""
Cleanup old UppASD wheels in a directory.

Rules:
- Group wheels by (python interpreter, platform)
- Keep ALL stable (non-prerelease) wheels
- Keep only the latest N prerelease/dev wheels per group
"""

from pathlib import Path
from packaging.version import Version
from packaging.utils import parse_wheel_filename
import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description="Cleanup old UppASD wheels")
    parser.add_argument(
        "--wheels-dir",
        default="wheels",
        help="Directory containing wheels (default: wheels)",
    )
    parser.add_argument(
        "--keep-dev",
        type=int,
        default=3,
        help="Number of dev/prerelease wheels to keep per (python, platform)",
    )
    args = parser.parse_args()

    wheels_dir = Path(args.wheels_dir)
    keep_dev = args.keep_dev

    if not wheels_dir.exists():
        print(f"ERROR: wheels directory not found: {wheels_dir}", file=sys.stderr)
        sys.exit(1)

    # (py, platform) -> list of (Version, Path)
    groups: dict[tuple[str, str], list[tuple[Version, Path]]] = {}

    for whl in wheels_dir.glob("*.whl"):
        try:
            name, version_str, build, tags = parse_wheel_filename(whl.name)
        except Exception as e:
            print(f"Skipping unparseable wheel: {whl.name} ({e})")
            continue

        version = Version(version_str)

        # A wheel may have multiple tags; treat each independently
        for tag in tags:
            key = (tag.interpreter, tag.platform)
            groups.setdefault(key, []).append((version, whl))

    keep: set[Path] = set()

    for (py, plat), items in groups.items():
        # Sort by version (PEP 440)
        items.sort(key=lambda x: x[0])

        stable = [w for v, w in items if not v.is_prerelease]
        dev = [(v, w) for v, w in items if v.is_prerelease]

        # Keep all stable wheels
        keep.update(stable)

        # Keep latest N dev wheels
        keep.update(w for _, w in dev[-keep_dev:])

        print(
            f"[{py} | {plat}] "
            f"keeping {len(stable)} stable + {min(len(dev), keep_dev)} dev"
        )

    # Remove everything else
    removed = 0
    for whl in wheels_dir.glob("*.whl"):
        if whl not in keep:
            print(f"Removing old wheel: {whl.name}")
            whl.unlink()
            removed += 1

    print(f"Cleanup complete: removed {removed} wheel(s)")


if __name__ == "__main__":
    main()
