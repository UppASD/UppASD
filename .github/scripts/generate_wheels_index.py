#!/usr/bin/env python3
"""
Generate a simple HTML index for pip --find-links.

Sorting rules (simple and robust):
  1) release        (string)
  2) dev/tag        (string)
  3) python tag:
       - cp3X  before cp3XX
       - alphabetical within each group
  4) platform       (string)
"""

from pathlib import Path
import argparse
import sys
from datetime import datetime


def parse_wheel_name(filename: str):
    """
    Parse wheel filename into sortable components.

    Expected format:
      uppasd-<release>.<tag>-<py>-<abi>-<platform>.whl
    """
    name = filename.removesuffix(".whl")
    parts = name.split("-")

    if len(parts) < 5:
        raise ValueError(f"Unrecognized wheel name: {filename}")

    _, version, py_tag, abi_tag, platform = parts[0], parts[1], parts[2], parts[3], "-".join(parts[4:])

    # Split version into release + dev/tag (string split, no semantics)
    if ".dev" in version:
        release, dev = version.split(".dev", 1)
        dev = "dev" + dev
    else:
        release = version
        dev = ""

    return release, dev, py_tag, platform


def python_sort_key(py_tag: str):
    """
    Sort cp3X before cp3XX, then alphabetically.
    """
    return (len(py_tag), py_tag)


def sort_key(path: Path):
    if path.suffix == ".whl":
        try:
            release, dev, py_tag, platform = parse_wheel_name(path.name)
            return (
                0,                # wheels first
                release,
                dev,
                python_sort_key(py_tag),
                platform,
            )
        except Exception:
            return (0, "", "", (99, ""), "")
    else:
        # sdists last
        return (1, path.name)


def main():
    parser = argparse.ArgumentParser(
        description="Generate simple HTML index for pip --find-links"
    )
    parser.add_argument(
        "--wheels-dir",
        default="wheels",
        help="Directory containing wheels and sdists (default: wheels)",
    )
    args = parser.parse_args()

    root = Path(args.wheels_dir)
    index_file = root / "index.html"

    if not root.exists():
        print(f"ERROR: '{root}' directory does not exist", file=sys.stderr)
        sys.exit(1)

    artifacts = sorted(
        list(root.glob("*.whl")) + list(root.glob("*.tar.gz")),
        key=sort_key,
    )

    if not artifacts:
        print("ERROR: No distributions found", file=sys.stderr)
        sys.exit(1)

    now = datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        "<!DOCTYPE html>",
        "<html>",
        "  <head>",
        '    <meta charset="utf-8">',
        "    <title>UppASD distributions</title>",
        "  </head>",
        "  <body>",
        "    <h1>UppASD distributions</h1>",
        f"    <p>Generated {now}</p>",
        "    <ul>",
    ]

    for path in artifacts:
        lines.append(f'      <li><a href="{path.name}">{path.name}</a></li>')

    lines.extend([
        "    </ul>",
        "  </body>",
        "</html>",
        "",
    ])

    index_file.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {index_file}")
    print(f"Indexed {len(artifacts)} file(s):")
    for a in artifacts:
        print(f"  - {a.name}")


if __name__ == "__main__":
    main()
