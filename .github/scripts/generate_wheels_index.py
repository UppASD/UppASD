#!/usr/bin/env python3
"""
Generate a simple HTML index for pip --find-links.

Scans the 'wheels/' directory for *.whl files and writes wheels/index.html.
"""

from pathlib import Path
import sys

WHEELS_DIR = Path("wheels")
INDEX_FILE = WHEELS_DIR / "index.html"

def main():
    if not WHEELS_DIR.exists():
        print(f"ERROR: '{WHEELS_DIR}' directory does not exist", file=sys.stderr)
        sys.exit(1)

    wheels = sorted(WHEELS_DIR.glob("*.whl"))

    print(f"Found {len(wheels)} .whl files in {WHEELS_DIR}/")
    if not wheels:
        print("ERROR: No .whl files found in wheels/", file=sys.stderr)
        sys.exit(1)

    lines = [
        "<!DOCTYPE html>",
        "<html>",
        "  <head>",
        '    <meta charset="utf-8">',
        "    <title>UppASD wheels</title>",
        "  </head>",
        "  <body>",
        "    <h1>UppASD wheels</h1>",
        "    <ul>",
    ]

    for whl in wheels:
        name = whl.name
        lines.append(f'      <li><a href="{name}">{name}</a></li>')

    lines.extend([
        "    </ul>",
        "  </body>",
        "</html>",
        "",
    ])

    INDEX_FILE.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {INDEX_FILE} with {len(wheels)} wheel(s):")
    for w in wheels:
        print(f"  - {w.name}")

if __name__ == "__main__":
    main()
