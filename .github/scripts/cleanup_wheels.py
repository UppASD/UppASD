from pathlib import Path
from packaging.version import Version
import re

# How many dev wheels to keep per (python, platform)
KEEP_DEV = 3

WHEELS_DIR = Path("wheels")

if not WHEELS_DIR.exists():
    raise SystemExit("No wheels directory found")

wheel_re = re.compile(
    r"""
    uppasd-
    (?P<version>[^-]+)-
    (?P<py>cp\d+)-
    (?P<abi>[^-]+)-
    (?P<plat>[^.]+)
    """,
    re.VERBOSE,
)

groups = {}

for whl in WHEELS_DIR.glob("*.whl"):
    m = wheel_re.search(whl.name)
    if not m:
        continue

    key = (m["py"], m["plat"])
    groups.setdefault(key, []).append((Version(m["version"]), whl))

keep = set()

for (py, plat), items in groups.items():
    # sort by version (PEP 440)
    items.sort(key=lambda x: x[0])

    stable = [w for v, w in items if not v.is_prerelease]
    dev = [(v, w) for v, w in items if v.is_prerelease]

    # keep all stable wheels
    keep.update(stable)

    # keep only latest N dev wheels
    keep.update(w for _, w in dev[-KEEP_DEV:])

for whl in WHEELS_DIR.glob("*.whl"):
    if whl not in keep:
        print(f"Removing old wheel: {whl.name}")
        whl.unlink()
