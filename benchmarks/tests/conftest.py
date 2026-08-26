"""Test configuration for the benchmark data-contract tests."""

import pathlib
import sys

BENCHMARKS_DIR = pathlib.Path(__file__).resolve().parent.parent

# The harness is not installed as a package; make it importable from the tests.
if str(BENCHMARKS_DIR) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS_DIR))
