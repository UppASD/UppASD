#!/usr/bin/env python3
"""Validate UppASD benchmark records against the versioned result contract.

Two layers of checking are applied:

1. JSON Schema (``benchmarks/schema/``), which owns field presence, types,
   enumerations and the conditional rules that can be expressed structurally.
2. Cross-field consistency checks, for the handful of invariants that JSON
   Schema cannot express (products, cardinalities, orderings).

Usage:

    python3 -m benchmarks.harness.schema_validate record1.json record2.json
    python3 benchmarks/harness/schema_validate.py results/**/*.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import pathlib
import sys

import jsonschema

SCHEMA_DIR = pathlib.Path(__file__).resolve().parent.parent / "schema"

SCHEMA_FILES = {
    "raw_sample": SCHEMA_DIR / "benchmark_record.v1.schema.json",
    "aggregate": SCHEMA_DIR / "benchmark_aggregate.v1.schema.json",
}


class RecordInvalid(ValueError):
    """A record violates the benchmark result contract."""


def load_schema(record_kind):
    """Return the parsed JSON Schema for ``record_kind``."""
    try:
        path = SCHEMA_FILES[record_kind]
    except KeyError:
        raise RecordInvalid(
            "unknown record_kind %r (expected one of %s)"
            % (record_kind, ", ".join(sorted(SCHEMA_FILES)))
        )
    with open(path) as handle:
        return json.load(handle)


def _cross_field_errors(record):
    """Return contract violations that JSON Schema cannot express."""
    errors = []
    kind = record.get("record_kind")

    if kind == "raw_sample":
        stamp = record.get("timestamp_utc")
        if isinstance(stamp, str):
            try:
                datetime.datetime.fromisoformat(stamp.replace("Z", "+00:00"))
            except ValueError:
                errors.append("timestamp_utc (%r) is not an ISO 8601 timestamp" % stamp)
        grid = record.get("fft_grid_padded") or record.get("fft_grid")
        points = record.get("fft_grid_points")
        if grid and points is not None:
            expected = math.prod(grid)
            if expected != points:
                errors.append(
                    "fft_grid_points (%d) is not the product of the transformed "
                    "grid %s (%d)" % (points, grid, expected)
                )
        mean_n = record.get("mean_neighbors")
        max_n = record.get("max_neighbors")
        if mean_n is not None and max_n is not None and mean_n > max_n:
            errors.append(
                "mean_neighbors (%s) exceeds max_neighbors (%s)" % (mean_n, max_n)
            )
        natom = record.get("natom")
        directed = record.get("directed_interactions")
        if natom and directed is not None and mean_n is not None:
            if not math.isclose(directed, natom * mean_n, rel_tol=1e-6):
                errors.append(
                    "directed_interactions (%s) is inconsistent with "
                    "natom * mean_neighbors (%s)" % (directed, natom * mean_n)
                )

    elif kind == "aggregate":
        run_ids = record.get("run_ids") or []
        count = record.get("sample_count")
        if count is not None and count != len(run_ids):
            errors.append(
                "sample_count (%s) does not match the number of run_ids (%d)"
                % (count, len(run_ids))
            )
        lo, med, hi = (
            record.get("minimum"),
            record.get("median"),
            record.get("maximum"),
        )
        if None not in (lo, med, hi) and not lo <= med <= hi:
            errors.append(
                "statistics are not ordered: minimum (%s) <= median (%s) <= "
                "maximum (%s) is violated" % (lo, med, hi)
            )

    return errors


def validate_record(record):
    """Validate one parsed record. Raises :class:`RecordInvalid` on failure."""
    if not isinstance(record, dict):
        raise RecordInvalid("record must be a JSON object")
    schema = load_schema(record.get("record_kind"))
    validator = jsonschema.Draft202012Validator(schema)
    failures = [
        "%s: %s" % ("/".join(str(p) for p in error.absolute_path) or "<root>", error.message)
        for error in validator.iter_errors(record)
    ]
    failures.extend(_cross_field_errors(record))
    if failures:
        raise RecordInvalid("; ".join(failures))
    return True


def validate_file(path):
    """Validate one record file."""
    with open(path) as handle:
        try:
            record = json.load(handle)
        except json.JSONDecodeError as exc:
            raise RecordInvalid("not valid JSON: %s" % exc)
    return validate_record(record)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", help="record files to validate")
    args = parser.parse_args(argv)

    failed = 0
    for path in args.paths:
        try:
            validate_file(path)
        except RecordInvalid as exc:
            failed += 1
            print("INVALID %s: %s" % (path, exc), file=sys.stderr)
        else:
            print("ok      %s" % path)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
