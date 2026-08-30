"""Post-hoc aggregation over WP-10's raw sample records.

wp10_driver.py streamed one steady_step_seconds aggregate per config as the
campaign ran, but never a matching setup_seconds aggregate (needed for
section H's setup-vs-steady-state economics) -- both are derivable from the
same FitResults it already computed and threw away. Rather than re-run any
UppASD executable, this module re-reads the raw sample records the driver
already wrote (every (nstep, process_wall_seconds) point is there), regroups
them exactly the way the driver did (by cell + config + sample_index),
re-fits via steady_state.fit_multi_nstep, and writes a complete aggregate set
-- both metrics, every config -- to a fresh directory. This is pure
post-processing: no new measurement.
"""

from __future__ import annotations

import argparse
import collections
import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from harness import aggregate as aggregate_mod
from harness import steady_state

CELL_FIELDS = (
    "campaign_id", "case_id", "variant_id", "size_id", "machine_id",
    "build_id", "backend", "gpu_backend", "omp_threads", "requested_precision",
    "measurement_profile",
)


def load_raw_records(results_dir):
    results_dir = pathlib.Path(results_dir)
    records = []
    for path in results_dir.glob("*.json"):
        if path.name == "driver_log.jsonl":
            continue
        records.append(json.loads(path.read_text()))
    return records


_RUN_ID_NSTEP_RE = re.compile(r"__n(\d+)$")


def _true_nstep(record):
    """The nstep this sample actually ran, from wp10_driver.py's own run_id
    suffix (``..._n{nstep}``) -- not ``record["nstep"]``.

    A now-fixed case-sensitivity bug in cases.py's `_apply_overrides`
    (differently-cased template keywords, e.g. this project's lowercase
    `nstep` templates vs. the harness's `{"Nstep": ...}` override dict) meant
    every override was appended as a second line rather than replacing the
    first. UppASD's own parser lowercases keywords before dispatch and
    applies whatever it reads last, so the *simulation* still ran the
    correct, intended nstep -- but `read_keyword`/`runner.run_sample`'s
    `record["nstep"]` field read the now-stale, unreachable first line, so
    every raw sample in this campaign has an honest `process_wall_seconds`
    behind a misleading `nstep` value. run_id's suffix reflects the value
    wp10_driver.py's own loop variable actually requested, independent of
    that bug, so it is used here instead.
    """
    match = _RUN_ID_NSTEP_RE.search(record["run_id"])
    if not match:
        raise ValueError(f"run_id has no __n<nstep> suffix: {record['run_id']!r}")
    return int(match.group(1))


def group_by_config(records):
    """{(case_id, variant_id, size_id, backend, omp_threads, requested_precision):
    {sample_index: [(nstep, process_wall_seconds, record), ...]}}"""
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for record in records:
        if record.get("run_status") != "COMPLETED":
            continue
        key = (
            record["case_id"], record["variant_id"], record["size_id"], record["backend"],
            record["omp_threads"], record["requested_precision"],
        )
        groups[key][record["sample_index"]].append(
            (_true_nstep(record), record["process_wall_seconds"], record)
        )
    return groups


def rebuild_aggregates(results_dir, out_dir):
    records = load_raw_records(results_dir)
    groups = group_by_config(records)

    written = []
    skipped = []
    for config_key, by_sample in groups.items():
        case_id, variant_id, size_id, backend, omp_threads, precision = config_key
        fits = []
        fit_run_ids = []
        contributing_records = []
        for sample_index, points in sorted(by_sample.items()):
            distinct = {nstep for nstep, _, _ in points}
            if len(distinct) < steady_state.MIN_DISTINCT_NSTEP_FOR_REGRESSION:
                continue
            try:
                fit = steady_state.fit_multi_nstep([(n, t) for n, t, _ in points])
            except steady_state.SteadyStateError:
                continue
            fits.append(fit)
            lowest_nstep_record = min(points, key=lambda p: p[0])[2]
            fit_run_ids.append(lowest_nstep_record["run_id"])
            contributing_records.extend(r for _, _, r in points)

        if not fits:
            skipped.append(config_key)
            continue

        cell = {field: contributing_records[-1][field] for field in CELL_FIELDS}
        quality_flags_union = set()
        for r in contributing_records:
            quality_flags_union.update(r.get("quality_flags", []))
        all_numerical_valid = all(r.get("numerical_valid") for r in contributing_records)
        all_environment_valid = all(r.get("environment_valid") for r in contributing_records)

        for metric in ("steady_step_seconds", "setup_seconds"):
            agg_id = (
                f"{cell['campaign_id']}__{backend}__{case_id}__{variant_id}__{size_id}__"
                f"{'t' + str(omp_threads) if backend == 'CPU' else precision}__{metric}"
            )
            agg = aggregate_mod.build_fit_aggregate(
                fits, metric=metric, run_ids=fit_run_ids, cell=cell, aggregate_id=agg_id,
                all_numerical_valid=all_numerical_valid, all_environment_valid=all_environment_valid,
                quality_flags_union=quality_flags_union,
            )
            out_path = pathlib.Path(out_dir) / f"{agg_id}.json"
            out_path.parent.mkdir(parents=True, exist_ok=True)
            out_path.write_text(json.dumps(agg, indent=2, sort_keys=True) + "\n")
            written.append(agg_id)

    return written, skipped


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    written, skipped = rebuild_aggregates(args.results_dir, args.out_dir)
    print(f"wrote {len(written)} aggregates ({len(written)//2} configs x 2 metrics)")
    print(f"skipped {len(skipped)} configs with zero valid fits:")
    for key in skipped:
        print(" ", key)
