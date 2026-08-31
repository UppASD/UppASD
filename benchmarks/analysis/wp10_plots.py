"""WP-10 campaign plots: renders every figure `analysis.scaling_plots`
(GPU-crossover plots 1-2, 5-9) and `analysis.omp_scaling_report`
(OpenMP scaling/efficiency plots 3-4) already know how to draw, fed from
wp10_2026-08-28's aggregate records via `analysis.wp10_summary`'s loaders.

Consumes records (raw samples + aggregates); produces PNGs under
`--out-dir` (gitignored per benchmarks/.gitignore's "analysis/out/" rule,
same as every other generated analysis product) plus one small
`plots_manifest.json` listing what was written and why anything was
skipped (a case/variant with no GPU data has nothing for
`plot_speedup_vs_natom` to draw, for instance -- skipped, not errored).
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from analysis import omp_scaling_report
from analysis import scaling_plots
from analysis import wp10_summary


def build_size_families(aggregates, case_id, variant_id, wm_by_cell):
    """One `omp_scaling_report.build_scaling_summary` size-family dict per
    size of (case_id, variant_id): that size's CPU steady_step_seconds
    thread-count sweep plus its natom."""
    sizes = sorted({
        a["size_id"] for a in aggregates
        if a["case_id"] == case_id and a["variant_id"] == variant_id and a["backend"] == "CPU"
        and a["metric"] == "steady_step_seconds"
    })
    families = []
    for size_id in sizes:
        cpu_aggs = [
            a for a in aggregates
            if a["case_id"] == case_id and a["variant_id"] == variant_id and a["size_id"] == size_id
            and a["backend"] == "CPU" and a["metric"] == "steady_step_seconds"
        ]
        wm = wm_by_cell.get((case_id, variant_id, size_id))
        if not cpu_aggs or not wm:
            continue
        families.append({"size_id": size_id, "natom": wm["natom"], "aggregates": cpu_aggs})
    return families


def generate_all_plots(results_dir, out_dir):
    aggregates = wp10_summary.load_aggregates(results_dir)
    records = wp10_summary.load_raw_records(results_dir)
    wm_by_cell_full = wp10_summary.workload_metadata_by_cell(records)
    wm_by_cell = {key: wm for key, wm in wm_by_cell_full.items()}

    case_variants = sorted({(a["case_id"], a["variant_id"]) for a in aggregates})
    manifest = {"omp_scaling": {}, "crossover": {}, "skipped": []}
    out_dir = pathlib.Path(out_dir)

    for case_id, variant_id in case_variants:
        slug = f"{case_id}__{variant_id}"

        # OpenMP scaling/efficiency plots (3-4): CPU-only, no GPU needed.
        # plot_omp_scaling names files by size_id alone (e.g.
        # "omp_scaling__20x20x20.png", "p_best_vs_natom.png"); several cases
        # in this campaign share size_ids (B01's own two variants both use
        # 20x20x20/32x32x32/40x40x40; B03 and B04 both use 25x25x25; every
        # case/variant shares the literal "p_best_vs_natom.png"), which
        # silently overwrote each other when they all wrote into one shared
        # directory. Each case/variant gets its own subdirectory instead.
        families = build_size_families(aggregates, case_id, variant_id, wm_by_cell)
        if families:
            try:
                omp_summary = omp_scaling_report.build_scaling_summary(families)
                paths = omp_scaling_report.plot_omp_scaling(omp_summary, out_dir / "omp_scaling" / slug)
                manifest["omp_scaling"][slug] = [str(p) for p in paths]
            except Exception as exc:  # noqa: BLE001
                manifest["skipped"].append({"slug": slug, "stage": "omp_scaling", "error": str(exc)})
        else:
            manifest["skipped"].append({"slug": slug, "stage": "omp_scaling", "error": "no CPU size families"})

    # GPU-crossover plots (1-2, 5-9): reuse wp10_summary's own
    # campaign_summary.build_case_summary output directly.
    full_summary = wp10_summary.build_full_summary(results_dir)
    for slug, cv_summary in full_summary["case_variant_summaries"].items():
        if "error" in cv_summary:
            manifest["skipped"].append({"slug": slug, "stage": "crossover", "error": cv_summary["error"]})
            continue
        try:
            paths = scaling_plots.plot_campaign_report(cv_summary, out_dir / "crossover")
            manifest["crossover"][slug] = [str(p) for p in paths]
        except Exception as exc:  # noqa: BLE001
            manifest["skipped"].append({"slug": slug, "stage": "crossover", "error": str(exc)})

    manifest_path = out_dir / "plots_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    manifest = generate_all_plots(args.results_dir, args.out_dir)
    n_omp = sum(len(v) for v in manifest["omp_scaling"].values())
    n_cross = sum(len(v) for v in manifest["crossover"].values())
    print(f"wrote {n_omp} OMP-scaling PNGs across {len(manifest['omp_scaling'])} case/variants")
    print(f"wrote {n_cross} crossover PNGs across {len(manifest['crossover'])} case/variants")
    if manifest["skipped"]:
        print(f"skipped {len(manifest['skipped'])} (case/variant, stage) pairs:")
        for entry in manifest["skipped"]:
            print(" ", entry)
