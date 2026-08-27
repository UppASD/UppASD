"""Concise Markdown report generator (WP-07 section H).

Turns one case/variant's `analysis.campaign_summary.build_case_summary`
output into Markdown, split into three clearly labelled kinds of claim so a
reader (or a later automated consumer) never has to guess which is which:

* **FACT** -- a measured number: an aggregate's median (with its MAD), or a
  ratio of two measured medians (`S_GPU_BESTCPU`, `S_GPU_1T`, `R_GPU_32_64`).
  Nothing under FACT is interpolated or extrapolated.
* **INTERPOLATION** -- a crossover threshold located strictly between two
  measured neighbouring points (`analysis.crossover.find_crossover`,
  `interpolated: True`). Never a number outside the tested range; the report
  says `below_tested_range`/`above_tested_range` instead of fabricating one
  when that is what was found.
* **HYPOTHESIS** -- a plain-language interpretation (e.g. "consistent with a
  bandwidth-limited DOUBLE-precision path") offered *because* the measured or
  interpolated facts above support raising it, always phrased as an
  interpretation and never restated as a measured claim.

`generate_markdown_report` is pure text generation over already-computed
values -- it does not itself compute anything `campaign_summary`/`crossover`/
`precision_metrics` did not already compute.
"""

from __future__ import annotations

import pathlib

from analysis import crossover as crossover_mod

# R_GPU_32_64 (T_DOUBLE / T_SINGLE) above this is called out as a hypothesis
# candidate for a bandwidth-limited DOUBLE-precision path. Chosen as "enough
# above the trivial 2x data-volume ratio to be worth a sentence", not derived
# from any measurement -- purely a threshold for when to *mention* the
# hypothesis in the generated text.
_PRECISION_PENALTY_HYPOTHESIS_THRESHOLD = 1.8


def _fmt(value, digits=4):
    if value is None:
        return "n/a"
    return f"{value:.{digits}g}"


def _fact_measurement_table(summary):
    lines = [
        "| precision | size | N_atom | T_GPU (s) | T_CPU-BEST (s) | T_CPU-1T (s) | S_GPU/BESTCPU | S_GPU/1T |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for precision in sorted(summary["by_precision"]):
        for row in summary["by_precision"][precision]["rows"]:
            gpu_mad = f" ± {_fmt(row.get('t_gpu_mad'))}" if row.get("t_gpu_mad") is not None else ""
            lines.append(
                f"| {precision} | {row['size_id']} | {row['natom']} | "
                f"{_fmt(row['t_gpu_seconds'])}{gpu_mad} | {_fmt(row['t_cpu_best_seconds'])} | "
                f"{_fmt(row['t_cpu_1t_seconds'])} | {_fmt(row['s_gpu_bestcpu'])} | {_fmt(row['s_gpu_1t'])} |"
            )
    return "\n".join(lines)


def _fact_precision_penalty_table(summary):
    penalties = summary.get("precision_penalty") or []
    if not penalties:
        return None
    lines = [
        "| size | R_GPU_32_64 (T_DOUBLE / T_SINGLE) | T_SINGLE (s) | T_DOUBLE (s) |",
        "|---|---|---|---|",
    ]
    for penalty in penalties:
        lines.append(
            f"| {penalty['size_id']} | {_fmt(penalty['r_gpu_32_64'])} | "
            f"{_fmt(penalty['t_single'])} | {_fmt(penalty['t_double'])} |"
        )
    return "\n".join(lines)


def _fact_section(summary):
    parts = [
        "## FACT",
        "",
        "Measured values only: aggregate medians (± MAD where available) and ratios of measured medians. "
        "No interpolation or extrapolation appears in this section.",
        "",
        _fact_measurement_table(summary),
    ]
    penalty_table = _fact_precision_penalty_table(summary)
    if penalty_table:
        parts += ["", "**Precision penalty (R_GPU_32_64):**", "", penalty_table]
    quality_lines = []
    for precision in sorted(summary["by_precision"]):
        status = summary["by_precision"][precision]["quality_status"]
        quality_lines.append(
            f"- {precision}: authoritative={status['authoritative']}, "
            f"quality_flags_union={status['quality_flags_union'] or 'none'}"
        )
    parts += ["", "**Quality status:**", ""] + quality_lines
    return "\n".join(parts)


def _interpolation_section(summary):
    lines = [
        "## INTERPOLATION",
        "",
        "Crossover thresholds located strictly between two measured neighbouring points "
        "(log-log interpolation; see `analysis.crossover`). A threshold outside the tested "
        "size range is reported as `below_tested_range`/`above_tested_range`, never as a "
        "fabricated crossover size.",
        "",
        "| precision | threshold | status | crossover N_atom | tested range |",
        "|---|---|---|---|---|",
    ]
    for precision in sorted(summary["by_precision"]):
        block = summary["by_precision"][precision]
        results = sorted(block["crossover"].values(), key=lambda result: result["threshold"])
        for result in results:
            if result["status"] == crossover_mod.WITHIN_TESTED_RANGE:
                crossover_display = (
                    f"~{result['crossover_natom']:g} (interpolated)" if result["interpolated"]
                    else f"{result['crossover_natom']:g} (measured)"
                )
            else:
                crossover_display = "-"
            lo, hi = result["tested_range"]
            lines.append(
                f"| {precision} | {result['threshold']}x | {result['status']} | "
                f"{crossover_display} | {lo:g} .. {hi:g} |"
            )
    return "\n".join(lines)


def _hypothesis_section(summary):
    sentences = []
    for precision in sorted(summary["by_precision"]):
        block = summary["by_precision"][precision]
        one_x = block["crossover"].get("1.0")
        if one_x is None:
            continue
        if one_x["status"] == crossover_mod.BELOW_TESTED_RANGE:
            sentences.append(
                f"- {precision}: GPU already exceeds CPU-BEST across the entire tested range "
                f"(1.0x crossover falls below `{one_x['tested_range'][0]:g}` atoms). This is "
                "*consistent with* per-launch/setup overhead already being amortized at the "
                "smallest tested size; the true crossover point has not itself been measured."
            )
        elif one_x["status"] == crossover_mod.ABOVE_TESTED_RANGE:
            sentences.append(
                f"- {precision}: GPU has not yet reached CPU-BEST parity within the tested range "
                f"(largest tested size `{one_x['tested_range'][1]:g}` atoms). This *may* indicate "
                "fixed per-step overhead (kernel launch, host/device synchronization) that is not "
                "yet amortized at the sizes tested -- a hypothesis, not a measured attribution."
            )
        five_x = block["crossover"].get("5.0")
        if five_x is not None and five_x["status"] == crossover_mod.ABOVE_TESTED_RANGE:
            sentences.append(
                f"- {precision}: the 5.0x speedup threshold is not reached within the tested "
                "range. Larger sizes would be needed to determine whether the speedup curve "
                "continues rising or has begun to plateau."
            )

    for penalty in summary.get("precision_penalty") or []:
        if penalty["r_gpu_32_64"] >= _PRECISION_PENALTY_HYPOTHESIS_THRESHOLD:
            sentences.append(
                f"- size {penalty['size_id']}: R_GPU_32_64 = {_fmt(penalty['r_gpu_32_64'])} is "
                "well above the 2x data-volume ratio DOUBLE precision implies on its own, which "
                "*is consistent with* a memory-bandwidth-limited DOUBLE-precision path -- an "
                "interpretation to investigate with kernel-level profiling, not a measured "
                "bandwidth figure."
            )

    if not sentences:
        sentences = [
            "- No hypothesis is raised here beyond the FACT/INTERPOLATION sections above; "
            "the measured and interpolated evidence did not surface a pattern worth calling out."
        ]

    return "\n".join([
        "## HYPOTHESIS",
        "",
        "Interpretive statements only. Every sentence below is explicitly hedged and must "
        "never be quoted as a measured or interpolated result.",
        "",
        *sentences,
    ])


def generate_markdown_report(summary):
    """Render `summary` (one case/variant's `campaign_summary.build_case_summary`
    output) as a Markdown report with FACT / INTERPOLATION / HYPOTHESIS sections.
    """
    identity = summary["identity"]
    header = (
        f"# GPU crossover report: {identity['case_id']} / {identity['variant_id']}\n\n"
        f"Campaign: `{identity['campaign_id']}` · Machine: `{identity['machine_id']}`"
    )
    return "\n\n".join([
        header,
        _fact_section(summary),
        _interpolation_section(summary),
        _hypothesis_section(summary),
    ]) + "\n"


def write_markdown_report(summary, path):
    """Write `generate_markdown_report(summary)` to `path`. Returns the path."""
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(generate_markdown_report(summary), encoding="utf-8")
    return path
