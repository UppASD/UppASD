"""Tests for the WP-07 Markdown report generator (section H)."""

from __future__ import annotations

from analysis import campaign_summary, markdown_report


def _aggregate(backend, median, *, size_id, aggregate_id, omp_threads=1, mad=0.05, **overrides):
    base = {
        "aggregate_id": aggregate_id,
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": size_id,
        "machine_id": "host1",
        "measurement_profile": "DYNAMICS_ONLY",
        "metric": "steady_step_seconds",
        "backend": backend,
        "omp_threads": omp_threads,
        "median": median,
        "mad": mad,
        "authoritative": True,
        "quality_flags_union": [],
        "all_samples_numerical_valid": True,
        "all_samples_environment_valid": True,
    }
    base.update(overrides)
    return base


def _entry(size_id, natom, requested_precision, *, t_gpu, t_cpu_best, t_cpu_1t, directed_interactions=None):
    return {
        "size_id": size_id,
        "natom": natom,
        "directed_interactions": directed_interactions,
        "fft_grid_points": None,
        "workload_class": "NEIGHBOR_LIST",
        "requested_precision": requested_precision,
        "gpu_aggregate": _aggregate(
            "GPU", t_gpu, size_id=size_id, aggregate_id=f"gpu-{requested_precision}-{size_id}",
            gpu_backend="CUDA", requested_precision=requested_precision,
        ),
        "cpu_best_aggregate": _aggregate(
            "CPU", t_cpu_best, size_id=size_id, aggregate_id=f"cpubest-{size_id}",
            omp_threads=8, requested_precision="DOUBLE",
        ),
        "cpu_1t_aggregate": _aggregate(
            "CPU", t_cpu_1t, size_id=size_id, aggregate_id=f"cpu1t-{size_id}",
            omp_threads=1, requested_precision="DOUBLE",
        ),
    }


def _summary_above_range():
    """GPU never reaches CPU-BEST within the tested range -- exercises ABOVE_TESTED_RANGE."""
    entries = [
        _entry("4096", 4096, "DOUBLE", t_gpu=4.0, t_cpu_best=2.0, t_cpu_1t=8.0, directed_interactions=200000),
        _entry("32768", 32768, "DOUBLE", t_gpu=8.0, t_cpu_best=6.0, t_cpu_1t=24.0, directed_interactions=1600000),
    ]
    return campaign_summary.build_case_summary(entries)


def _summary_below_range():
    """GPU already exceeds CPU-BEST at the smallest tested size -- exercises BELOW_TESTED_RANGE."""
    entries = [
        _entry("4096", 4096, "DOUBLE", t_gpu=0.1, t_cpu_best=1.0, t_cpu_1t=4.0, directed_interactions=200000),
        _entry("32768", 32768, "DOUBLE", t_gpu=0.2, t_cpu_best=4.0, t_cpu_1t=16.0, directed_interactions=1600000),
    ]
    return campaign_summary.build_case_summary(entries)


def _summary_with_precision_penalty():
    entries = [
        _entry("4096", 4096, "DOUBLE", t_gpu=2.0, t_cpu_best=1.0, t_cpu_1t=4.0, directed_interactions=200000),
        _entry("4096", 4096, "SINGLE", t_gpu=0.5, t_cpu_best=1.0, t_cpu_1t=4.0, directed_interactions=200000),
    ]
    return campaign_summary.build_case_summary(entries)


def test_report_contains_all_three_sections_in_order():
    text = markdown_report.generate_markdown_report(_summary_above_range())
    fact_index = text.index("## FACT")
    interp_index = text.index("## INTERPOLATION")
    hypothesis_index = text.index("## HYPOTHESIS")
    assert fact_index < interp_index < hypothesis_index


def test_report_header_names_case_and_variant():
    text = markdown_report.generate_markdown_report(_summary_above_range())
    assert "B01_bccFe / T0" in text
    assert "camp1" in text
    assert "host1" in text


def test_fact_section_reports_measured_medians():
    text = markdown_report.generate_markdown_report(_summary_above_range())
    fact_section = text.split("## FACT")[1].split("## INTERPOLATION")[0]
    assert "4096" in fact_section
    assert "32768" in fact_section


def test_interpolation_section_reports_above_tested_range_without_fabricated_natom():
    summary = _summary_above_range()
    text = markdown_report.generate_markdown_report(summary)
    interp_section = text.split("## INTERPOLATION")[1].split("## HYPOTHESIS")[0]
    assert "above_tested_range" in interp_section


def test_interpolation_section_reports_below_tested_range():
    summary = _summary_below_range()
    text = markdown_report.generate_markdown_report(summary)
    interp_section = text.split("## INTERPOLATION")[1].split("## HYPOTHESIS")[0]
    assert "below_tested_range" in interp_section


def test_hypothesis_section_flags_precision_penalty():
    text = markdown_report.generate_markdown_report(_summary_with_precision_penalty())
    hypothesis_section = text.split("## HYPOTHESIS")[1]
    assert "bandwidth" in hypothesis_section.lower()


def test_hypothesis_section_is_hedged_language():
    text = markdown_report.generate_markdown_report(_summary_above_range())
    hypothesis_section = text.split("## HYPOTHESIS")[1]
    assert "interpretation" in hypothesis_section.lower() or "hypothesis" in hypothesis_section.lower()


def test_write_markdown_report(tmp_path):
    summary = _summary_above_range()
    path = markdown_report.write_markdown_report(summary, tmp_path / "report.md")
    assert path.is_file()
    assert path.read_text(encoding="utf-8") == markdown_report.generate_markdown_report(summary)


# ---------------------------------------------------------------------------
# WP-09: LEAN campaign non-authoritative banner
# ---------------------------------------------------------------------------

_LEAN_BANNER = "LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS"


def test_no_banner_by_default():
    text = markdown_report.generate_markdown_report(_summary_above_range())
    assert _LEAN_BANNER not in text


def test_banner_renders_as_first_line_when_given():
    text = markdown_report.generate_markdown_report(_summary_above_range(), report_banner=_LEAN_BANNER)
    assert text.startswith(f"> # ⚠ {_LEAN_BANNER} ⚠")
    assert text.index(_LEAN_BANNER) < text.index("# GPU crossover report")


def test_write_markdown_report_carries_banner_through(tmp_path):
    summary = _summary_above_range()
    path = markdown_report.write_markdown_report(summary, tmp_path / "lean_report.md", report_banner=_LEAN_BANNER)
    assert _LEAN_BANNER in path.read_text(encoding="utf-8")
