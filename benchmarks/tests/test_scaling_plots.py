"""Tests for the WP-07 campaign plots (section E, plots 1-2 and 5-9).

Like `test_omp_scaling_report.py`, plots are checked for existence/
non-emptiness, not pixel content.
"""

from __future__ import annotations

from analysis import campaign_summary, scaling_plots


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


def _entry(size_id, natom, requested_precision, *, t_gpu, t_cpu_best, t_cpu_1t,
           directed_interactions=None, fft_grid_points=None, workload_class="NEIGHBOR_LIST"):
    return {
        "size_id": size_id,
        "natom": natom,
        "directed_interactions": directed_interactions,
        "fft_grid_points": fft_grid_points,
        "workload_class": workload_class,
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


def _neighbor_summary():
    entries = [
        _entry("4096", 4096, "DOUBLE", t_gpu=1.0, t_cpu_best=0.5, t_cpu_1t=2.0, directed_interactions=200000),
        _entry("32768", 32768, "DOUBLE", t_gpu=2.0, t_cpu_best=4.0, t_cpu_1t=16.0, directed_interactions=1600000),
        _entry("4096", 4096, "SINGLE", t_gpu=0.6, t_cpu_best=0.5, t_cpu_1t=2.0, directed_interactions=200000),
        _entry("32768", 32768, "SINGLE", t_gpu=1.0, t_cpu_best=4.0, t_cpu_1t=16.0, directed_interactions=1600000),
    ]
    return campaign_summary.build_case_summary(entries)


def _fft_summary():
    entries = [
        _entry("16x16x16", 4096, "DOUBLE", t_gpu=0.02, t_cpu_best=0.05, t_cpu_1t=0.2,
               fft_grid_points=29791, workload_class="FFT_DIPOLE"),
        _entry("32x32x32", 32768, "DOUBLE", t_gpu=0.08, t_cpu_best=0.4, t_cpu_1t=1.6,
               fft_grid_points=250047, workload_class="FFT_DIPOLE"),
    ]
    return campaign_summary.build_case_summary(entries)


def _assert_valid_png(path):
    assert path is not None
    assert path.is_file()
    assert path.stat().st_size > 0


def test_plot_step_time_vs_natom(tmp_path):
    _assert_valid_png(scaling_plots.plot_step_time_vs_natom(_neighbor_summary(), tmp_path))


def test_plot_step_time_vs_ndirected(tmp_path):
    _assert_valid_png(scaling_plots.plot_step_time_vs_ndirected(_neighbor_summary(), tmp_path))


def test_plot_speedup_vs_natom(tmp_path):
    _assert_valid_png(scaling_plots.plot_speedup_vs_natom(_neighbor_summary(), tmp_path))


def test_plot_speedup_vs_ndirected(tmp_path):
    _assert_valid_png(scaling_plots.plot_speedup_vs_ndirected(_neighbor_summary(), tmp_path))


def test_plot_single_vs_double(tmp_path):
    _assert_valid_png(scaling_plots.plot_single_vs_double(_neighbor_summary(), tmp_path))


def test_plot_crossover_summary(tmp_path):
    _assert_valid_png(scaling_plots.plot_crossover_summary(_neighbor_summary(), tmp_path))


def test_plot_fft_vs_grid_size_returns_none_for_pure_neighbor_case(tmp_path):
    assert scaling_plots.plot_fft_vs_grid_size(_neighbor_summary(), tmp_path) is None


def test_plot_fft_vs_grid_size_plots_fft_case(tmp_path):
    _assert_valid_png(scaling_plots.plot_fft_vs_grid_size(_fft_summary(), tmp_path))


def test_plot_step_time_vs_ndirected_returns_none_for_fft_only_case(tmp_path):
    assert scaling_plots.plot_step_time_vs_ndirected(_fft_summary(), tmp_path) is None


def test_plot_campaign_report_writes_expected_plots_for_neighbor_case(tmp_path):
    written = scaling_plots.plot_campaign_report(_neighbor_summary(), tmp_path)
    names = {path.name for path in written}
    assert any(name.startswith("step_time_vs_natom__") for name in names)
    assert any(name.startswith("step_time_vs_ndirected__") for name in names)
    assert any(name.startswith("speedup_vs_natom__") for name in names)
    assert any(name.startswith("speedup_vs_ndirected__") for name in names)
    assert any(name.startswith("single_vs_double__") for name in names)
    assert any(name.startswith("crossover_summary__") for name in names)
    assert not any(name.startswith("fft_time_vs_grid__") for name in names)
    for path in written:
        _assert_valid_png(path)
