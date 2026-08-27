"""Tests for the WP-07 throughput metrics (sections C-D)."""

from __future__ import annotations

import pytest

from analysis import throughput


# ---------------------------------------------------------------------------
# spin_steps_per_second
# ---------------------------------------------------------------------------

def test_spin_steps_per_second():
    assert throughput.spin_steps_per_second(65536, 0.001) == pytest.approx(65536000.0)


def test_spin_steps_per_second_requires_positive_inputs():
    with pytest.raises(throughput.ThroughputError, match="natom"):
        throughput.spin_steps_per_second(0, 0.001)
    with pytest.raises(throughput.ThroughputError, match="t_step"):
        throughput.spin_steps_per_second(65536, 0.0)


# ---------------------------------------------------------------------------
# directed_interaction_visits_per_second
# ---------------------------------------------------------------------------

def test_directed_interaction_visits_per_second():
    value = throughput.directed_interaction_visits_per_second(
        3801088, 0.001, workload_class="NEIGHBOR_LIST"
    )
    assert value == pytest.approx(3801088000.0)


def test_directed_interaction_visits_per_second_rejects_fft_workload():
    with pytest.raises(throughput.ThroughputError, match="only applies to"):
        throughput.directed_interaction_visits_per_second(100, 0.001, workload_class="FFT_DIPOLE")


def test_directed_interaction_visits_per_second_requires_value():
    with pytest.raises(throughput.ThroughputError, match="required"):
        throughput.directed_interaction_visits_per_second(None, 0.001, workload_class="NEIGHBOR_LIST")


# ---------------------------------------------------------------------------
# fft_grid_points_per_second
# ---------------------------------------------------------------------------

def test_fft_grid_points_per_second():
    value = throughput.fft_grid_points_per_second(262144, 0.01, workload_class="FFT_DIPOLE")
    assert value == pytest.approx(26214400.0)


def test_fft_grid_points_per_second_rejects_neighbor_workload():
    with pytest.raises(throughput.ThroughputError, match="only applies to"):
        throughput.fft_grid_points_per_second(262144, 0.01, workload_class="NEIGHBOR_LIST")


# ---------------------------------------------------------------------------
# compute_throughput
# ---------------------------------------------------------------------------

def test_compute_throughput_neighbor_list_workload():
    result = throughput.compute_throughput(
        natom=65536, t_step=0.001, workload_class="NEIGHBOR_LIST", directed_interactions=3801088,
    )
    assert set(result) == {"spin_steps_per_second", "directed_interaction_visits_per_second"}


def test_compute_throughput_fft_dipole_workload():
    result = throughput.compute_throughput(
        natom=4096, t_step=0.01, workload_class="FFT_DIPOLE", fft_grid_points=262144,
    )
    assert set(result) == {"spin_steps_per_second", "fft_grid_points_per_second"}


def test_compute_throughput_neighbor_plus_fft_workload_gets_both():
    result = throughput.compute_throughput(
        natom=4096, t_step=0.01, workload_class="NEIGHBOR_LIST_PLUS_FFT_DIPOLE",
        directed_interactions=200000, fft_grid_points=262144,
    )
    assert set(result) == {
        "spin_steps_per_second", "directed_interaction_visits_per_second", "fft_grid_points_per_second",
    }


def test_compute_throughput_rejects_unknown_workload_class():
    with pytest.raises(throughput.ThroughputError, match="unknown workload_class"):
        throughput.compute_throughput(natom=1, t_step=1.0, workload_class="BOGUS")


# ---------------------------------------------------------------------------
# compute_throughput_for_aggregate
# ---------------------------------------------------------------------------

def test_compute_throughput_for_aggregate_uses_median_as_t_step():
    aggregate = {"metric": "steady_step_seconds", "median": 0.002}
    result = throughput.compute_throughput_for_aggregate(
        aggregate, natom=1000, workload_class="NEIGHBOR_LIST", directed_interactions=5000,
    )
    assert result["spin_steps_per_second"] == pytest.approx(1000 / 0.002)


def test_compute_throughput_for_aggregate_requires_steady_step_metric():
    aggregate = {"metric": "process_wall_seconds", "median": 0.002}
    with pytest.raises(throughput.ThroughputError, match="steady_step_seconds"):
        throughput.compute_throughput_for_aggregate(aggregate, natom=1000, workload_class="NEIGHBOR_LIST")
