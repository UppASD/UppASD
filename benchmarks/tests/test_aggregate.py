"""WP-04 section F: sample-statistics aggregation tests."""

from __future__ import annotations

import copy

import pytest

from harness import aggregate
from harness import steady_state


def _base_record(**overrides):
    record = {
        "run_id": "run-0",
        "campaign_id": "camp1",
        "case_id": "B01_bccFe",
        "variant_id": "T0",
        "size_id": "16x16x16",
        "machine_id": "m1",
        "build_id": "b1",
        "backend": "CPU",
        "gpu_backend": None,
        "omp_threads": 4,
        "requested_precision": "DOUBLE",
        "measurement_profile": "DYNAMICS_ONLY",
        "nstep": 1000,
        "run_status": "COMPLETED",
        "process_wall_seconds": 1.0,
        "numerical_valid": True,
        "environment_valid": True,
        "quality_flags": [],
    }
    record.update(overrides)
    return record


def test_compute_sample_statistics_known_values():
    stats = aggregate.compute_sample_statistics([1.0, 2.0, 3.0, 100.0])
    assert stats.count == 4
    assert stats.median == 2.5
    assert stats.minimum == 1.0
    assert stats.maximum == 100.0
    assert stats.mad == pytest.approx(1.0)


def test_compute_sample_statistics_rejects_empty():
    with pytest.raises(aggregate.AggregateError):
        aggregate.compute_sample_statistics([])


def test_build_process_wall_aggregate_happy_path():
    records = [
        _base_record(run_id="r0", process_wall_seconds=1.0),
        _base_record(run_id="r1", process_wall_seconds=1.2),
        _base_record(run_id="r2", process_wall_seconds=0.9),
    ]
    result = aggregate.build_process_wall_aggregate(records, aggregate_id="agg-1")

    assert result["record_kind"] == "aggregate"
    assert result["metric"] == "process_wall_seconds"
    assert result["sample_count"] == 3
    assert result["run_ids"] == ["r0", "r1", "r2"]
    assert result["minimum"] == 0.9
    assert result["maximum"] == 1.2
    assert result["median"] == 1.0
    assert result["nstep"] == 1000
    assert result["all_samples_numerical_valid"] is True
    assert result["all_samples_environment_valid"] is True
    assert result["quality_flags_union"] == []
    assert result["authoritative"] is True


def test_build_process_wall_aggregate_excludes_failed_samples():
    records = [
        _base_record(run_id="r0", process_wall_seconds=1.0),
        _base_record(
            run_id="r1", run_status="FAILED", process_wall_seconds=None,
            numerical_valid=False, environment_valid=False,
        ),
    ]
    result = aggregate.build_process_wall_aggregate(records, aggregate_id="agg-2")
    assert result["run_ids"] == ["r0"]
    assert result["sample_count"] == 1


def test_build_process_wall_aggregate_all_failed_raises():
    records = [
        _base_record(run_id="r0", run_status="FAILED", process_wall_seconds=None,
                      numerical_valid=False, environment_valid=False),
    ]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_process_wall_aggregate(records, aggregate_id="agg-3")


def test_build_process_wall_aggregate_rejects_heterogeneous_cell():
    records = [
        _base_record(run_id="r0", omp_threads=4),
        _base_record(run_id="r1", omp_threads=8),
    ]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_process_wall_aggregate(records, aggregate_id="agg-4")


def test_build_process_wall_aggregate_rejects_mixed_nstep():
    records = [
        _base_record(run_id="r0", nstep=1000),
        _base_record(run_id="r1", nstep=2000),
    ]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_process_wall_aggregate(records, aggregate_id="agg-5")


def test_build_process_wall_aggregate_quality_flag_blocks_authoritative():
    records = [
        _base_record(run_id="r0", process_wall_seconds=1.0),
        _base_record(
            run_id="r1", process_wall_seconds=1.1,
            environment_valid=False, quality_flags=["gpu_busy"],
        ),
    ]
    result = aggregate.build_process_wall_aggregate(records, aggregate_id="agg-6")
    assert result["quality_flags_union"] == ["gpu_busy"]
    assert result["all_samples_environment_valid"] is False
    assert result["authoritative"] is False
    # Every individual sample is still preserved in run_ids -- never just the fastest.
    assert result["run_ids"] == ["r0", "r1"]


_CELL = {
    "campaign_id": "camp1", "case_id": "B01_bccFe", "variant_id": "T0", "size_id": "16x16x16",
    "machine_id": "m1", "build_id": "b1", "backend": "CPU", "gpu_backend": None,
    "omp_threads": 4, "requested_precision": "DOUBLE", "measurement_profile": "DYNAMICS_ONLY",
}


def _fit(steady_step_seconds, intercept=0.01):
    return steady_state.FitResult(
        setup_or_fixed_intercept=intercept,
        steady_step_seconds=steady_step_seconds,
        fit_method=steady_state.MULTI_NSTEP_LEAST_SQUARES,
        nstep_fit_points=[1000, 2000, 3000],
    )


def test_build_fit_aggregate_steady_step_seconds():
    fits = [_fit(0.01), _fit(0.011), _fit(0.0095)]
    result = aggregate.build_fit_aggregate(
        fits, metric="steady_step_seconds", run_ids=["rep-0", "rep-1", "rep-2"],
        cell=_CELL, aggregate_id="fit-agg-1", all_numerical_valid=True, all_environment_valid=True,
    )
    assert result["metric"] == "steady_step_seconds"
    assert result["sample_count"] == 3
    assert result["fit_method"] == steady_state.MULTI_NSTEP_LEAST_SQUARES
    assert result["nstep_fit_points"] == [1000, 2000, 3000]
    assert result["nstep"] is None
    assert result["authoritative"] is True


def test_build_fit_aggregate_setup_seconds_allows_nonpositive():
    fits = [_fit(0.01, intercept=-0.002), _fit(0.011, intercept=0.0), _fit(0.0095, intercept=-0.001)]
    result = aggregate.build_fit_aggregate(
        fits, metric="setup_seconds", run_ids=["a", "b", "c"],
        cell=_CELL, aggregate_id="fit-agg-2", all_numerical_valid=True, all_environment_valid=True,
    )
    assert result["metric"] == "setup_seconds"
    assert result["median"] <= 0.0


def test_build_fit_aggregate_rejects_mismatched_fit_method():
    two_point = steady_state.FitResult(
        setup_or_fixed_intercept=0.0, steady_step_seconds=0.01,
        fit_method=steady_state.TWO_POINT_SUBTRACTION, nstep_fit_points=[1000, 2000, 3000],
    )
    fits = [_fit(0.01), two_point]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            fits, metric="steady_step_seconds", run_ids=["a", "b"],
            cell=_CELL, aggregate_id="fit-agg-3", all_numerical_valid=True, all_environment_valid=True,
        )


def test_build_fit_aggregate_rejects_mismatched_nstep_fit_points():
    other = steady_state.FitResult(
        setup_or_fixed_intercept=0.0, steady_step_seconds=0.01,
        fit_method=steady_state.MULTI_NSTEP_LEAST_SQUARES, nstep_fit_points=[500, 1500, 2500],
    )
    fits = [_fit(0.01), other]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            fits, metric="steady_step_seconds", run_ids=["a", "b"],
            cell=_CELL, aggregate_id="fit-agg-4", all_numerical_valid=True, all_environment_valid=True,
        )


def test_build_fit_aggregate_rejects_empty_run_ids():
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            [_fit(0.01)], metric="steady_step_seconds", run_ids=[],
            cell=_CELL, aggregate_id="fit-agg-5", all_numerical_valid=True, all_environment_valid=True,
        )


def test_build_fit_aggregate_requires_one_run_id_per_fit():
    fits = [_fit(0.01), _fit(0.011)]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            fits, metric="steady_step_seconds", run_ids=["only-one"],
            cell=_CELL, aggregate_id="fit-agg-5b", all_numerical_valid=True, all_environment_valid=True,
        )


def test_build_fit_aggregate_rejects_duplicate_run_ids():
    fits = [_fit(0.01), _fit(0.011)]
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            fits, metric="steady_step_seconds", run_ids=["same", "same"],
            cell=_CELL, aggregate_id="fit-agg-5c", all_numerical_valid=True, all_environment_valid=True,
        )


def test_build_fit_aggregate_rejects_unknown_metric():
    with pytest.raises(aggregate.AggregateError):
        aggregate.build_fit_aggregate(
            [_fit(0.01)], metric="process_wall_seconds", run_ids=["a"],
            cell=_CELL, aggregate_id="fit-agg-6", all_numerical_valid=True, all_environment_valid=True,
        )


def test_aggregate_never_mutates_input_records():
    records = [_base_record(run_id="r0"), _base_record(run_id="r1")]
    snapshot = copy.deepcopy(records)
    aggregate.build_process_wall_aggregate(records, aggregate_id="agg-7")
    assert records == snapshot
