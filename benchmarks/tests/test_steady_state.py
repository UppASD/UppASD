"""Unit tests for the WP-03 steady-state extraction (sections B, C, D).

Pure logic: no executable is invoked. `calibrate_step_span` is exercised
against an injected synthetic `measure_fn`.
"""

import pytest

from harness import steady_state


# ---------------------------------------------------------------------------
# B. Multi-nstep regression
# ---------------------------------------------------------------------------

def test_multi_nstep_fit_recovers_exact_line():
    t_fixed, t_step = 0.25, 0.002
    samples = [(n, t_fixed + n * t_step) for n in (1_000, 5_000, 20_000)]
    fit = steady_state.fit_multi_nstep(samples)
    assert fit.setup_or_fixed_intercept == pytest.approx(t_fixed, abs=1e-9)
    assert fit.steady_step_seconds == pytest.approx(t_step, abs=1e-12)
    assert fit.fit_method == steady_state.MULTI_NSTEP_LEAST_SQUARES
    assert fit.nstep_fit_points == [1_000, 5_000, 20_000]


def test_multi_nstep_fit_handles_noisy_repeated_samples():
    t_fixed, t_step = 0.1, 0.0015
    samples = []
    for n in (2_000, 8_000, 32_000):
        for noise in (-0.001, 0.0, 0.001):
            samples.append((n, t_fixed + n * t_step + noise))
    fit = steady_state.fit_multi_nstep(samples)
    assert fit.steady_step_seconds == pytest.approx(t_step, abs=1e-4)
    assert fit.nstep_fit_points == [2_000, 8_000, 32_000]


def test_multi_nstep_fit_requires_at_least_three_distinct_nstep():
    samples = [(1_000, 1.0), (2_000, 1.5)]
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.fit_multi_nstep(samples)


def test_multi_nstep_fit_rejects_a_flat_line():
    samples = [(1_000, 1.0), (2_000, 1.0), (3_000, 1.0)]
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.fit_multi_nstep(samples)


def test_multi_nstep_fit_rejects_a_negative_slope():
    samples = [(1_000, 3.0), (2_000, 2.0), (3_000, 1.0)]
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.fit_multi_nstep(samples)


# ---------------------------------------------------------------------------
# C. Lean two-point mode
# ---------------------------------------------------------------------------

def test_two_point_estimate_matches_the_exact_slope():
    fit = steady_state.two_point_estimate(1_000, 1.0, 5_000, 5.0)
    assert fit.steady_step_seconds == pytest.approx(0.001)
    assert fit.setup_or_fixed_intercept == pytest.approx(0.0)
    assert fit.fit_method == steady_state.TWO_POINT_SUBTRACTION
    assert fit.nstep_fit_points == [1_000, 5_000]


def test_two_point_estimate_is_order_independent():
    forward = steady_state.two_point_estimate(1_000, 1.0, 5_000, 5.0)
    backward = steady_state.two_point_estimate(5_000, 5.0, 1_000, 1.0)
    assert forward.steady_step_seconds == pytest.approx(backward.steady_step_seconds)
    assert forward.setup_or_fixed_intercept == pytest.approx(backward.setup_or_fixed_intercept)


def test_two_point_estimate_rejects_equal_nstep():
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.two_point_estimate(1_000, 1.0, 1_000, 1.2)


def test_two_point_estimate_rejects_a_non_positive_slope():
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.two_point_estimate(1_000, 5.0, 5_000, 5.0)


# ---------------------------------------------------------------------------
# D. Calibration / pilot step
# ---------------------------------------------------------------------------

def test_calibration_accepts_an_already_resolvable_span():
    t_fixed, t_step = 0.1, 0.01

    def measure(n):
        return t_fixed + n * t_step

    result = steady_state.calibrate_step_span(measure, n1=100, n2=200, max_nstep=10_000)
    assert result.n1 == 100
    assert result.n2 == 200
    assert result.iterations == 0
    assert result.jitter_estimate == pytest.approx(0.05 * measure(100))  # single pilot sample fallback


def test_calibration_widens_the_span_until_resolvable():
    t_fixed, t_step = 0.1, 1e-7  # tiny slope: n1/n2 start too close together to resolve
    calls = []

    def measure(n):
        calls.append(n)
        return t_fixed + n * t_step

    result = steady_state.calibrate_step_span(
        measure, n1=10, n2=20, max_nstep=10_000_000, min_separation_ratio=5.0,
        growth_factor=4.0, max_iterations=20,
    )
    assert result.n2 > 20
    assert result.iterations > 0
    # The widened span must actually be resolvable against the measured jitter.
    t1 = min(result.samples[result.n1])
    t2 = min(result.samples[result.n2])
    assert abs(t2 - t1) >= 5.0 * result.jitter_estimate


def test_calibration_raises_when_max_nstep_cannot_resolve_the_span():
    def measure(_n):
        return 1.0  # perfectly flat: no span will ever resolve

    with pytest.raises(steady_state.SteadyStateError):
        steady_state.calibrate_step_span(measure, n1=10, n2=20, max_nstep=1_000, max_iterations=4)


def test_calibration_requires_n2_greater_than_n1():
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.calibrate_step_span(lambda n: 1.0, n1=100, n2=100, max_nstep=1_000)


def test_estimate_jitter_uses_median_absolute_deviation():
    jitter = steady_state.estimate_jitter([1.0, 1.02, 0.98, 1.05, 0.95])
    assert jitter == pytest.approx(0.02, abs=1e-9)


def test_estimate_jitter_falls_back_for_a_single_sample():
    assert steady_state.estimate_jitter([2.0]) == pytest.approx(0.1)


def test_estimate_jitter_requires_at_least_one_sample():
    with pytest.raises(steady_state.SteadyStateError):
        steady_state.estimate_jitter([])
