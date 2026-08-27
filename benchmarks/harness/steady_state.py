"""Steady-state timestep extraction (WP-03 sections B, C, D).

Implements the production timing model

    T(n) = T_fixed + n * t_step

from complete production executions at several `Nstep` values, plus a lean
two-point subtraction and a pilot/calibration step that widens the sampled
step-count spread when it is too small relative to run-to-run jitter.

Never derives `t_step` from a profiler phase timing or a kernel
microbenchmark: every point fed in here is a `process_wall_seconds` from one
complete production run (see `runner.py`).
"""

from __future__ import annotations

import statistics

MULTI_NSTEP_LEAST_SQUARES = "MULTI_NSTEP_LEAST_SQUARES"
TWO_POINT_SUBTRACTION = "TWO_POINT_SUBTRACTION"

MIN_DISTINCT_NSTEP_FOR_REGRESSION = 3


class SteadyStateError(ValueError):
    """A steady-state fit could not be computed from the given samples."""


class FitResult:
    """A fitted ``T(n) = T_fixed + n*t_step`` model.

    `setup_or_fixed_intercept` is deliberately named for what it is -- a
    fitted intercept -- not assumed to be pure startup time (blueprint
    section 6.2). `timing_method` on the samples that fed the fit remains
    `EXTERNAL_PROCESS_WALL_CLOCK`/`INTERNAL_REPORTED`; `fit_method` here is
    the separate, aggregate-level provenance of the derived quantities.
    """

    def __init__(self, setup_or_fixed_intercept, steady_step_seconds, fit_method, nstep_fit_points):
        self.setup_or_fixed_intercept = setup_or_fixed_intercept
        self.steady_step_seconds = steady_step_seconds
        self.fit_method = fit_method
        self.nstep_fit_points = nstep_fit_points

    def __repr__(self):
        return (
            f"FitResult(setup_or_fixed_intercept={self.setup_or_fixed_intercept!r}, "
            f"steady_step_seconds={self.steady_step_seconds!r}, "
            f"fit_method={self.fit_method!r}, nstep_fit_points={self.nstep_fit_points!r})"
        )


def _ordinary_least_squares(xs, ys):
    """Return (intercept, slope) minimizing sum((y - (intercept + slope*x))**2)."""
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx == 0.0:
        raise SteadyStateError("all samples share one nstep; cannot fit a slope")
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    return intercept, slope


def fit_multi_nstep(samples, *, min_distinct_nstep=MIN_DISTINCT_NSTEP_FOR_REGRESSION):
    """Fit ``T(n) = T_fixed + n*t_step`` by least squares over several runs.

    ``samples`` is an iterable of ``(nstep, process_wall_seconds)`` pairs
    from individually valid, completed production runs. Repeated samples at
    the same `nstep` are accepted and simply weight that point in the fit.
    Raises :class:`SteadyStateError` if fewer than ``min_distinct_nstep``
    distinct step counts are present -- authoritative campaigns require at
    least 3 (blueprint section 6.3) -- or if the fitted slope is not
    positive, which real production timesteps must never yield.
    """
    samples = list(samples)
    distinct = sorted({nstep for nstep, _ in samples})
    if len(distinct) < min_distinct_nstep:
        raise SteadyStateError(
            f"multi-nstep fit needs at least {min_distinct_nstep} distinct nstep "
            f"values, got {len(distinct)}: {distinct}"
        )
    xs = [float(nstep) for nstep, _ in samples]
    ys = [float(seconds) for _, seconds in samples]
    intercept, slope = _ordinary_least_squares(xs, ys)
    if slope <= 0.0:
        raise SteadyStateError(
            f"fitted steady_step_seconds is not positive ({slope!r}); the "
            "samples do not support a valid steady-state estimate"
        )
    return FitResult(
        setup_or_fixed_intercept=intercept,
        steady_step_seconds=slope,
        fit_method=MULTI_NSTEP_LEAST_SQUARES,
        nstep_fit_points=distinct,
    )


def two_point_estimate(n1, t1, n2, t2):
    """Lean developer-tier estimate: ``t_step = (T(n2) - T(n1)) / (n2 - n1)``.

    Labelled ``timing_method = TWO_POINT_SUBTRACTION`` rather than
    regression (blueprint section 6.3/C). Both points must come from
    complete production runs, exactly like the regression path -- this
    function only ever does the subtraction, never the execution.
    """
    if n1 == n2:
        raise SteadyStateError(f"two-point estimate needs two distinct nstep values, got {n1!r} twice")
    n1, n2, t1, t2 = float(n1), float(n2), float(t1), float(t2)
    if n2 < n1:
        n1, n2, t1, t2 = n2, n1, t2, t1
    slope = (t2 - t1) / (n2 - n1)
    if slope <= 0.0:
        raise SteadyStateError(
            f"two-point steady_step_seconds is not positive ({slope!r}) for "
            f"n1={n1}, t1={t1}, n2={n2}, t2={t2}"
        )
    intercept = t1 - n1 * slope
    return FitResult(
        setup_or_fixed_intercept=intercept,
        steady_step_seconds=slope,
        fit_method=TWO_POINT_SUBTRACTION,
        nstep_fit_points=[int(n1), int(n2)],
    )


# ---------------------------------------------------------------------------
# D. Calibration / pilot step
# ---------------------------------------------------------------------------

class CalibrationResult:
    """Outcome of widening the sampled step-count spread until it is
    resolvable above run-to-run jitter."""

    def __init__(self, n1, n2, samples, jitter_estimate, iterations):
        self.n1 = n1
        self.n2 = n2
        self.samples = samples  # {nstep: [process_wall_seconds, ...]}
        self.jitter_estimate = jitter_estimate
        self.iterations = iterations


def estimate_jitter(repeated_samples):
    """Median absolute deviation of repeated same-``nstep`` timing samples.

    Falls back to a small fraction of the sample magnitude when only one
    pilot sample is available, so calibration still has a usable (if
    conservative) jitter floor rather than dividing by zero.
    """
    repeated_samples = list(repeated_samples)
    if len(repeated_samples) < 2:
        if not repeated_samples:
            raise SteadyStateError("need at least one pilot sample to estimate jitter")
        return 0.05 * repeated_samples[0]
    med = statistics.median(repeated_samples)
    return statistics.median(abs(s - med) for s in repeated_samples) or (0.05 * med)


def calibrate_step_span(
    measure_fn,
    *,
    n1,
    n2,
    max_nstep,
    pilot_repeats=2,
    min_separation_ratio=5.0,
    growth_factor=2.0,
    max_iterations=6,
):
    """Widen ``[n1, n2]`` until ``T(n2) - T(n1)`` clears run-to-run jitter.

    ``measure_fn(nstep)`` must return one real ``process_wall_seconds`` from
    a complete production run at that step count (``runner.run_sample`` in
    practice); this function never fabricates a timing value. Jitter is
    estimated from ``pilot_repeats`` repeated measurements at ``n1``. The
    span doubles (``growth_factor```) up to ``max_nstep`` and
    ``max_iterations`` attempts; :class:`SteadyStateError` is raised if no
    resolvable span is found within those configured bounds -- this function
    never silently accepts an unresolvable span.
    """
    if n2 <= n1:
        raise SteadyStateError(f"calibration needs n2 > n1, got n1={n1}, n2={n2}")
    if pilot_repeats < 1:
        raise SteadyStateError("pilot_repeats must be at least 1")

    samples = {n1: [measure_fn(n1) for _ in range(pilot_repeats)]}
    jitter = estimate_jitter(samples[n1])

    iterations = 0
    while True:
        if n2 not in samples:
            samples[n2] = [measure_fn(n2)]
        t1 = statistics.median(samples[n1])
        t2 = statistics.median(samples[n2])
        separation = abs(t2 - t1)
        if separation >= min_separation_ratio * jitter:
            return CalibrationResult(n1=n1, n2=n2, samples=samples, jitter_estimate=jitter, iterations=iterations)

        iterations += 1
        if iterations >= max_iterations:
            raise SteadyStateError(
                f"could not find a resolvable step-count span within "
                f"{max_iterations} iterations (n1={n1}, n2={n2}, jitter={jitter!r}, "
                f"last separation={separation!r}); widen max_nstep or investigate jitter"
            )
        next_n2 = min(max_nstep, int(n1 + (n2 - n1) * growth_factor))
        if next_n2 == n2:
            raise SteadyStateError(
                f"reached max_nstep={max_nstep} without a resolvable step-count "
                f"span (n1={n1}, n2={n2}, jitter={jitter!r}, separation={separation!r})"
            )
        n2 = next_n2
