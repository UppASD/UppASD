"""Throughput metrics (WP-07 sections C-D, blueprint section 12).

Three normalizations, none forced onto a workload it does not describe
(blueprint: "Do not force every workload into one normalization" / "Do not
force FFT workloads into the directed-interaction metric"):

* `spin_steps_per_second` -- `natom / t_step`. Reported for every workload.
* `directed_interaction_visits_per_second` -- `directed_interactions / t_step`,
  for `NEIGHBOR_LIST` and `NEIGHBOR_LIST_PLUS_FFT_DIPOLE` workloads only. The
  name deliberately says "visits", not "bond operations" or similar: it counts
  how many directed (i, j) neighbour entries the workload's own struct-file
  dump enumerates per second of steady-state time
  (`harness.workload_metadata.neighbor_list_from_struct_output`), not a claim
  about what a GPU kernel does internally with each one.
* `fft_grid_points_per_second` -- `fft_grid_points / t_step`, for `FFT_DIPOLE`
  and `NEIGHBOR_LIST_PLUS_FFT_DIPOLE` workloads only.

`compute_throughput` dispatches on `workload_class` and returns exactly the
metrics that apply -- a pure `NEIGHBOR_LIST` case never gets a spurious
`fft_grid_points_per_second` key, and vice versa.
"""

from __future__ import annotations

_NEIGHBOR_WORKLOAD_CLASSES = ("NEIGHBOR_LIST", "NEIGHBOR_LIST_PLUS_FFT_DIPOLE")
_FFT_WORKLOAD_CLASSES = ("FFT_DIPOLE", "NEIGHBOR_LIST_PLUS_FFT_DIPOLE")
_KNOWN_WORKLOAD_CLASSES = ("NEIGHBOR_LIST", "FFT_DIPOLE", "NEIGHBOR_LIST_PLUS_FFT_DIPOLE")


class ThroughputError(ValueError):
    """A throughput metric could not be computed from the given inputs."""


def _require_positive(value, name):
    if value is None or not value > 0:
        raise ThroughputError(f"{name} must be a positive number, got {value!r}")


def spin_steps_per_second(natom, t_step):
    """`R_spin = N_atom / t_step` (blueprint section 12). Applies to every workload."""
    _require_positive(natom, "natom")
    _require_positive(t_step, "t_step")
    return natom / t_step


def directed_interaction_visits_per_second(directed_interactions, t_step, *, workload_class):
    """`directed_interactions / t_step`, for neighbour-list workloads only.

    Raises for an `FFT_DIPOLE` workload rather than silently returning a
    meaningless number, and for `directed_interactions is None` (which the
    schema requires for `NEIGHBOR_LIST`/`NEIGHBOR_LIST_PLUS_FFT_DIPOLE`
    records, so `None` here means the caller passed the wrong metadata, not a
    legitimate zero).
    """
    if workload_class not in _NEIGHBOR_WORKLOAD_CLASSES:
        raise ThroughputError(
            f"directed_interaction_visits_per_second only applies to {_NEIGHBOR_WORKLOAD_CLASSES}, "
            f"got workload_class={workload_class!r}"
        )
    if directed_interactions is None:
        raise ThroughputError("directed_interactions is required for a neighbour-list workload")
    if directed_interactions < 0:
        raise ThroughputError(f"directed_interactions must be non-negative, got {directed_interactions!r}")
    _require_positive(t_step, "t_step")
    return directed_interactions / t_step


def fft_grid_points_per_second(fft_grid_points, t_step, *, workload_class):
    """`fft_grid_points / t_step`, for FFT/dipole workloads only."""
    if workload_class not in _FFT_WORKLOAD_CLASSES:
        raise ThroughputError(
            f"fft_grid_points_per_second only applies to {_FFT_WORKLOAD_CLASSES}, "
            f"got workload_class={workload_class!r}"
        )
    _require_positive(fft_grid_points, "fft_grid_points")
    _require_positive(t_step, "t_step")
    return fft_grid_points / t_step


def compute_throughput(
    *, natom, t_step, workload_class, directed_interactions=None, fft_grid_points=None,
):
    """The full set of throughput metrics that apply to `workload_class`.

    Always includes `spin_steps_per_second`; adds
    `directed_interaction_visits_per_second` for neighbour-list workloads and
    `fft_grid_points_per_second` for FFT/dipole workloads (both, for
    `NEIGHBOR_LIST_PLUS_FFT_DIPOLE`).
    """
    if workload_class not in _KNOWN_WORKLOAD_CLASSES:
        raise ThroughputError(f"unknown workload_class {workload_class!r}; expected one of {_KNOWN_WORKLOAD_CLASSES}")

    result = {"spin_steps_per_second": spin_steps_per_second(natom, t_step)}
    if workload_class in _NEIGHBOR_WORKLOAD_CLASSES:
        result["directed_interaction_visits_per_second"] = directed_interaction_visits_per_second(
            directed_interactions, t_step, workload_class=workload_class
        )
    if workload_class in _FFT_WORKLOAD_CLASSES:
        result["fft_grid_points_per_second"] = fft_grid_points_per_second(
            fft_grid_points, t_step, workload_class=workload_class
        )
    return result


def compute_throughput_for_aggregate(
    aggregate, *, natom, workload_class, directed_interactions=None, fft_grid_points=None,
):
    """`compute_throughput` using a `steady_step_seconds` aggregate's median as `t_step`.

    `aggregate` itself carries no `natom`/`directed_interactions`/
    `fft_grid_points` (only raw sample records do -- same caveat as
    `analysis.omp_scaling_report`), so the caller supplies them from that
    cell's own raw records or the case's size ladder.
    """
    if aggregate.get("metric") != "steady_step_seconds":
        raise ThroughputError(
            f"throughput requires a steady_step_seconds aggregate, got metric={aggregate.get('metric')!r}"
        )
    return compute_throughput(
        natom=natom, t_step=aggregate["median"], workload_class=workload_class,
        directed_interactions=directed_interactions, fft_grid_points=fft_grid_points,
    )
