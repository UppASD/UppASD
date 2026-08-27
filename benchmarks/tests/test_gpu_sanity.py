"""Tests for the WP-06 T=0 CPU/GPU sanity check (section G).

Restart-file text mirrors `test_omp_sanity.py`'s fixture -- see that module's
docstring for the real production format this represents.
"""

from __future__ import annotations

import pytest

from harness import gpu_sanity
from harness import omp_sanity

_HEADER = (
    "################################################################################\n"
    "# File type: R\n"
    "# Simulation type: S\n"
    "# Number of atoms:         2\n"
    "# Number of ensembles:         1\n"
    "################################################################################\n"
    "  # iter     ens   iatom           |Mom|             M_x             M_y             M_z\n"
)


def _restart_text(rows):
    lines = [_HEADER]
    for iteration, ens, iatom, mmom, mx, my, mz in rows:
        lines.append(f"   {iteration}       {ens}       {iatom}    {mmom:.8E}  {mx:.8E}  {my:.8E}  {mz:.8E}\n")
    return "".join(lines)


def _write_restart(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text(_restart_text(rows))
    return path


_BASE_ROWS = [
    (1000, 1, 1, 1.3, 0.0, 0.0, 1.0),
    (1000, 1, 2, 2.5, 0.0, 0.0, 1.0),
]


def _record(**overrides):
    base = {
        "run_id": "r1", "run_status": "COMPLETED", "temperature": 0.0,
        "case_id": "B01_bccFe", "variant_id": "T0", "size_id": "16x16x16",
        "backend": "CPU", "effective_gpu_precision": None,
    }
    base.update(overrides)
    return base


def _cpu_record(**overrides):
    return _record(backend="CPU", effective_gpu_precision=None, **overrides)


def _gpu_record(effective_gpu_precision="DOUBLE", **overrides):
    return _record(backend="GPU", effective_gpu_precision=effective_gpu_precision, **overrides)


# ---------------------------------------------------------------------------
# default_magnitude_rtol
# ---------------------------------------------------------------------------

def test_default_rtol_is_looser_for_single_precision():
    assert gpu_sanity.default_magnitude_rtol("SINGLE") == gpu_sanity.DEFAULT_MAGNITUDE_RTOL_SINGLE
    assert gpu_sanity.default_magnitude_rtol("DOUBLE") == gpu_sanity.DEFAULT_MAGNITUDE_RTOL_DOUBLE
    assert gpu_sanity.DEFAULT_MAGNITUDE_RTOL_SINGLE > gpu_sanity.DEFAULT_MAGNITUDE_RTOL_DOUBLE


# ---------------------------------------------------------------------------
# sanity_check_cpu_vs_gpu -- preconditions
# ---------------------------------------------------------------------------

def test_rejects_non_completed_run(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    with pytest.raises(gpu_sanity.GpuSanityError, match="COMPLETED"):
        gpu_sanity.sanity_check_cpu_vs_gpu(
            _cpu_record(run_status="FAILED"), path_a, _gpu_record(), path_b
        )


def test_rejects_finite_temperature(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    with pytest.raises(gpu_sanity.GpuSanityError, match="T=0"):
        gpu_sanity.sanity_check_cpu_vs_gpu(
            _cpu_record(temperature=300.0), path_a, _gpu_record(temperature=300.0), path_b
        )


def test_rejects_wrong_backend_for_cpu_record(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    with pytest.raises(gpu_sanity.GpuSanityError, match="backend=CPU"):
        gpu_sanity.sanity_check_cpu_vs_gpu(_gpu_record(), path_a, _gpu_record(), path_b)


def test_rejects_wrong_backend_for_gpu_record(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    with pytest.raises(gpu_sanity.GpuSanityError, match="backend=GPU"):
        gpu_sanity.sanity_check_cpu_vs_gpu(_cpu_record(), path_a, _cpu_record(), path_b)


def test_rejects_mismatched_cell(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    with pytest.raises(gpu_sanity.GpuSanityError, match="same cell"):
        gpu_sanity.sanity_check_cpu_vs_gpu(
            _cpu_record(size_id="16x16x16"), path_a, _gpu_record(size_id="32x32x32"), path_b
        )


# ---------------------------------------------------------------------------
# sanity_check_cpu_vs_gpu -- comparison behaviour
# ---------------------------------------------------------------------------

def test_identical_states_are_consistent(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    result = gpu_sanity.sanity_check_cpu_vs_gpu(_cpu_record(), path_a, _gpu_record(), path_b)
    assert result.consistent is True


def test_single_precision_gpu_uses_the_looser_default_tolerance(tmp_path):
    # A magnitude difference that would fail the DOUBLE-vs-DOUBLE tolerance
    # (1e-3) but should pass the SINGLE-precision GPU default (1e-2).
    perturbed = [
        (1000, 1, 1, 1.3 * 1.005, 0.0, 0.0, 1.0),
        (1000, 1, 2, 2.5, 0.0, 0.0, 1.0),
    ]
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", perturbed)

    with_double_default = gpu_sanity.sanity_check_cpu_vs_gpu(
        _cpu_record(), path_a, _gpu_record(effective_gpu_precision="DOUBLE"), path_b
    )
    assert with_double_default.consistent is False

    with_single_default = gpu_sanity.sanity_check_cpu_vs_gpu(
        _cpu_record(), path_a, _gpu_record(effective_gpu_precision="SINGLE"), path_b
    )
    assert with_single_default.consistent is True


def test_gross_divergence_still_fails_even_with_single_precision_tolerance(tmp_path):
    diverged = [
        (1000, 1, 1, 1.3, 0.0, 0.0, 1.0),
        (1000, 1, 2, 9.9, 0.0, 0.0, 1.0),
    ]
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", diverged)
    result = gpu_sanity.sanity_check_cpu_vs_gpu(
        _cpu_record(), path_a, _gpu_record(effective_gpu_precision="SINGLE"), path_b
    )
    assert result.consistent is False


def test_explicit_tolerance_overrides_the_precision_default(tmp_path):
    path_a = _write_restart(tmp_path, "restart.cpu.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.gpu.out", _BASE_ROWS)
    result = gpu_sanity.sanity_check_cpu_vs_gpu(
        _cpu_record(), path_a, _gpu_record(), path_b, magnitude_rtol=1e-9
    )
    assert result.consistent is True


def test_gpu_sanity_error_is_an_omp_sanity_error():
    assert issubclass(gpu_sanity.GpuSanityError, omp_sanity.OmpSanityError)
