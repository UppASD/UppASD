"""Tests for the WP-05 T=0 OpenMP thread-count sanity check.

Restart-file text mirrors the real production format written by
`source/System/restart.f90::prn_mag_conf_iter` (see a worked real example
under any `build*/**/restart.*.out` in this repository) -- never a
fabricated format of its own.
"""

from __future__ import annotations

import pytest

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
        lines.append(
            f"   {iteration}       {ens}       {iatom}    {mmom:.8E}  {mx:.8E}  {my:.8E}  {mz:.8E}\n"
        )
    return "".join(lines)


def _write_restart(tmp_path, name, rows):
    path = tmp_path / name
    path.write_text(_restart_text(rows))
    return path


_BASE_ROWS = [
    (1000, 1, 1, 1.3, 0.0, 0.0, 1.0),
    (1000, 1, 2, 2.5, 0.0, 0.0, 1.0),
]


# ---------------------------------------------------------------------------
# parse_restart_moments
# ---------------------------------------------------------------------------

def test_parse_restart_moments_reads_final_iteration(tmp_path):
    rows = [
        (500, 1, 1, 1.3, 0.1, 0.0, 0.99),
        (500, 1, 2, 2.5, 0.1, 0.0, 0.99),
        (1000, 1, 1, 1.3, 0.0, 0.0, 1.0),
        (1000, 1, 2, 2.5, 0.0, 0.0, 1.0),
    ]
    path = _write_restart(tmp_path, "restart.mixed.out", rows)
    final, natom = omp_sanity.parse_restart_moments(path)
    assert natom == 2
    assert set(final) == {(1, 1), (1, 2)}
    assert final[(1, 1)].iteration == 1000
    assert final[(1, 1)].mmom == pytest.approx(1.3)


def test_parse_restart_moments_rejects_file_without_header(tmp_path):
    path = tmp_path / "restart.bad.out"
    path.write_text("not a restart file\n")
    with pytest.raises(omp_sanity.OmpSanityError):
        omp_sanity.parse_restart_moments(path)


def test_parse_restart_moments_rejects_file_without_data_rows(tmp_path):
    path = tmp_path / "restart.empty.out"
    path.write_text(_HEADER)
    with pytest.raises(omp_sanity.OmpSanityError):
        omp_sanity.parse_restart_moments(path)


# ---------------------------------------------------------------------------
# compare_restart_moments
# ---------------------------------------------------------------------------

def test_identical_states_are_consistent(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS)
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is True
    assert result.max_magnitude_rel_diff == pytest.approx(0.0)
    assert result.reasons == []


def test_small_reduction_order_style_perturbation_is_still_consistent(tmp_path):
    perturbed = [
        (1000, 1, 1, 1.3000001, 1e-9, 0.0, 0.9999999),
        (1000, 1, 2, 2.4999998, 1e-9, 0.0, 0.9999999),
    ]
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", perturbed)
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is True


def test_gross_magnitude_divergence_is_not_consistent(tmp_path):
    diverged = [
        (1000, 1, 1, 1.3, 0.0, 0.0, 1.0),
        (1000, 1, 2, 9.9, 0.0, 0.0, 1.0),  # wildly different magnitude
    ]
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", diverged)
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is False
    assert any("magnitude" in reason for reason in result.reasons)


def test_flipped_moment_direction_is_not_consistent(tmp_path):
    flipped = [
        (1000, 1, 1, 1.3, 0.0, 0.0, -1.0),  # direction reversed
        (1000, 1, 2, 2.5, 0.0, 0.0, 1.0),
    ]
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", flipped)
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is False
    assert any("cosine" in reason for reason in result.reasons)


def test_mismatched_atom_set_is_flagged(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS[:1])
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is False
    assert any("cover the same" in reason for reason in result.reasons)


def test_mismatched_declared_atom_count_is_flagged(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = tmp_path / "restart.b.out"
    path_b.write_text(
        _restart_text(_BASE_ROWS).replace("Number of atoms:         2", "Number of atoms:         3")
    )
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is False
    assert any("atom count" in reason for reason in result.reasons)


def test_nan_moment_is_flagged(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = tmp_path / "restart.b.out"
    path_b.write_text(
        _HEADER
        + "   1000       1       1    NaN  0.0E+00  0.0E+00  1.0E+00\n"
        + "   1000       1       2    2.5000000E+00  0.0000000E+00  0.0000000E+00  1.0000000E+00\n"
    )
    result = omp_sanity.compare_restart_moments(path_a, path_b)
    assert result.consistent is False
    assert any("nan" in reason.lower() for reason in result.reasons)


# ---------------------------------------------------------------------------
# sanity_check_thread_counts -- record-level preconditions
# ---------------------------------------------------------------------------

def _record(**overrides):
    base = {
        "run_id": "r1", "run_status": "COMPLETED", "temperature": 0.0,
        "case_id": "B01_bccFe", "variant_id": "T0", "size_id": "16x16x16",
    }
    base.update(overrides)
    return base


def test_sanity_check_rejects_non_completed_run(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS)
    with pytest.raises(omp_sanity.OmpSanityError, match="COMPLETED"):
        omp_sanity.sanity_check_thread_counts(
            _record(run_status="FAILED"), path_a, _record(), path_b
        )


def test_sanity_check_rejects_finite_temperature(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS)
    with pytest.raises(omp_sanity.OmpSanityError, match="T=0"):
        omp_sanity.sanity_check_thread_counts(
            _record(temperature=300.0), path_a, _record(temperature=300.0), path_b
        )


def test_sanity_check_rejects_mismatched_cell(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS)
    with pytest.raises(omp_sanity.OmpSanityError, match="same cell"):
        omp_sanity.sanity_check_thread_counts(
            _record(size_id="16x16x16"), path_a, _record(size_id="32x32x32"), path_b
        )


def test_sanity_check_passes_through_to_comparison_for_a_legal_pair(tmp_path):
    path_a = _write_restart(tmp_path, "restart.a.out", _BASE_ROWS)
    path_b = _write_restart(tmp_path, "restart.b.out", _BASE_ROWS)
    result = omp_sanity.sanity_check_thread_counts(_record(), path_a, _record(), path_b)
    assert result.consistent is True
