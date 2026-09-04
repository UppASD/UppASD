"""Small numerical regressions for the bosonic Chern/OAM formulas.

These tests intentionally mirror only the algebra used in the Fortran
implementation.  They require no UppASD build and catch accidental changes
back to particle-only overlaps or Berry-flux-as-OAM.
"""

import cmath
import math
from pathlib import Path


def boson_overlap(v1, v2, n_a):
    return sum(a.conjugate() * b for a, b in zip(v1[:n_a], v2[:n_a])) - sum(
        a.conjugate() * b for a, b in zip(v1[n_a:], v2[n_a:])
    )


def unit_link(v1, v2, n_a):
    overlap = boson_overlap(v1, v2, n_a)
    return overlap / abs(overlap)


def paraunitarity_error(matrix, n_a):
    expected = [[0j for _ in matrix] for _ in matrix]
    for i in range(n_a):
        expected[i][i] = 1.0 + 0j
        expected[n_a + i][n_a + i] = -1.0 + 0j
    return max(
        abs(boson_overlap([row[i] for row in matrix], [row[j] for row in matrix], n_a) - expected[i][j])
        for i in range(2 * n_a)
        for j in range(2 * n_a)
    )


def test_bosonic_wilson_loop_is_phase_invariant():
    """The metric plaquette link product is unchanged by arbitrary gauges."""

    n_a = 2
    vertices = [
        [1.2 + 0.1j, 0.2 - 0.1j, 0.3 + 0.2j, -0.1 + 0.05j],
        [1.0 - 0.2j, 0.3 + 0.1j, 0.2 - 0.1j, 0.1 + 0.15j],
        [0.9 + 0.3j, 0.4 - 0.2j, 0.1 + 0.2j, 0.2 - 0.05j],
        [1.1 + 0.2j, 0.1 + 0.3j, 0.25 - 0.1j, -0.05 + 0.2j],
    ]
    # 00 -> 10 -> 11 -> 01 -> 00, matching a plaquette Wilson loop.
    loop = (
        unit_link(vertices[0], vertices[1], n_a)
        * unit_link(vertices[1], vertices[2], n_a)
        * unit_link(vertices[2], vertices[3], n_a)
        * unit_link(vertices[3], vertices[0], n_a)
    )
    phases = [cmath.exp(1j * angle) for angle in (0.31, -1.2, 2.4, -2.0)]
    gauged = [[phase * value for value in vector] for phase, vector in zip(phases, vertices)]
    gauged_loop = (
        unit_link(gauged[0], gauged[1], n_a)
        * unit_link(gauged[1], gauged[2], n_a)
        * unit_link(gauged[2], gauged[3], n_a)
        * unit_link(gauged[3], gauged[0], n_a)
    )
    assert abs(loop - gauged_loop) < 1.0e-14


def test_metric_overlap_includes_hole_sector():
    v1 = [1.0 + 0j, 0.0 + 0j, 0.4 + 0j, 0.0 + 0j]
    v2 = [0.8 + 0j, 0.0 + 0j, 0.2 + 0j, 0.0 + 0j]
    assert abs(boson_overlap(v1, v2, 2) - 0.72) < 1.0e-14
    assert abs(sum(a.conjugate() * b for a, b in zip(v1[:2], v2[:2])) - 0.8) < 1.0e-14


def test_paraunitary_squeeze_has_the_bosonic_metric():
    n_a = 2
    r = 0.37
    c = math.cosh(r)
    s = math.sinh(r)
    matrix = [
        [c, 0j, s, 0j],
        [0j, 1.0 + 0j, 0j, 0j],
        [s, 0j, c, 0j],
        [0j, 0j, 0j, 1.0 + 0j],
    ]
    assert paraunitarity_error(matrix, n_a) < 1.0e-14


def test_raw_pointwise_oam_changes_under_a_gauge_phase():
    """A phase with angular dependence changes the pointwise diagnostic."""

    n_a = 1
    mode = [math.cosh(0.4) + 0j, math.sinh(0.4) + 0j]
    qx, qy = 0.4, 0.7
    alpha = 0.9
    dtdx = [1j * alpha * qy * value for value in mode]
    dtdy = [1j * alpha * qx * value for value in mode]
    angular = [qx * y - qy * x for x, y in zip(dtdx, dtdy)]
    gauged_oam = 0.5 * boson_overlap(mode, angular, n_a).imag
    assert abs(gauged_oam) > 1.0e-3


def test_fortran_oam_is_not_the_berry_flux_proxy():
    source = (Path(__file__).parents[2] / "source" / "SpinWaves" / "chern_number.f90").read_text()
    assert "aimag( Berry_cuv ) * dkx * dky" not in source
    assert "Lz_band" not in source


def test_oblique_directional_derivatives_recover_cartesian_derivative():
    """The reciprocal-basis inverse Jacobian used by OAM is correct."""

    dq1 = (0.4, 0.2)
    dq2 = (0.1, 0.5)
    derivative_x = 1.7
    derivative_y = -0.8
    d1 = dq1[0] * derivative_x + dq1[1] * derivative_y
    d2 = dq2[0] * derivative_x + dq2[1] * derivative_y
    detq = dq1[0] * dq2[1] - dq1[1] * dq2[0]
    recovered_x = (dq2[1] * d1 - dq1[1] * d2) / detq
    recovered_y = (-dq2[0] * d1 + dq1[0] * d2) / detq
    assert math.isclose(recovered_x, derivative_x, rel_tol=0.0, abs_tol=1.0e-14)
    assert math.isclose(recovered_y, derivative_y, rel_tol=0.0, abs_tol=1.0e-14)


if __name__ == "__main__":
    test_bosonic_wilson_loop_is_phase_invariant()
    test_metric_overlap_includes_hole_sector()
    test_paraunitary_squeeze_has_the_bosonic_metric()
    test_raw_pointwise_oam_changes_under_a_gauge_phase()
    test_fortran_oam_is_not_the_berry_flux_proxy()
    test_oblique_directional_derivatives_recover_cartesian_derivative()
    print("Fishman OAM algebra regressions passed")
