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


def ring_oam(eigenvectors):
    """Return the pointwise Fishman OAM for a periodic ring."""

    n_phi = len(eigenvectors)
    dphi = 2.0 * math.pi / n_phi
    values = []
    for j, vector in enumerate(eigenvectors):
        plus = eigenvectors[(j + 1) % n_phi]
        minus = eigenvectors[(j - 1) % n_phi]
        derivative = [(a - b) / (2.0 * dphi) for a, b in zip(plus, minus)]
        values.append(-0.5 * boson_overlap(vector, derivative, 1).imag)
    return values


def parallel_transport_ring(raw_vectors):
    """Mirror the production spanning-ring gauge for one isolated band."""

    n_phi = len(raw_vectors)
    fixed = [list(raw_vectors[0])]
    for j in range(1, n_phi):
        current = list(raw_vectors[j])
        overlap = boson_overlap(fixed[-1], current, 1)
        assert abs(overlap) > 1.0e-8
        phase = overlap.conjugate() / abs(overlap)
        fixed.append([phase * value for value in current])
    closure = boson_overlap(fixed[-1], fixed[0], 1)
    theta = cmath.phase(closure)
    return [
        [value * cmath.exp(1j * j * theta / n_phi) for value in vector]
        for j, vector in enumerate(fixed)
    ]


def test_fishman_minus_sign_and_periodic_gauge_invariance():
    """The corrected T=X^-1 sign and periodic-gauge average are checked."""

    n_phi = 256
    squeeze = 0.4
    mode = [math.cosh(squeeze) + 0j, math.sinh(squeeze) + 0j]
    winding_mode = [
        [mode[0], mode[1] * cmath.exp(1j * 1.0 * 2.0 * math.pi * j / n_phi)]
        for j in range(n_phi)
    ]
    pointwise = ring_oam(winding_mode)
    expected = 0.5 * math.sinh(squeeze) ** 2
    assert math.isclose(sum(pointwise) / n_phi, expected, rel_tol=0.0, abs_tol=5.0e-5)

    phases = [
        0.7 * math.sin(2.0 * math.pi * j / n_phi)
        + 0.2 * math.cos(4.0 * math.pi * j / n_phi)
        + 0.1 * math.sin(6.0 * math.pi * j / n_phi)
        for j in range(n_phi)
    ]
    gauged = [
        [value * cmath.exp(-1j * phase) for value in vector]
        for vector, phase in zip(winding_mode, phases)
    ]
    gauged_pointwise = ring_oam(gauged)
    assert max(abs(a - b) for a, b in zip(pointwise, gauged_pointwise)) > 1.0e-3
    phase_derivative = [
        (phases[(j + 1) % n_phi] - phases[(j - 1) % n_phi])
        / (4.0 * math.pi / n_phi)
        for j in range(n_phi)
    ]
    # The phase convention above is exp(-i lambda): O' - O = +1/2 d(lambda)
    # for the hole weight of this squeezed mode.
    expected_shift = [0.5 * value for value in phase_derivative]
    assert max(
        abs((after - before) - shift)
        for before, after, shift in zip(pointwise, gauged_pointwise, expected_shift)
    ) < 2.0e-4
    # lambda=m*phi is intentionally not used here: it is not single-valued
    # on the ring and changes the allowed absolute Fishman branch.
    assert abs(sum(pointwise) / n_phi - sum(gauged_pointwise) / n_phi) < 2.0e-5


def test_random_raw_phases_are_removed_by_ring_gauge():
    n_phi = 128
    squeeze = 0.31
    raw = [
        [
            math.cosh(squeeze) + 0j,
            math.sinh(squeeze) * cmath.exp(1j * 2.0 * math.pi * j / n_phi),
        ]
        for j in range(n_phi)
    ]
    reference = sum(ring_oam(parallel_transport_ring(raw))) / n_phi
    randomized = [
        [value * cmath.exp(1j * phase) for value in vector]
        for vector, phase in zip(raw, [0.1 * (j + 1) ** 2 for j in range(n_phi)])
    ]
    recovered = sum(ring_oam(parallel_transport_ring(randomized))) / n_phi
    assert math.isclose(reference, recovered, rel_tol=0.0, abs_tol=2.0e-5)


def test_fortran_oam_is_not_the_berry_flux_proxy():
    source = (Path(__file__).parents[2] / "source" / "SpinWaves" / "chern_number.f90").read_text()
    assert "aimag( Berry_cuv ) * dkx * dky" not in source
    assert "Lz_band" not in source
    assert "f_oam(i)=f_oam(i)-0.5_dblprec*aimag" in source
    assert "f_oam_nphi must be an even integer >= 8" in source


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


def reciprocal_basis(c1, c2, c3):
    r1 = (
        c2[1] * c3[2] - c2[2] * c3[1],
        c2[2] * c3[0] - c2[0] * c3[2],
        c2[0] * c3[1] - c2[1] * c3[0],
    )
    volume = sum(a * b for a, b in zip(c1, r1))
    r2 = (
        c3[1] * c1[2] - c3[2] * c1[1],
        c3[2] * c1[0] - c3[0] * c1[2],
        c3[0] * c1[1] - c3[1] * c1[0],
    )
    r3 = (
        c1[1] * c2[2] - c1[2] * c2[1],
        c1[2] * c2[0] - c1[0] * c2[2],
        c1[0] * c2[1] - c1[1] * c2[0],
    )
    return [tuple(x / volume for x in r) for r in (r1, r2, r3)]


def test_cartesian_reduced_round_trip_for_oblique_cell():
    c1 = (2.0, 0.0, 0.0)
    c2 = (1.0, math.sqrt(3.0), 0.0)
    c3 = (0.0, 0.0, 1.0)
    b1, b2, b3 = reciprocal_basis(c1, c2, c3)
    basis = (b1, b2, b3)
    k = (0.37, -0.51, 0.0)
    q = tuple(sum(k[i] * c[i] for i in range(3)) / (2.0 * math.pi) for c in (c1, c2, c3))
    recovered = tuple(
        2.0 * math.pi * sum(q[j] * basis[j][i] for j in range(3))
        for i in range(3)
    )
    assert max(abs(a - b) for a, b in zip(k, recovered)) < 1.0e-14


if __name__ == "__main__":
    test_bosonic_wilson_loop_is_phase_invariant()
    test_metric_overlap_includes_hole_sector()
    test_paraunitary_squeeze_has_the_bosonic_metric()
    test_raw_pointwise_oam_changes_under_a_gauge_phase()
    test_fishman_minus_sign_and_periodic_gauge_invariance()
    test_random_raw_phases_are_removed_by_ring_gauge()
    test_fortran_oam_is_not_the_berry_flux_proxy()
    test_oblique_directional_derivatives_recover_cartesian_derivative()
    test_cartesian_reduced_round_trip_for_oblique_cell()
    print("Fishman OAM algebra regressions passed")
