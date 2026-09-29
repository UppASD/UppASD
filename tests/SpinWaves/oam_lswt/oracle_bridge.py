"""Independent all-site oracle for the LSWT/trajectory OAM bridge.

This module is deliberately separate from the Fortran mesh and from
``run_bridge_checks.py``.  It provides:

* a minimum-image, all-site weighted-least-squares gradient (A--B neighbours
  are included);
* the fixed-spinor packet used by the first bridge probe; and
* a k-space annular packet whose band spinor varies around the angular ring.
* a reciprocal-supercell packet with analytic spatial derivatives, suitable
  for separating lattice/WLS error from packet-construction error.

The fixed-spinor packet is useful as a control, but it must not be expected to
produce the Berry/Fishman term: it contains ``u(k0)`` but not ``d u/d phi``.
The annular packet is the development oracle for that term.  For the
particle-only HP field and UppASD's ``exp(+i k.r)`` reconstruction, the bridge
target is ``l - 2 F_n``.  The factor of two is part of the field-to-band
normalization; it is not a fitting parameter.

Examples::

    PYTHONPATH=tests/SpinWaves/oam_lswt python3 \\
        tests/SpinWaves/oam_lswt/oracle_bridge.py --packet fixed

    PYTHONPATH=tests/SpinWaves/oam_lswt python3 \\
        tests/SpinWaves/oam_lswt/oracle_bridge.py --packet ring --n 24
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np

import oracle_honey as H


S3 = np.sqrt(3.0)
C1 = np.array([1.0, 0.0])
C2 = np.array([0.5, S3 / 2.0])
BASIS = np.array([[0.0, 0.0], [0.5, 0.5 / S3]])
# The oriented dmfile written by mkhoney.py maps to the opposite algebraic
# Haldane-arrow sign in oracle_honey.py.  nu=-1 matches the confirmed positive
# UppASD-D / positive-Fishman convention used by the bridge test.
NU = -1


def build_coords(n: int) -> np.ndarray:
    """Honeycomb coordinates in UppASD cell-x-fastest, basis-inner order."""

    return np.array([
        ix * C1 + iy * C2 + basis
        for iy in range(n)
        for ix in range(n)
        for basis in BASIS
    ])


def minimum_image(displacement: np.ndarray, n: int) -> np.ndarray:
    cell = np.column_stack((n * C1, n * C2))
    reduced = np.linalg.solve(cell, displacement)
    reduced -= np.round(reduced)
    return cell @ reduced


def all_site_wls_gradient(coords: np.ndarray, psi: np.ndarray, n: int,
                          max_neighbors: int = 12) -> tuple[np.ndarray, np.ndarray]:
    """Return complex dpsi/dx,dpsi/dy using every sublattice's neighbours.

    The fit is local and linear:

        psi_j - psi_i = dx_ij * dpsi/dx + dy_ij * dpsi/dy.

    Periodic displacements use the minimum image in the oblique supercell.
    The implementation is intentionally simple and independent of the
    production FEM mesh; the oracle is run on modest validation cells.
    """

    natom = len(psi)
    if coords.shape != (natom, 2):
        raise ValueError("coords must have shape (natom, 2)")
    if max_neighbors < 2 or max_neighbors >= natom:
        raise ValueError("max_neighbors must be between 2 and natom-1")

    cell = np.column_stack((n * C1, n * C2))
    inv_cell = np.linalg.inv(cell)
    dx = coords[None, :, :] - coords[:, None, :]
    reduced = dx @ inv_cell.T
    reduced -= np.round(reduced)
    dx = reduced @ cell.T
    distance2 = np.einsum("ijk,ijk->ij", dx, dx)
    np.fill_diagonal(distance2, np.inf)

    grad_x = np.zeros(natom, dtype=complex)
    grad_y = np.zeros(natom, dtype=complex)
    for i in range(natom):
        neighbours = np.argpartition(distance2[i], max_neighbors - 1)[:max_neighbors]
        design = dx[i, neighbours]
        weights = 1.0 / np.maximum(np.einsum("ij,ij->i", design, design), 1.0e-14)
        lhs = (design.T * weights) @ design
        rhs = (design.T * weights) @ (psi[neighbours] - psi[i])
        try:
            gradient = np.linalg.solve(lhs, rhs)
        except np.linalg.LinAlgError as exc:
            raise ValueError(f"singular all-site gradient fit at site {i}") from exc
        grad_x[i], grad_y[i] = gradient
    return grad_x, grad_y


def oam_from_gradient(coords: np.ndarray, psi: np.ndarray,
                      grad_x: np.ndarray, grad_y: np.ndarray) -> float:
    """Evaluate centroid-referenced dimensionless OAM for an all-site field."""

    weight = np.abs(psi) ** 2
    norm = float(np.sum(weight))
    if norm <= 1.0e-14:
        return float("nan")
    centroid = np.sum(coords * weight[:, None], axis=0) / norm
    lever = coords - centroid
    density = np.imag(np.conjugate(psi) *
                      (lever[:, 0] * grad_y - lever[:, 1] * grad_x))
    return float(np.sum(density) / norm)


def exact_kspace_reference(k0: float, ell: int, *, D: float = 0.3,
                           band: int = 0, nphi: int = 2048) -> OracleResult:
    """Return the gauge-invariant k-space bridge reference.

    ``F_n`` is the contract quantity, obtained from the Wilson-loop phase.
    With the documented ``exp(+i k.r)`` reconstruction and
    ``lambda = Im[psi* (x d_y-y d_x) psi] / |psi|^2``, the particle-only
    Fourier-field value is
    ``ell + Im[u^dagger d_phi u] = ell - 2 F_n``.
    """

    points = [np.array([k0 * np.cos(angle), k0 * np.sin(angle)])
              for angle in 2.0 * np.pi * np.arange(nphi) / nphi]
    gamma = float(H.ring_berry_phase(points, band, J=1.0, D=D, S=1.0, nu=NU))
    fishman_f = gamma / (4.0 * np.pi)
    return OracleResult("kspace", ell, k0, gamma, fishman_f, float("nan"),
                        ell - 2.0 * fishman_f)


def ring_spinors(k0: float, nphi: int, *, D: float = 0.3,
                 band: int = 0) -> tuple[np.ndarray, np.ndarray, float]:
    """Return a continuously transported, periodic-gauge ring of spinors."""

    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    points = [np.array([k0 * np.cos(angle), k0 * np.sin(angle)]) for angle in phi]
    raw = [H.evecs(point, J=1.0, D=D, S=1.0, nu=NU)[1][:, band] for point in points]

    continuous = [raw[0]]
    for vector in raw[1:]:
        overlap = np.vdot(continuous[-1], vector)
        continuous.append(vector * np.exp(-1j * np.angle(overlap)))

    gamma = H.ring_berry_phase(points, band, J=1.0, D=D, S=1.0, nu=NU)
    # Spread the holonomy uniformly so the sampled gauge is periodic.  This
    # makes the Fourier packet reproducible; the ring Berry phase itself is
    # gauge invariant and remains the target quantity.
    periodic = np.array([
        vector * np.exp(-1j * gamma * index / nphi)
        for index, vector in enumerate(continuous)
    ])
    return phi, periodic, float(gamma)


def fixed_spinor_packet(coords: np.ndarray, k0: float, ell: int, *, D: float = 0.3,
                        sigma_r: float = 4.0, amplitude: float = 1.0) -> np.ndarray:
    """The control packet: one band spinor u(k0) times a real-space envelope."""

    _, vectors = H.evecs(np.array([k0, 0.0]), J=1.0, D=D, S=1.0, nu=NU)
    mode = vectors[:, 0]
    center = coords.mean(axis=0)
    displacement = coords - center
    theta = np.arctan2(displacement[:, 1], displacement[:, 0])
    envelope = np.exp(-0.5 * np.sum(displacement ** 2, axis=1) / sigma_r ** 2)
    phase = np.exp(1j * (k0 * coords[:, 0] + ell * theta))
    return amplitude * mode[np.arange(len(coords)) % 2] * envelope * phase


def ring_packet(coords: np.ndarray, k0: float, ell: int, *, D: float = 0.3,
                sigma_k: float = 0.12, nphi: int = 128, nrad: int = 9,
                amplitude: float = 1.0) -> np.ndarray:
    """Synthesize a narrow annular band packet in real space.

    The radial envelope is sampled around k0 and the angular band spinor is
    sampled in a periodic gauge.  This is a development reference, not a
    periodic-boundary initializer: the returned field is evaluated on a
    finite honeycomb patch and should be checked for boundary leakage.
    """

    phi, spinors, _ = ring_spinors(k0, nphi, D=D)
    radial = np.linspace(max(1.0e-4, k0 - 3.0 * sigma_k),
                         k0 + 3.0 * sigma_k, nrad)
    weights = np.exp(-0.5 * ((radial - k0) / sigma_k) ** 2) * radial
    field = np.zeros(len(coords), dtype=complex)
    for radius, radial_weight in zip(radial, weights):
        for angle, spinor in zip(phi, spinors):
            wavevector = radius * np.array([np.cos(angle), np.sin(angle)])
            field += radial_weight * spinor[np.arange(len(coords)) % 2] \
                * np.exp(1j * ell * angle) \
                * np.exp(1j * (coords @ wavevector))
    field *= amplitude / (nphi * nrad)
    return field


def periodic_band_packet(coords: np.ndarray, n: int, k0: float, ell: int, *,
                         D: float = 0.3, sigma_k: float = 0.18,
                         amplitude: float = 1.0):
    """Build a packet from reciprocal-supercell modes and exact derivatives.

    Only reciprocal vectors compatible with the periodic ``n x n`` cell are
    used.  The particle and hole terms are written separately even though the
    present Haldane FM has ``v=0``.  The returned analytic gradient is an
    oracle for the sampled field; the WLS gradient remains an independent
    approximation to compare against it.
    """

    positions = coords - coords.mean(axis=0)
    cell = np.column_stack((n * C1, n * C2))
    reciprocal = 2.0 * np.pi * np.linalg.inv(cell).T
    field = np.zeros(len(coords), dtype=complex)
    grad_x = np.zeros(len(coords), dtype=complex)
    grad_y = np.zeros(len(coords), dtype=complex)
    sublattice = np.arange(len(coords)) % 2
    modes = []
    for m1 in range(-n // 2, n // 2 + 1):
        for m2 in range(-n // 2, n // 2 + 1):
            wavevector = reciprocal @ np.array([m1, m2], dtype=float)
            radius = float(np.linalg.norm(wavevector))
            if radius <= 1.0e-12:
                continue
            radial_weight = np.exp(-0.5 * ((radius - k0) / sigma_k) ** 2)
            if radial_weight < 1.0e-7:
                continue
            angle = float(np.arctan2(wavevector[1], wavevector[0]))
            _, vectors = H.evecs(wavevector, J=1.0, D=D, S=1.0, nu=NU)
            u = vectors[:, 0]
            # A deterministic local gauge.  The exact Wilson-loop target is
            # gauge invariant; this gauge only defines the packet amplitude.
            if abs(u[0]) > 1.0e-12:
                u = u * np.exp(-1j * np.angle(u[0]))
            v = np.zeros_like(u)  # particle-conserving Haldane control
            beta = radial_weight * np.exp(1j * ell * angle)
            phase = np.exp(1j * (positions @ wavevector))
            particle = beta * u[sublattice] * phase
            hole = np.conjugate(beta) * np.conjugate(v[sublattice]) * np.conjugate(phase)
            field += particle + hole
            grad_x += 1j * wavevector[0] * particle - 1j * wavevector[0] * hole
            grad_y += 1j * wavevector[1] * particle - 1j * wavevector[1] * hole
            modes.append(radius)
    if not modes:
        raise ValueError("reciprocal packet has no modes; increase sigma_k or n")
    scale = amplitude / len(modes)
    return field * scale, grad_x * scale, grad_y * scale, len(modes)


@dataclass
class OracleResult:
    packet: str
    ell: int
    k0: float
    gamma: float
    fishman_f: float
    lambda_all_site: float
    bridge_target: float
    lambda_exact_gradient: float = float("nan")


def evaluate(packet_name: str, n: int, k0: float, ell: int, *, D: float = 0.3,
             sigma_k: float = 0.18, max_neighbors: int = 12) -> OracleResult:
    coords = build_coords(n)
    if packet_name == "fixed":
        psi = fixed_spinor_packet(coords, k0, ell, D=D)
        grad_exact_x = grad_exact_y = None
    elif packet_name == "ring":
        psi = ring_packet(coords, k0, ell, D=D)
        grad_exact_x = grad_exact_y = None
    elif packet_name == "periodic":
        psi, grad_exact_x, grad_exact_y, _ = periodic_band_packet(
            coords, n, k0, ell, D=D, sigma_k=sigma_k
        )
    else:
        raise ValueError(f"unknown packet type: {packet_name}")
    grad_x, grad_y = all_site_wls_gradient(coords, psi, n, max_neighbors)
    lam = oam_from_gradient(coords, psi, grad_x, grad_y)
    lam_exact = (oam_from_gradient(coords, psi, grad_exact_x, grad_exact_y)
                 if grad_exact_x is not None else float("nan"))
    reference = exact_kspace_reference(k0, ell, D=D)
    gamma = reference.gamma
    fishman_f = reference.fishman_f
    return OracleResult(packet_name, ell, k0, gamma, fishman_f, lam,
                        ell - 2.0 * fishman_f,
                        lambda_exact_gradient=lam_exact)


def selftest() -> None:
    coords = build_coords(8)
    constant = np.ones(len(coords), dtype=complex)
    grad_x, grad_y = all_site_wls_gradient(coords, constant, 8)
    if not (np.max(np.abs(grad_x)) == 0.0 and np.max(np.abs(grad_y)) == 0.0):
        raise AssertionError("constant field has a nonzero all-site gradient")
    if not np.isfinite(oam_from_gradient(coords, constant, grad_x, grad_y)):
        raise AssertionError("constant-field OAM is not finite")
    points = [np.array([1.7 * np.cos(angle), 1.7 * np.sin(angle)])
              for angle in 2.0 * np.pi * np.arange(128) / 128]
    vectors = [H.evecs(point, J=1.0, D=0.3, S=1.0, nu=NU)[1][:, 0]
               for point in points]
    rng = np.random.default_rng(7)
    phases = np.exp(1j * rng.uniform(0.0, 2.0 * np.pi, len(vectors)))
    phased = [vector * phase for vector, phase in zip(vectors, phases)]
    def loop_phase(ring):
        product = 1.0 + 0.0j
        for index, vector in enumerate(ring):
            product *= np.vdot(vector, ring[(index + 1) % len(ring)])
        return -np.angle(product)
    if abs(loop_phase(vectors) - loop_phase(phased)) > 1.0e-12:
        raise AssertionError("Wilson-loop phase changed under a gauge transform")
    print("oracle_bridge selftest: PASS")


def validate_particle_bridge(n: int, k0: float, *, D: float = 0.3,
                             sigma_k: float = 0.18,
                             tolerance: float = 0.03) -> bool:
    """Validate the settled particle-field bridge convention.

    The reciprocal packet has an analytic spatial derivative, so this check
    does not depend on the production FEM/WLS discretisation.  The finite
    radial width accounts for the small residual relative to the single-ring
    target.  The check covers both the intrinsic offset and the +1 envelope
    winding shift.
    """

    results = [evaluate("periodic", n, k0, ell, D=D, sigma_k=sigma_k)
               for ell in (0, 1)]
    errors = [result.lambda_exact_gradient - result.bridge_target
              for result in results]
    shift = (results[1].lambda_exact_gradient -
             results[0].lambda_exact_gradient)
    passed = max(abs(error) for error in errors) <= tolerance and abs(shift - 1.0) <= tolerance
    print(
        f"bridge oracle: F_n={results[0].fishman_f:.8f} "
        f"lambda0={results[0].lambda_exact_gradient:.8f} "
        f"target0={results[0].bridge_target:.8f} "
        f"lambda1={results[1].lambda_exact_gradient:.8f} "
        f"target1={results[1].bridge_target:.8f} "
        f"shift={shift:.8f} tolerance={tolerance:.3f} "
        f"status={'PASS' if passed else 'FAIL'}"
    )
    return passed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selftest", action="store_true")
    parser.add_argument("--validate", action="store_true",
                        help="validate the settled particle-field bridge convention")
    parser.add_argument("--reference-only", action="store_true",
                        help="print only the exact gauge-invariant k-space reference")
    parser.add_argument("--packet", choices=("fixed", "ring", "periodic"), default="ring")
    parser.add_argument("--n", type=int, default=24, help="cells per Bravais direction")
    parser.add_argument("--k0", type=float, default=2.9)
    parser.add_argument("--ell", type=int, action="append", default=None,
                        help="winding to evaluate; may be repeated (default: 0,1)")
    parser.add_argument("--d", type=float, default=0.3, dest="D")
    parser.add_argument("--max-neighbors", type=int, default=12)
    parser.add_argument("--sigma-k", type=float, default=0.18)
    parser.add_argument("--convergence", action="store_true",
                        help="run periodic packets at n=12,16,20,24")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return 0
    if args.validate:
        return 0 if validate_particle_bridge(args.n, args.k0, D=args.D,
                                              sigma_k=args.sigma_k) else 1
    ells = args.ell if args.ell is not None else [0, 1]
    if args.convergence:
        for n in (12, 16, 20, 24):
            for ell in ells:
                coords = build_coords(n)
                psi, gx, gy, nmodes = periodic_band_packet(
                    coords, n, args.k0, ell, D=args.D, sigma_k=args.sigma_k
                )
                wx, wy = all_site_wls_gradient(coords, psi, n, args.max_neighbors)
                exact = oam_from_gradient(coords, psi, gx, gy)
                wls = oam_from_gradient(coords, psi, wx, wy)
                reference = exact_kspace_reference(args.k0, ell, D=args.D)
                print(f"convergence n={n} ell={ell} modes={nmodes} "
                     f"lambda_exact={exact:.8f} lambda_wls={wls:.8f} "
                     f"bridge_target={reference.bridge_target:.8f} "
                     f"exact_error={exact - reference.bridge_target:+.8f}")
        return 0

    for ell in ells:
        result = (exact_kspace_reference(args.k0, ell, D=args.D)
                  if args.reference_only else
                  evaluate(args.packet, args.n, args.k0, ell,
                           D=args.D, sigma_k=args.sigma_k,
                           max_neighbors=args.max_neighbors))
        print(
            f"packet={result.packet} ell={result.ell} k0={result.k0:.8f} "
            f"gamma={result.gamma:.8f} F_n={result.fishman_f:.8f} "
            f"lambda_all_site={result.lambda_all_site:.8f} "
            f"lambda_exact_gradient={result.lambda_exact_gradient:.8f} "
            f"bridge_target_l_minus_2F={result.bridge_target:.8f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
