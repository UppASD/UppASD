"""RCG-04E independent nontriviality check for the all-coarse block state.

Scope and explicit limitation (read before using this module for anything
beyond a nonzero-torque nontriviality gate)
-----------------------------------------------------------------------------
RCG-04D's ``torque_oracle.py`` independently computes the *atomistic*
initial effective field/torque and calibrates its Hamiltonian convention
against real production ``localenergy.<simid>.out`` output on a uniform
state (two structurally different configurations, matched to 9 significant
figures -- see that module's docstring). No equivalent live calibration is
available for the *coarse-block* exchange-stiffness tensor: production's own
determination (``source/SpinWaves/stiffness.f90:fit_coarse_material``) fits
a least-squares-regularized tensor over several regularization scales
(``eta``) from the atomistic ``jfile`` shells, and every existing coarse
fixture prior to RCG-04E uses a *uniform* block texture, for which the
fitted exchange-stiffness tensor's value is exactly unobservable (a uniform
texture has zero spatial gradient regardless of the tensor's magnitude, so
no existing production output can calibrate it against).

Rather than reverse-engineer that regularized fit from source alone (the
RCG-04D governing rule's explicit warning against "presenting a misleading
formula" when "the analytic mapping would omit enabled terms or lattice
details"), this module instead:

1. Computes an **unregularized, standard second-moment** coarse
   exchange-stiffness estimate directly from the atomistic ``jfile`` shells
   (``D_xx = sum_bonds J(delta) * dx**2 / (2 * V_cell)``, the textbook
   long-wavelength-limit coefficient of ``q**2`` in the atomistic magnon
   dispersion ``E(q) - E(0) = sum_delta J(delta) * (1 - cos(q.delta))``,
   Taylor-expanded to leading order). This coincides with production's fit
   in the zero-regularization limit for this repository's short-ranged,
   Bravais-equivalent ``jfile`` shells, but a quantitative match to
   production's actual regularized value is **not claimed or tested**.
2. Uses that estimate only to build a **discrete-Laplacian torque proxy**
   for the periodic block chain (central second difference across
   neighbouring blocks, the same finite-difference structure
   ``coarsetensoroperator.f90:physical_forward_gradient``/
   ``add_physical_gradient_transpose`` implement, but not a byte-for-byte
   reproduction of that code).
3. Reports **structural, qualitative nonzero-torque/non-collinearity
   evidence** -- that block-average directions genuinely differ between
   neighbouring blocks (a purely geometric fact, independent of any tensor
   value) and that the resulting torque proxy is nonzero. This is
   sufficient to demonstrate the all-coarse state is not a stationary
   special case (RCG-04E checklist: "the all-coarse state has independently
   demonstrated nonzero initial torque"), but it is **not** RCG-04E's
   accuracy oracle. Per the governing rule, that role belongs to the
   accepted all-fine atomistic trajectory (``moving_all_fine_wide``),
   compared directly against each all-coarse fixture's own reconstructed
   trajectory in ``run_moving_all_coarse.py``.

Every input here -- geometry, ``jfile`` shells, and the generator's own
``direction_by_atom`` -- is independent of any UppASD production diagnostic,
matching the RCG-04D/E governing rule that the nontriviality oracle must not
be derived from the production output it gates.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from moving_state_generator import Geometry
from torque_oracle import ExchangeShell, build_geometric_bonds

Vector = tuple[float, float, float]


def _norm(vector: Vector) -> float:
    return math.sqrt(sum(component * component for component in vector))


def _cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def _minimum_image_component(delta: float, length: float) -> float:
    return delta - round(delta / length) * length


def block_average_directions(
    geometry: Geometry, block_size_x: int, direction_by_atom: dict[int, Vector],
) -> tuple[dict[int, Vector], dict[int, Vector], dict[int, float]]:
    """Average per-atom directions into blocks grouped purely by ``i1 // block_size_x``.

    Every RCG-04E all-coarse fixture sets ``block_size_y=n2``,
    ``block_size_z=n3`` (blocks fully span y/z, matching the state's
    uniformity in those directions), so grouping by ``i1`` alone identifies
    the same atom set production's ``BlockTopology`` assigns to each block
    -- but this function computes that grouping itself from ``geometry``,
    never from a production block-id readback, so it remains independent.

    Returns ``(raw_mean, normalized_mean, raw_norm)`` keyed by 0-based block
    index along x. ``raw_norm`` (the pre-normalization vector length) is the
    intra-block averaging reduction factor: it is < 1 whenever atoms in the
    same block point in different directions, and would be exactly 0 only at
    an exact aliasing null (see module docstring point 3).
    """
    if geometry.n1 % block_size_x != 0:
        raise ValueError(
            f"block_size_x={block_size_x} does not divide n1={geometry.n1}"
        )
    nblocks = geometry.n1 // block_size_x
    sums: dict[int, list[float]] = {b: [0.0, 0.0, 0.0] for b in range(nblocks)}
    counts: dict[int, int] = {b: 0 for b in range(nblocks)}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        block = i1 // block_size_x
        direction = direction_by_atom[atom]
        for component in range(3):
            sums[block][component] += direction[component]
        counts[block] += 1
    raw_mean: dict[int, Vector] = {
        b: (sums[b][0] / counts[b], sums[b][1] / counts[b], sums[b][2] / counts[b])
        for b in sums
    }
    raw_norm: dict[int, float] = {b: _norm(raw_mean[b]) for b in raw_mean}
    normalized_mean: dict[int, Vector] = {
        b: (
            (raw_mean[b][0] / raw_norm[b], raw_mean[b][1] / raw_norm[b], raw_mean[b][2] / raw_norm[b])
            if raw_norm[b] > 0.0 else (0.0, 0.0, 0.0)
        )
        for b in raw_mean
    }
    return raw_mean, normalized_mean, raw_norm


def estimate_unregularized_stiffness_xx(
    geometry: Geometry, shells: tuple[ExchangeShell, ...], *, alat: float = 1.0,
    cell_volume: float = 1.0, representative_atom: int = 1,
) -> float:
    """Standard unregularized second-moment ``D_xx`` estimate (see module docstring point 1).

    Computed from the full geometrically-expanded neighbour list of one
    representative atom (reusing ``torque_oracle.build_geometric_bonds``,
    itself independent of production's own neighbour-list code), not merely
    the raw ``jfile`` shell lines -- the ``jfile`` lists one representative
    bond per crystal-symmetry-equivalent shell (``Sym 1``), and a correct
    second-moment sum requires every symmetry-expanded bond, each with its
    own signed ``dx``.
    """
    bonds = build_geometric_bonds(geometry, shells, alat=alat)
    positions: dict[int, Vector] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        positions[atom] = geometry.cartesian(i0, i1, i2, i3)
    box_x = geometry.n1 * alat
    total = 0.0
    for neighbour, jij in bonds[representative_atom]:
        dx = _minimum_image_component(
            positions[neighbour][0] - positions[representative_atom][0], box_x,
        )
        total += jij * dx * dx
    return total / (2.0 * cell_volume)


@dataclass(frozen=True)
class CoarseChainTorqueReport:
    """Discrete-Laplacian torque proxy for a periodic 1D block chain."""

    max_torque_proxy: float
    rms_torque_proxy: float
    max_neighbor_angle_deg: float
    min_block_average_norm: float

    def as_dict(self) -> dict:
        return {
            "max_torque_proxy": self.max_torque_proxy,
            "rms_torque_proxy": self.rms_torque_proxy,
            "max_neighbor_angle_deg": self.max_neighbor_angle_deg,
            "min_block_average_norm": self.min_block_average_norm,
        }


def coarse_chain_nontriviality(
    normalized_block_directions: dict[int, Vector], raw_block_norms: dict[int, float],
    *, block_length_x: float, stiffness_xx: float,
) -> CoarseChainTorqueReport:
    """Central-difference discrete-Laplacian torque proxy (module docstring point 2/3)."""
    nblocks = len(normalized_block_directions)
    torque_magnitudes: list[float] = []
    max_neighbor_angle_deg = 0.0
    for block in range(nblocks):
        d0 = normalized_block_directions[block]
        dp = normalized_block_directions[(block + 1) % nblocks]
        dm = normalized_block_directions[(block - 1) % nblocks]
        laplacian = tuple(
            (dp[c] - 2.0 * d0[c] + dm[c]) / (block_length_x * block_length_x)
            for c in range(3)
        )
        field = tuple(stiffness_xx * component for component in laplacian)
        torque = _cross(d0, field)
        torque_magnitudes.append(_norm(torque))
        cos_angle = max(-1.0, min(1.0, _dot(d0, dp)))
        max_neighbor_angle_deg = max(max_neighbor_angle_deg, math.degrees(math.acos(cos_angle)))
    max_torque = max(torque_magnitudes)
    rms_torque = math.sqrt(sum(t * t for t in torque_magnitudes) / len(torque_magnitudes))
    return CoarseChainTorqueReport(
        max_torque_proxy=max_torque,
        rms_torque_proxy=rms_torque,
        max_neighbor_angle_deg=max_neighbor_angle_deg,
        min_block_average_norm=min(raw_block_norms.values()),
    )
