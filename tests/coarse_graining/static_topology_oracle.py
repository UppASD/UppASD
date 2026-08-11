"""RCG-04F independent static fine/buffer/coarse block-ownership oracle.

This module computes the *expected* per-block and per-atom ownership
classification (``FINE``/``BUFFER``["interface"]/``COARSE``) for a static
AdaptiveCG mask, purely from (a) ``moving_state_generator.Geometry``'s own
atom-index/position formulas, (b) an exchange ``jfile``'s declared shell
distances, and (c) the caller's requested block shape and FINE seed block
ids -- never by reading back any production diagnostic. It independently
reproduces the *algorithm* documented and read from
``source/CoarseGraining/statichybridoperator.f90:rebuild_static_hybrid_ownership``
(periodic block dilation of every FINE seed block by a per-axis
``buffer_width_blocks`` vector, ``ceil(radius*|inverse_block_vectors(axis,:)| -
64*eps) + cg_buffer_blocks`` blocks on each axis, Chebyshev/box distance
under periodic wrap -- a genuine box test, ``all(periodic_delta <= width)``,
not a spherical one) and ``source/CoarseGraining/blocktopology.f90`` (block
numbering ``block_id = 1 + bx + grid_x*(by + grid_y*bz)``, atom-to-block via
floor division of each atom's 0-based cell index by the block shape, and
``block_vectors(:,axis) = block_shape(axis)*cell_vectors(:,axis)`` for the
buffer-width metric), so that RCG-04F/RCG-05 can assert the *runtime*
``AdaptiveCG: resolution_state``/``interface_atoms=`` diagnostics agree with
an independently derived expectation, per the RCG-04F governing rule ("do
not infer ownership solely from the input mask... compare the runtime map
or diagnostics with the expected map").

RCG-04F originally restricted this module to the ``block_grid_y ==
block_grid_z == 1`` case every RCG-04 fixture used (``block_size_y=n2``,
``block_size_z=n3``, i.e. the state and the mask both vary only along the
x/modulation axis), raising ``NotImplementedError`` for a genuine 3D block
grid rather than returning a silently wrong answer. RCG-05B removes that
restriction: ``compute_expected_topology`` now computes a genuine 3D
``(grid_x, grid_y, grid_z)`` block partition and a per-axis buffer-width
vector for an arbitrary (including non-orthogonal/skew) ``cell_vectors``,
using the exact general block-vector-inverse metric
``statichybridoperator.f90:180-187``/``inverse3`` uses (not the simplified
``radius/block_size_x`` shortcut, which is only the identity-cell special
case of that same formula). Every existing RCG-04 call site passes
``block_size_y=geometry.n2``, ``block_size_z=geometry.n3`` and no
``cell_vectors`` (implying the default identity cell): this generalization
reduces to the exact old 1D formula and old ``nblocks_x`` semantics in that
case (``nblocks_x`` now holds the *total* block count, numerically
identical to the old per-axis count whenever ``grid_y == grid_z == 1``), so
no existing caller's behaviour changes -- checked here by regression test,
not merely asserted (see ``test_static_topology_oracle.py``).
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

from moving_state_generator import Geometry, content_hash
from torque_oracle import ExchangeShell, MalformedJfileError, build_geometric_bonds

# Matches STATIC_HYBRID_COARSE/BUFFER/FINE (source/CoarseGraining/statichybridoperator.f90:45-47).
COARSE = 0
BUFFER = 1
FINE = 2
_STATE_NAME = {COARSE: "coarse", BUFFER: "interface", FINE: "fine"}

# Matches Fortran ``64.0_dblprec*epsilon(1.0_dblprec)`` (blocktopology.f90/
# statichybridoperator.f90 both use this exact guard band before ceiling).
_EPS_GUARD = 64.0 * sys.float_info.epsilon

_IDENTITY_CELL: tuple[tuple[float, float, float], ...] = (
    (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
)

Vector = tuple[float, float, float]
Matrix3 = tuple[Vector, Vector, Vector]


class UnsupportedTopologyError(NotImplementedError):
    """No longer raised by this module; kept only for import-compatibility."""


class DegenerateGeometryError(ValueError):
    """Raised when ``cell_vectors``/``block_shape`` do not span a finite nonzero volume."""


def max_interaction_radius_m(shells: tuple[ExchangeShell, ...]) -> float:
    return max(shell.distance for shell in shells)


def _validate_cell_vectors(cell_vectors: Matrix3) -> None:
    if len(cell_vectors) != 3 or any(len(vector) != 3 for vector in cell_vectors):
        raise ValueError(f"cell_vectors must be exactly 3 vectors of 3 components each, got {cell_vectors!r}")
    if not all(math.isfinite(component) for vector in cell_vectors for component in vector):
        raise ValueError(f"cell_vectors must contain only finite values, got {cell_vectors!r}")


def _inverse3(matrix: Matrix3) -> tuple[Matrix3, float]:
    """General 3x3 matrix inverse, term-for-term identical to Fortran's ``inverse3``
    (``source/CoarseGraining/statichybridoperator.f90:414-430`` and its sibling
    copies), so this is the *same* formula production evaluates, not merely an
    equivalent one. ``matrix[row][col]`` matches Fortran's ``matrix(row,col)``
    (1-based there, 0-based here); returns ``(inverse, determinant)``.
    """
    m = matrix
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
           - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
           + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    if not math.isfinite(det) or abs(det) <= sys.float_info.min:
        raise DegenerateGeometryError(f"matrix {matrix!r} is singular (determinant {det!r})")
    inv = (
        (
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
        ),
        (
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
        ),
        (
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
        ),
    )
    return inv, det


def _block_vectors_matrix(cell_vectors: Matrix3, block_shape: tuple[int, int, int]) -> Matrix3:
    """Fortran ``block_vectors(:,axis) = block_shape(axis)*cell_vectors(:,axis)``
    (``blocktopology.f90:244-247``), returned as ``matrix[component][axis]`` so
    its *columns* are the three (possibly skew) block vectors -- the same
    row/column convention ``_inverse3`` and production's own ``inverse3`` use.
    """
    scaled = [tuple(block_shape[axis] * component for component in cell_vectors[axis]) for axis in range(3)]
    return tuple(tuple(scaled[axis][component] for axis in range(3)) for component in range(3))


def buffer_width_blocks(
    shells: tuple[ExchangeShell, ...], cell_vectors: Matrix3,
    block_shape: tuple[int, int, int], cg_buffer_blocks: int = 0,
) -> tuple[int, int, int]:
    """General per-axis dilation width, independently reproducing
    ``setup_static_hybrid_operator``'s ``buffer_width_blocks(block) =
    ceiling(max(0,radius*|inverse3(block_vectors)(block,:)| - 64*eps)) +
    safety_dilation_blocks`` (``statichybridoperator.f90:180-187``) for an
    arbitrary, possibly non-orthogonal, ``block_vectors`` metric -- not just
    the identity-cell ``radius/block_size`` shortcut ``buffer_width_blocks_x``
    still provides for the 1-axis case.
    """
    if any(size <= 0 for size in block_shape):
        raise ValueError(f"block_shape must be positive on every axis, got {block_shape}")
    if cg_buffer_blocks < 0:
        raise ValueError(f"cg_buffer_blocks must be nonnegative, got {cg_buffer_blocks}")
    _validate_cell_vectors(cell_vectors)
    radius = max_interaction_radius_m(shells)
    matrix = _block_vectors_matrix(cell_vectors, block_shape)
    inverse, _det = _inverse3(matrix)
    widths = []
    for axis in range(3):
        row = inverse[axis]
        fractional = radius * math.sqrt(sum(component * component for component in row))
        widths.append(math.ceil(max(0.0, fractional - _EPS_GUARD)) + cg_buffer_blocks)
    return (widths[0], widths[1], widths[2])


def buffer_width_blocks_x(
    shells: tuple[ExchangeShell, ...], block_size_x: int, cg_buffer_blocks: int = 0,
) -> int:
    """Independent reproduction of ``setup_static_hybrid_operator``'s per-axis dilation width,
    for the identity-cell, single-axis case (kept for every existing caller).

    Defined as the axis-0 component of the general ``buffer_width_blocks`` on
    an identity cell with ``block_shape=(block_size_x, 1, 1)`` -- so this is
    now a regression check that the general formula reduces to the original
    one, not a second independent implementation.
    """
    if block_size_x <= 0:
        raise ValueError(f"block_size_x must be positive, got {block_size_x}")
    if cg_buffer_blocks < 0:
        raise ValueError(f"cg_buffer_blocks must be nonnegative, got {cg_buffer_blocks}")
    return buffer_width_blocks(shells, _IDENTITY_CELL, (block_size_x, 1, 1), cg_buffer_blocks)[0]


def _block_id(coord: tuple[int, int, int], grid: tuple[int, int, int]) -> int:
    """1-based block id, matching ``regular_spatial_block_id``
    (``blocktopology.f90:92-97``): ``1 + coord(1) + grid(1)*(coord(2) + grid(2)*coord(3))``."""
    return 1 + coord[0] + grid[0] * (coord[1] + grid[1] * coord[2])


def _block_coord(block_id: int, grid: tuple[int, int, int]) -> tuple[int, int, int]:
    """Inverse of ``_block_id``."""
    zero_based = block_id - 1
    c0 = zero_based % grid[0]
    remainder = zero_based // grid[0]
    c1 = remainder % grid[1]
    c2 = remainder // grid[1]
    return (c0, c1, c2)


@dataclass(frozen=True)
class ExpectedTopology:
    geometry: Geometry
    block_size: tuple[int, int, int]
    block_grid: tuple[int, int, int]
    nblocks_x: int  # total block count (grid[0]*grid[1]*grid[2]); name kept for every existing caller
    buffer_width: tuple[int, int, int]
    fine_block_ids: frozenset[int]
    block_state_by_id: dict[int, int]  # 1-based block id -> COARSE/BUFFER/FINE
    atom_block_by_atom: dict[int, int]  # 1-based atom -> 1-based block id
    distance_from_boundary_blocks: dict[int, int]  # 1-based block id -> Chebyshev (L-infinity, grid-index) blocks to nearest differently-classed block

    @property
    def nblocks(self) -> int:
        """Preferred name for the total block count; ``nblocks_x`` is kept only
        because every RCG-04 caller already depends on that exact attribute name."""
        return self.nblocks_x

    @property
    def buffer_width_x(self) -> int:
        """Axis-0 component of ``buffer_width``, kept for every existing caller
        that only ever used the 1D (identity-cell, x-only) case."""
        return self.buffer_width[0]

    def block_ids_in_state(self, state: int) -> list[int]:
        return sorted(bid for bid, s in self.block_state_by_id.items() if s == state)

    def atoms_per_block(self) -> int:
        return self.geometry.na * self.block_size[0] * self.block_size[1] * self.block_size[2]

    def block_counts(self) -> dict[str, int]:
        counts = {"fine": 0, "interface": 0, "coarse": 0}
        for state in self.block_state_by_id.values():
            counts[_STATE_NAME[state]] += 1
        return counts

    def atom_counts(self) -> dict[str, int]:
        per_block = self.atoms_per_block()
        return {name: count * per_block for name, count in self.block_counts().items()}

    def atom_state(self, atom: int) -> int:
        return self.block_state_by_id[self.atom_block_by_atom[atom]]

    def atom_distance_from_boundary_blocks(self, atom: int) -> int:
        return self.distance_from_boundary_blocks[self.atom_block_by_atom[atom]]

    def resolution_state_values(self) -> tuple[int, ...]:
        """Per-block state vector in ``block_id`` order (1..nblocks_x), matching the
        exact ordering ``print_resolution_state`` prints (block loop 1..n_spatial_blocks)."""
        return tuple(self.block_state_by_id[bid] for bid in range(1, self.nblocks_x + 1))


def _dilate_periodic_box(
    fine_block_ids: set[int], grid: tuple[int, int, int], width: tuple[int, int, int],
) -> dict[int, int]:
    """Core of ``rebuild_static_hybrid_ownership``'s periodic Chebyshev/box
    dilation (``statichybridoperator.f90:196-256``): every block within
    ``width`` blocks (per axis, independently, under periodic wrap) of any
    FINE seed becomes BUFFER; the seeds themselves stay FINE; everything else
    is COARSE. Factored out of ``compute_expected_topology`` so the isotropic
    (``width=(m,m,m)`` for ``m=max(width)``) cross-check
    ``compute_isotropic_dilation_topology`` below shares the exact same
    periodic-box logic and cannot silently diverge from it (RCG-05C).
    """
    nblocks = grid[0] * grid[1] * grid[2]
    fine_coords = [_block_coord(bid, grid) for bid in fine_block_ids]
    atomistic: set[int] = set(fine_block_ids)
    for bid in range(1, nblocks + 1):
        if bid in atomistic:
            continue
        coord = _block_coord(bid, grid)
        for seed_coord in fine_coords:
            delta = tuple(abs(coord[axis] - seed_coord[axis]) for axis in range(3))
            periodic_delta = tuple(min(delta[axis], grid[axis] - delta[axis]) for axis in range(3))
            if all(periodic_delta[axis] <= width[axis] for axis in range(3)):
                atomistic.add(bid)
                break

    block_state_by_id: dict[int, int] = {}
    for bid in range(1, nblocks + 1):
        if bid in fine_block_ids:
            block_state_by_id[bid] = FINE
        elif bid in atomistic:
            block_state_by_id[bid] = BUFFER
        else:
            block_state_by_id[bid] = COARSE
    return block_state_by_id


def _atom_block_map(
    geometry: Geometry, block_shape: tuple[int, int, int], grid: tuple[int, int, int],
) -> dict[int, int]:
    atom_block_by_atom: dict[int, int] = {}
    for atom, _i0, i1, i2, i3 in geometry.iter_atoms():
        coord = (i1 // block_shape[0], i2 // block_shape[1], i3 // block_shape[2])
        atom_block_by_atom[atom] = _block_id(coord, grid)
    return atom_block_by_atom


def _distance_from_boundary(
    block_state_by_id: dict[int, int], grid: tuple[int, int, int],
) -> dict[int, int]:
    nblocks = grid[0] * grid[1] * grid[2]
    distance_from_boundary_blocks: dict[int, int] = {}
    block_coords = {bid: _block_coord(bid, grid) for bid in range(1, nblocks + 1)}
    for bid in range(1, nblocks + 1):
        own_state = block_state_by_id[bid]
        coord = block_coords[bid]
        best = max(grid)
        for bid2 in range(1, nblocks + 1):
            if block_state_by_id[bid2] == own_state:
                continue
            coord2 = block_coords[bid2]
            delta = tuple(abs(coord[axis] - coord2[axis]) for axis in range(3))
            periodic_delta = tuple(min(delta[axis], grid[axis] - delta[axis]) for axis in range(3))
            chebyshev = max(periodic_delta)
            best = min(best, chebyshev)
        distance_from_boundary_blocks[bid] = best
    return distance_from_boundary_blocks


def _block_shape_and_grid(
    geometry: Geometry, *, block_size_x: int, block_size_y: int, block_size_z: int,
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    block_shape = (block_size_x, block_size_y, block_size_z)
    if any(size <= 0 for size in block_shape):
        raise ValueError(f"block_size_x/y/z must all be positive, got {block_shape}")
    repetitions = (geometry.n1, geometry.n2, geometry.n3)
    if any(repetitions[axis] % block_shape[axis] != 0 for axis in range(3)):
        raise ValueError(
            f"(n1,n2,n3)={repetitions} must each be divisible by the matching "
            f"block_size, got block_size={block_shape}"
        )
    grid = tuple(repetitions[axis] // block_shape[axis] for axis in range(3))
    return block_shape, grid


def compute_expected_topology(
    geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
    block_size_x: int, block_size_y: int, block_size_z: int,
    fine_block_ids: set[int], cg_buffer_blocks: int = 0,
    cell_vectors: Matrix3 = _IDENTITY_CELL,
) -> ExpectedTopology:
    """Compute the expected fine/buffer/coarse ownership map for a genuine 3D
    block grid and an arbitrary (including skew) ``cell_vectors`` metric.

    ``cell_vectors`` defaults to the identity cell every RCG-04 caller
    implicitly used; every RCG-04 call site also happens to pass
    ``block_size_y=geometry.n2, block_size_z=geometry.n3`` (a degenerate
    ``grid_y == grid_z == 1`` grid), for which this reduces to exactly the
    original 1D formula and ``nblocks_x``/``buffer_width_x`` semantics
    (checked by regression test against real RCG-04 fixture geometries, not
    merely asserted).
    """
    block_shape, grid = _block_shape_and_grid(
        geometry, block_size_x=block_size_x, block_size_y=block_size_y,
        block_size_z=block_size_z,
    )
    nblocks = grid[0] * grid[1] * grid[2]
    if not fine_block_ids or not fine_block_ids <= set(range(1, nblocks + 1)):
        raise ValueError(f"fine_block_ids {fine_block_ids} must be a nonempty subset of 1..{nblocks}")

    width = buffer_width_blocks(shells, cell_vectors, block_shape, cg_buffer_blocks)
    block_state_by_id = _dilate_periodic_box(fine_block_ids, grid, width)
    atom_block_by_atom = _atom_block_map(geometry, block_shape, grid)
    distance_from_boundary_blocks = _distance_from_boundary(block_state_by_id, grid)

    return ExpectedTopology(
        geometry=geometry, block_size=block_shape, block_grid=grid,
        nblocks_x=nblocks, buffer_width=width, fine_block_ids=frozenset(fine_block_ids),
        block_state_by_id=block_state_by_id, atom_block_by_atom=atom_block_by_atom,
        distance_from_boundary_blocks=distance_from_boundary_blocks,
    )


def compute_isotropic_dilation_topology(
    geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
    block_size_x: int, block_size_y: int, block_size_z: int,
    fine_block_ids: set[int], cg_buffer_blocks: int = 0,
    cell_vectors: Matrix3 = _IDENTITY_CELL,
) -> ExpectedTopology:
    """RCG-05C cross-check: the ownership map production's GPU/CUDA/HIP path
    *actually* computes, reproducing the isotropic-cube dilation defect
    (``gpu_buffer_dilation = int(maxval(buffer_width_blocks))``,
    ``adaptivecgproduction.f90:615-616``, consumed by ``dilateAdaptiveState``'s
    single-scalar-radius periodic box scan, ``gpuAdaptiveRuntime.cpp:322-349``)
    rather than the correct per-axis ``buffer_width_blocks(3)``
    ``compute_expected_topology`` above uses.

    Same signature and same ``_dilate_periodic_box`` core as
    ``compute_expected_topology`` -- the *only* difference is collapsing the
    correct per-axis width to ``(m,m,m)`` with ``m=max(width)`` before
    dilating, exactly mirroring the Fortran ``maxval`` collapse. This is a
    test-oracle cross-check, not a second production code path: it lets a
    comparator ask "does this backend's actual reported map match the
    *correct* directional dilation of its own reported FINE seeds, or does it
    match the *isotropic* one instead" -- the precise question RCG-05C's
    defect demonstration needs, tying an observed mismatch back to this exact
    scalarization rather than leaving it as an unexplained raw diff.
    """
    block_shape, grid = _block_shape_and_grid(
        geometry, block_size_x=block_size_x, block_size_y=block_size_y,
        block_size_z=block_size_z,
    )
    nblocks = grid[0] * grid[1] * grid[2]
    if not fine_block_ids or not fine_block_ids <= set(range(1, nblocks + 1)):
        raise ValueError(f"fine_block_ids {fine_block_ids} must be a nonempty subset of 1..{nblocks}")

    correct_width = buffer_width_blocks(shells, cell_vectors, block_shape, cg_buffer_blocks)
    isotropic = max(correct_width)
    width = (isotropic, isotropic, isotropic)
    block_state_by_id = _dilate_periodic_box(fine_block_ids, grid, width)
    atom_block_by_atom = _atom_block_map(geometry, block_shape, grid)
    distance_from_boundary_blocks = _distance_from_boundary(block_state_by_id, grid)

    return ExpectedTopology(
        geometry=geometry, block_size=block_shape, block_grid=grid,
        nblocks_x=nblocks, buffer_width=width, fine_block_ids=frozenset(fine_block_ids),
        block_state_by_id=block_state_by_id, atom_block_by_atom=atom_block_by_atom,
        distance_from_boundary_blocks=distance_from_boundary_blocks,
    )


def write_mask_file(topology: ExpectedTopology) -> str:
    """Render the ``cg_static_mask_file`` text this topology corresponds to.

    Only FINE block ids are listed (one-based, matching production's parser
    at ``adaptivecgproduction.f90:~1305-1328``); every omitted id defaults to
    COARSE, the same convention ``e2e/static_mixed/mask.dat`` and
    ``e2e/moving_all_coarse_bs*/mask.dat`` already use.
    """
    lines = ["# Omitted one-based block ids default to COARSE."]
    for bid in sorted(topology.fine_block_ids):
        lines.append(f"{bid} FINE")
    return "\n".join(lines) + "\n"


def interface_bond_count(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                          topology: ExpectedTopology) -> int:
    """Independent count of active exchange bonds crossing the atomistic/coarse boundary.

    A bond ``i -> j`` is counted when ``i`` is atomistic (FINE or BUFFER) and
    ``j`` is COARSE -- exactly the ``atomistic_bond_owner`` condition read
    from ``statichybridoperator.f90:254-258`` (``active_bond .and.
    (atomistic_atom(i) .or. atomistic_atom(j))``) restricted to the boundary
    case where the bond's atomistic endpoint's neighbour is *not* itself
    atomistic, i.e. the bond genuinely needs the ghost-prolongated coarse
    direction (``evaluate_static_hybrid_operator``'s ``effective_direction``
    for a coarse-owned atom) rather than a real neighbouring atomistic
    direction. A nonzero count is purely topological (geometry/jfile/mask
    derived) proof that the interface/buffer coupling path is structurally
    engaged for this fixture, independent of any single production run.

    Uses ``torque_oracle.build_geometric_bonds``'s plain per-axis minimum
    image, exactly as RCG-04F calibrated it -- correct only for an
    orthogonal (in particular, identity) cell; use
    ``interface_bond_count_general`` for a non-orthogonal/skew cell.
    """
    bonds = build_geometric_bonds(geometry, shells)
    count = 0
    for atom, neighbours in bonds.items():
        if topology.atom_state(atom) == COARSE:
            continue
        for neighbour, _jij in neighbours:
            if topology.atom_state(neighbour) == COARSE:
                count += 1
    return count


def _lattice_minimum_image(displacement: Vector, cell_vectors: Matrix3) -> Vector:
    """General periodic-image reduction under an arbitrary (including skew)
    lattice, matching production's own convention exactly
    (``source/Hamiltonian/neighbourmap.f90:266-290``: ``dfrac = dcart @
    inverse(cell_matrix)``, ``itrans = nint(dfrac)``, ``reduced = dcart -
    itrans @ cell_matrix``) -- not the per-axis orthogonal-box wrap
    ``torque_oracle._minimum_image`` uses, which is only correct when the
    cell is orthogonal (``validate_minimum_image_safe`` does not check
    orthogonality, only per-axis shell-vs-half-box-length safety, so it
    would not by itself catch a skew cell being handled incorrectly).
    """
    matrix = tuple(tuple(cell_vectors[axis][component] for axis in range(3)) for component in range(3))
    inverse, _det = _inverse3(matrix)
    fractional = tuple(
        sum(inverse[row][col] * displacement[col] for col in range(3)) for row in range(3)
    )
    translation = tuple(round(component) for component in fractional)
    reduced = tuple(
        displacement[row] - sum(matrix[row][col] * translation[col] for col in range(3))
        for row in range(3)
    )
    return reduced


def build_geometric_bonds_general(
    geometry: Geometry, shells: tuple[ExchangeShell, ...], cell_vectors: Matrix3, *,
    distance_tolerance: float = 1.0e-6,
) -> dict[int, list[tuple[int, float]]]:
    """Independent geometric neighbour list for an arbitrary (including
    non-orthogonal/skew) ``cell_vectors``, generalizing
    ``torque_oracle.build_geometric_bonds`` (which is restricted to an
    orthogonal cell by its own calibrated per-axis minimum-image convention)
    via ``_lattice_minimum_image``'s general lattice reduction. Same
    single-count-per-physical-neighbour convention (see
    ``torque_oracle.build_geometric_bonds`` docstring for the production
    calibration this convention was validated against).

    Periodicity is over the *simulation box* (``geometry.n1/n2/n3`` replicas
    of ``cell_vectors``), matching ``torque_oracle._minimum_image``'s own
    box-scale convention (``box = (n1*alat, n2*alat, n3*alat)``) exactly --
    not the bare unit cell, which would only find the nearest unit-cell
    replica rather than the nearest *periodic-box* image.
    """
    _validate_cell_vectors(cell_vectors)
    box_vectors = tuple(
        tuple(component * repetitions for component in cell_vectors[axis])
        for axis, repetitions in enumerate((geometry.n1, geometry.n2, geometry.n3))
    )
    positions: dict[int, Vector] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        positions[atom] = geometry.cartesian(i0, i1, i2, i3)

    shell_by_distance = [(shell.distance, shell.jij) for shell in shells]
    bonds: dict[int, list[tuple[int, float]]] = {atom: [] for atom in positions}
    natom = geometry.natom
    for i in range(1, natom + 1):
        pos_i = positions[i]
        for j in range(1, natom + 1):
            if i == j:
                continue
            raw = tuple(positions[j][c] - pos_i[c] for c in range(3))
            reduced = _lattice_minimum_image(raw, box_vectors)
            distance = math.sqrt(sum(component * component for component in reduced))
            matches = [
                jij for shell_distance, jij in shell_by_distance
                if abs(distance - shell_distance) <= distance_tolerance
            ]
            if len(matches) > 1:
                raise MalformedJfileError(
                    f"ambiguous shell match for atoms {i}->{j}: distance {distance} "
                    f"matches {len(matches)} shells within tolerance {distance_tolerance}"
                )
            if matches:
                bonds[i].append((j, matches[0]))
    return bonds


def interface_bond_count_general(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                                  topology: ExpectedTopology, cell_vectors: Matrix3) -> int:
    """``interface_bond_count``, but using ``build_geometric_bonds_general`` so it
    is correct for a non-orthogonal/skew ``cell_vectors`` too."""
    bonds = build_geometric_bonds_general(geometry, shells, cell_vectors)
    count = 0
    for atom, neighbours in bonds.items():
        if topology.atom_state(atom) == COARSE:
            continue
        for neighbour, _jij in neighbours:
            if topology.atom_state(neighbour) == COARSE:
                count += 1
    return count


# ---------------------------------------------------------------------------
# RCG-05B: deterministic block-geometry fixture generators
# ---------------------------------------------------------------------------

GENERATOR_VERSION = "rcg05b-block-geometry-generator-v1"

# The same 2-atom (BCC-like) basis every RCG-04 e2e fixture's shared
# ``posfile`` uses (``tests/coarse_graining/e2e/posfile``); reused rather
# than invented, since ``posfiletype`` defaults to Cartesian
# (``source/Input/inputdata.f90:335``), so this basis offset is valid
# unchanged under any ``cell_vectors``, skew or not
# (``source/Input/inputhandler_ext.f90:1079-1080``: ``posfiletype=='C'`` ->
# ``r_red=r_tmp`` directly, no cell-vector conversion of the basis).
_GENERATOR_BASIS: tuple[tuple[float, float, float], ...] = ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5))


def _manifest(generator_name: str, geometry: Geometry, parameters: dict,
              generated_text: str, extra: dict | None = None) -> dict:
    """Same provenance-record shape as ``moving_state_generator._manifest``
    (generator name/version, geometry, parameters, content hash, extra) --
    reimplemented locally, not imported, only because this module's
    ``GENERATOR_VERSION`` is its own and must not be shadowed by
    ``moving_state_generator``'s; ``content_hash`` itself
    (version-independent) is imported and reused directly.
    """
    record = {
        "generator_name": generator_name,
        "generator_version": GENERATOR_VERSION,
        "geometry": {
            "na": geometry.na, "n1": geometry.n1, "n2": geometry.n2, "n3": geometry.n3,
            "natom": geometry.natom,
            "basis": [list(coordinate) for coordinate in geometry.basis],
            "cell_vectors": [list(vector) for vector in geometry.cell_vectors],
        },
        "parameters": parameters,
        "content_sha256": content_hash(generated_text),
    }
    if extra:
        record["extra"] = extra
    return record


@dataclass(frozen=True)
class BlockGeometryFixture:
    """A deterministic, self-validated block-topology fixture: enough to
    write ``jfile``/``cg_static_mask_file``/the relevant ``inpsd.dat``
    records for a real UppASD run (RCG-05C/D), plus this module's own
    independently computed expected ownership map for the identical
    parameters (RCG-05B does not itself run the real executable or claim
    any CPU/GPU equivalence -- see module and prompt-pack scope notes)."""

    jfile_text: str
    mask_text: str
    inpsd_records: dict[str, str]
    geometry: Geometry
    shells: tuple[ExchangeShell, ...]
    expected_topology: ExpectedTopology
    manifest: dict


def _render_jfile(bonds: tuple[tuple[float, float, float, float], ...]) -> str:
    """Render a single-chemical-type ``jfile`` (``isite jsite dx dy dz Jij``),
    matching every existing fixture jfile under ``tests/coarse_graining/e2e``."""
    if not bonds:
        raise ValueError("at least one bond shell is required")
    lines = [f"1 1 {dx:.17g} {dy:.17g} {dz:.17g} {jij:.17g}" for dx, dy, dz, jij in bonds]
    return "\n".join(lines) + "\n"


def _render_cell_block(cell_vectors: Matrix3) -> str:
    return "\n".join(" ".join(f"{component:.17g}" for component in vector) for vector in cell_vectors)


def _require_nonempty_partition(topology: ExpectedTopology, *, context: str) -> None:
    counts = topology.block_counts()
    empty = [name for name, count in counts.items() if count == 0]
    if empty:
        raise ValueError(
            f"{context}: chosen parameters leave {empty} block state(s) empty "
            f"(counts={counts}); the host must be large enough for fine, buffer, "
            "and coarse regions to all be nonempty, per the RCG-05B prompt"
        )


def unequal_width_orthogonal_fixture(
    *, block_size: tuple[int, int, int], ncell: tuple[int, int, int],
    bond: tuple[float, float, float, float], fine_block_ids: set[int] = frozenset({1}),
    cg_buffer_blocks: int = 0, alat: float = 1.0, simid: str = "rcg05bunequal",
) -> BlockGeometryFixture:
    """Deterministic unequal-block-width (``block_size_x/y/z`` not all equal),
    orthogonal (identity) cell fixture, self-validated to leave fine, buffer,
    and coarse regions all nonempty (RCG-05B prompt requirement), by actually
    calling ``compute_expected_topology`` -- the same oracle this fixture's
    consumer would independently use to check it -- rather than asserting the
    host is "large enough" without checking.
    """
    if len(block_size) != 3 or any(not isinstance(v, int) or v <= 0 for v in block_size):
        raise ValueError(f"block_size must be 3 positive integers, got {block_size!r}")
    if len(set(block_size)) == 1:
        raise ValueError(
            f"block_size={block_size} is not anisotropic (unequal-width fixture "
            "requires at least two distinct axis sizes)"
        )
    if len(ncell) != 3 or any(not isinstance(v, int) or v <= 0 for v in ncell):
        raise ValueError(f"ncell must be 3 positive integers, got {ncell!r}")
    if alat <= 0.0:
        raise ValueError(f"alat must be positive, got {alat}")

    geometry = Geometry(na=2, n1=ncell[0], n2=ncell[1], n3=ncell[2], basis=_GENERATOR_BASIS)
    jfile_text = _render_jfile((bond,))
    shells = (ExchangeShell(dx=bond[0], dy=bond[1], dz=bond[2], jij=bond[3]),)

    topology = compute_expected_topology(
        geometry, shells, block_size_x=block_size[0], block_size_y=block_size[1],
        block_size_z=block_size[2], fine_block_ids=set(fine_block_ids),
        cg_buffer_blocks=cg_buffer_blocks, cell_vectors=_IDENTITY_CELL,
    )
    _require_nonempty_partition(topology, context="unequal_width_orthogonal_fixture")

    inpsd_records = {
        "simid": simid, "alat": f"{alat:.17g}", "ncell": " ".join(str(v) for v in ncell),
        "BC": "P P P", "cell": _render_cell_block(_IDENTITY_CELL),
        "block_size_x": str(block_size[0]), "block_size_y": str(block_size[1]),
        "block_size_z": str(block_size[2]), "cg_buffer_blocks": str(cg_buffer_blocks),
    }
    mask_text = write_mask_file(topology)
    parameters = {
        "block_size": list(block_size), "ncell": list(ncell), "bond": list(bond),
        "fine_block_ids": sorted(fine_block_ids), "cg_buffer_blocks": cg_buffer_blocks,
        "alat": alat, "simid": simid,
    }
    manifest = _manifest("unequal_width_orthogonal_fixture", geometry, parameters, jfile_text,
                          extra={"buffer_width_blocks": list(topology.buffer_width),
                                 "block_counts": topology.block_counts()})
    return BlockGeometryFixture(
        jfile_text=jfile_text, mask_text=mask_text, inpsd_records=inpsd_records,
        geometry=geometry, shells=shells, expected_topology=topology, manifest=manifest,
    )


def skew_cell_fixture(
    *, cell_vectors: Matrix3, block_size: tuple[int, int, int], ncell: tuple[int, int, int],
    bond: tuple[float, float, float, float], fine_block_ids: set[int] = frozenset({1}),
    cg_buffer_blocks: int = 0, alat: float = 1.0, simid: str = "rcg05bskew",
) -> BlockGeometryFixture:
    """Deterministic non-orthogonal (skew) ``cell_vectors`` fixture, self-
    validated (like ``unequal_width_orthogonal_fixture``) to leave fine,
    buffer, and coarse regions all nonempty for the given parameters.

    ``cell_vectors`` must actually be non-orthogonal (checked here): at
    least one pair of distinct lattice vectors must have a nonzero dot
    product, or this would silently degrade to an orthogonal-cell fixture
    that happens to call itself "skew".
    """
    _validate_cell_vectors(cell_vectors)
    pairwise_dots = [
        sum(cell_vectors[a][k] * cell_vectors[b][k] for k in range(3))
        for a, b in ((0, 1), (0, 2), (1, 2))
    ]
    if all(abs(dot) <= 1.0e-12 for dot in pairwise_dots):
        raise ValueError(
            f"cell_vectors={cell_vectors!r} is orthogonal (every pairwise dot product "
            "is ~0); a skew-cell fixture requires a genuinely non-orthogonal cell"
        )
    if len(block_size) != 3 or any(not isinstance(v, int) or v <= 0 for v in block_size):
        raise ValueError(f"block_size must be 3 positive integers, got {block_size!r}")
    if len(ncell) != 3 or any(not isinstance(v, int) or v <= 0 for v in ncell):
        raise ValueError(f"ncell must be 3 positive integers, got {ncell!r}")
    if alat <= 0.0:
        raise ValueError(f"alat must be positive, got {alat}")

    geometry = Geometry(na=2, n1=ncell[0], n2=ncell[1], n3=ncell[2], basis=_GENERATOR_BASIS,
                         cell_vectors=cell_vectors)
    jfile_text = _render_jfile((bond,))
    shells = (ExchangeShell(dx=bond[0], dy=bond[1], dz=bond[2], jij=bond[3]),)

    topology = compute_expected_topology(
        geometry, shells, block_size_x=block_size[0], block_size_y=block_size[1],
        block_size_z=block_size[2], fine_block_ids=set(fine_block_ids),
        cg_buffer_blocks=cg_buffer_blocks, cell_vectors=cell_vectors,
    )
    _require_nonempty_partition(topology, context="skew_cell_fixture")

    inpsd_records = {
        "simid": simid, "alat": f"{alat:.17g}", "ncell": " ".join(str(v) for v in ncell),
        "BC": "P P P", "cell": _render_cell_block(cell_vectors),
        "block_size_x": str(block_size[0]), "block_size_y": str(block_size[1]),
        "block_size_z": str(block_size[2]), "cg_buffer_blocks": str(cg_buffer_blocks),
    }
    mask_text = write_mask_file(topology)
    parameters = {
        "cell_vectors": [list(vector) for vector in cell_vectors], "block_size": list(block_size),
        "ncell": list(ncell), "bond": list(bond), "fine_block_ids": sorted(fine_block_ids),
        "cg_buffer_blocks": cg_buffer_blocks, "alat": alat, "simid": simid,
    }
    manifest = _manifest("skew_cell_fixture", geometry, parameters, jfile_text,
                          extra={"buffer_width_blocks": list(topology.buffer_width),
                                 "block_counts": topology.block_counts()})
    return BlockGeometryFixture(
        jfile_text=jfile_text, mask_text=mask_text, inpsd_records=inpsd_records,
        geometry=geometry, shells=shells, expected_topology=topology, manifest=manifest,
    )
