"""RCG-04F independent static fine/buffer/coarse block-ownership oracle.

This module computes the *expected* per-block and per-atom ownership
classification (``FINE``/``BUFFER``["interface"]/``COARSE``) for a static
AdaptiveCG mask, purely from (a) ``moving_state_generator.Geometry``'s own
atom-index/position formulas, (b) an exchange ``jfile``'s declared shell
distances, and (c) the caller's requested block shape and FINE seed block
ids -- never by reading back any production diagnostic. It independently
reproduces the *algorithm* documented and read from
``source/CoarseGraining/statichybridoperator.f90:rebuild_static_hybrid_ownership``
(periodic block dilation of every FINE seed block by
``buffer_width_blocks = ceil(max_interaction_radius_m/block_length - 64*eps) +
cg_buffer_blocks`` blocks, Chebyshev/box distance under periodic wrap) and
``source/CoarseGraining/blocktopology.f90`` (block numbering
``block_id = 1 + bx + grid_x*(by + grid_y*bz)``, atom-to-block via floor
division of each atom's 0-based cell index by the block shape), so that
RCG-04F can assert the *runtime* ``AdaptiveCG: resolution_state``/
``interface_atoms=`` diagnostics agree with an independently derived
expectation, per the RCG-04F governing rule ("do not infer ownership solely
from the input mask... compare the runtime map or diagnostics with the
expected map").

Restricted (like ``coarse_torque_oracle.py``) to the ``block_grid_y ==
block_grid_z == 1`` case every RCG-04 fixture uses (``block_size_y=n2``,
``block_size_z=n3``, i.e. the state and the mask both vary only along the
x/modulation axis) -- stated explicitly here, not silently assumed; a
non-trivial y/z block grid raises ``NotImplementedError`` rather than
returning a silently wrong answer.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass

from moving_state_generator import Geometry
from torque_oracle import ExchangeShell, build_geometric_bonds

# Matches STATIC_HYBRID_COARSE/BUFFER/FINE (source/CoarseGraining/statichybridoperator.f90:45-47).
COARSE = 0
BUFFER = 1
FINE = 2
_STATE_NAME = {COARSE: "coarse", BUFFER: "interface", FINE: "fine"}

# Matches Fortran ``64.0_dblprec*epsilon(1.0_dblprec)`` (blocktopology.f90/
# statichybridoperator.f90 both use this exact guard band before ceiling).
_EPS_GUARD = 64.0 * sys.float_info.epsilon


class UnsupportedTopologyError(NotImplementedError):
    pass


def max_interaction_radius_m(shells: tuple[ExchangeShell, ...]) -> float:
    return max(shell.distance for shell in shells)


def buffer_width_blocks_x(shells: tuple[ExchangeShell, ...], block_size_x: int,
                           cg_buffer_blocks: int = 0) -> int:
    """Independent reproduction of ``setup_static_hybrid_operator``'s per-axis dilation width."""
    if block_size_x <= 0:
        raise ValueError(f"block_size_x must be positive, got {block_size_x}")
    if cg_buffer_blocks < 0:
        raise ValueError(f"cg_buffer_blocks must be nonnegative, got {cg_buffer_blocks}")
    radius = max_interaction_radius_m(shells)
    fractional = radius / block_size_x
    return math.ceil(max(0.0, fractional - _EPS_GUARD)) + cg_buffer_blocks


@dataclass(frozen=True)
class ExpectedTopology:
    geometry: Geometry
    block_size: tuple[int, int, int]
    nblocks_x: int
    buffer_width_x: int
    fine_block_ids: frozenset[int]
    block_state_by_id: dict[int, int]  # 1-based block id -> COARSE/BUFFER/FINE
    atom_block_by_atom: dict[int, int]  # 1-based atom -> 1-based block id
    distance_from_boundary_blocks: dict[int, int]  # 1-based block id -> Chebyshev blocks to nearest differently-classed block

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


def compute_expected_topology(
    geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
    block_size_x: int, block_size_y: int, block_size_z: int,
    fine_block_ids: set[int], cg_buffer_blocks: int = 0,
) -> ExpectedTopology:
    if block_size_y != geometry.n2 or block_size_z != geometry.n3:
        raise UnsupportedTopologyError(
            "compute_expected_topology only supports a block grid that fully spans y/z "
            f"(block_size_y=n2={geometry.n2}, block_size_z=n3={geometry.n3}); got "
            f"block_size_y={block_size_y} block_size_z={block_size_z}"
        )
    if geometry.n1 % block_size_x != 0:
        raise ValueError(f"n1={geometry.n1} is not divisible by block_size_x={block_size_x}")
    nblocks_x = geometry.n1 // block_size_x
    if not fine_block_ids or not fine_block_ids <= set(range(1, nblocks_x + 1)):
        raise ValueError(f"fine_block_ids {fine_block_ids} must be a nonempty subset of 1..{nblocks_x}")

    width = buffer_width_blocks_x(shells, block_size_x, cg_buffer_blocks)

    fine_mask = [False] * nblocks_x
    for bid in fine_block_ids:
        fine_mask[bid - 1] = True

    atomistic_block = list(fine_mask)
    for seed_bx in range(nblocks_x):
        if not fine_mask[seed_bx]:
            continue
        for bx in range(nblocks_x):
            delta = abs(bx - seed_bx)
            periodic_delta = min(delta, nblocks_x - delta)
            if periodic_delta <= width:
                atomistic_block[bx] = True

    block_state_by_id: dict[int, int] = {}
    for bx in range(nblocks_x):
        bid = bx + 1
        if fine_mask[bx]:
            block_state_by_id[bid] = FINE
        elif atomistic_block[bx]:
            block_state_by_id[bid] = BUFFER
        else:
            block_state_by_id[bid] = COARSE

    atom_block_by_atom: dict[int, int] = {}
    for atom, _i0, i1, _i2, _i3 in geometry.iter_atoms():
        atom_block_by_atom[atom] = (i1 // block_size_x) + 1

    distance_from_boundary_blocks: dict[int, int] = {}
    for bx in range(nblocks_x):
        bid = bx + 1
        own_state = block_state_by_id[bid]
        best = nblocks_x
        for bx2 in range(nblocks_x):
            bid2 = bx2 + 1
            if block_state_by_id[bid2] == own_state:
                continue
            delta = abs(bx - bx2)
            periodic_delta = min(delta, nblocks_x - delta)
            best = min(best, periodic_delta)
        distance_from_boundary_blocks[bid] = best

    return ExpectedTopology(
        geometry=geometry, block_size=(block_size_x, block_size_y, block_size_z),
        nblocks_x=nblocks_x, buffer_width_x=width, fine_block_ids=frozenset(fine_block_ids),
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
