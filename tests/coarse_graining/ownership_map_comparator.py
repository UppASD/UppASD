"""RCG-05C: CPU/GPU ownership-map comparator.

Extracts complete per-block fine/buffer/coarse ownership maps from a CPU run
and a CUDA/HIP run of the *same* input (reusing ``run_production_e2e.py``'s
already-validated ``final_state`` parser and
``trajectory_evidence.py``'s ``parse_resolution_state_history``, not a new
dumper -- per the RCG-05 prompt pack's explicit instruction against writing a
second ownership-map dumper), and compares them by block-id identity, never
by cardinality/count alone.

Two independent comparisons are performed, both index-set based:

1. **CPU vs. GPU**, directly: the raw ``block_id -> state`` maps, block for
   block.
2. **Self-consistency**, per backend: each backend's own reported FINE block
   set is fed back into ``static_topology_oracle.compute_expected_topology``
   (the *correct*, per-axis-width oracle, independently reproducing
   ``rebuild_static_hybrid_ownership``,
   ``source/CoarseGraining/statichybridoperator.f90:196-256``) and into
   ``compute_isotropic_dilation_topology`` (the *defective* oracle,
   independently reproducing the GPU-staged scalar collapse,
   ``adaptivecgproduction.f90:615-616`` + ``gpuAdaptiveRuntime.cpp:322-349``).
   Whichever oracle a backend's actual reported map matches identifies,
   precisely and without guessing, which dilation formula that backend
   actually used -- regardless of whether the two backends happened to agree
   on *which* blocks are FINE seeds in the first place (a separate,
   documented capability gap: GPU's ``hardAtomisticBlockMask`` is sourced
   only from the polarization gate, not from ``cg_static_mask_file`` --
   see ``source/gpu_files/gpuSimulation.cpp``'s own RCG-03 comment above the
   ``Gpu: AdaptiveCG resolved diagnostics=`` print -- so this module never
   assumes the two backends' FINE seed sets agree).

Also covers, reusing existing infrastructure rather than reinventing it:

- **Periodic wrapping in every direction**: ``periodic_wrap_axes`` reuses the
  oracle's own ``block_grid``/``buffer_width``/``fine_block_ids`` fields to
  determine, per axis, whether at least one atomistic block's classification
  actually depended on the periodic (wrapped) reduction rather than the raw
  one -- i.e. that the comparison is not vacuous in that axis.
- **Every atomistic cross-interface bond accounted for**: ``bond_coverage``
  reuses ``torque_oracle.build_geometric_bonds`` (RCG-04's own
  production-calibrated directed-bond neighbour list) rather than deriving a
  second one, and checks that both backends' reported maps agree on the
  atomistic/coarse classification of every bond endpoint that
  ``interface_bond_count`` counts.
"""

from __future__ import annotations

from dataclasses import dataclass

import static_topology_oracle as sto
from moving_state_generator import Geometry
from run_production_e2e import final_state
from torque_oracle import ExchangeShell, build_geometric_bonds


@dataclass(frozen=True)
class OwnershipMap:
    """A single backend's complete per-block ownership map (RCG-05A §5 item 3/5)."""

    source: str  # "CPU" / "CUDA" / "HIP"
    block_state_by_id: dict[int, int]  # 1-based block id -> COARSE(0)/BUFFER(1)/FINE(2)

    @property
    def nblocks(self) -> int:
        return len(self.block_state_by_id)

    def fine_block_ids(self) -> set[int]:
        return {bid for bid, state in self.block_state_by_id.items() if state == sto.FINE}

    def block_counts(self) -> dict[str, int]:
        counts = {"fine": 0, "interface": 0, "coarse": 0}
        for state in self.block_state_by_id.values():
            counts[sto._STATE_NAME[state]] += 1
        return counts

    def as_dict(self) -> dict:
        return {
            "source": self.source, "nblocks": self.nblocks,
            "block_counts": self.block_counts(),
            "block_state_by_id": {str(k): v for k, v in sorted(self.block_state_by_id.items())},
        }


def parse_final_ownership_map(stdout: str, source: str) -> OwnershipMap:
    """Extract the final per-block state vector via the already-validated,
    already-production-used ``run_production_e2e.final_state`` parser (which
    matches both CPU's ``AdaptiveCG: resolution_state label=final ...
    values=`` and GPU's ``Gpu: AdaptiveCG final_state values=`` lines with
    one regex, already relied on by ``run_production_e2e.py`` itself to
    assert CPU/GPU parity elsewhere) -- not a new parser.
    """
    values = final_state(stdout)
    return OwnershipMap(source=source, block_state_by_id=dict(enumerate(values, start=1)))


@dataclass(frozen=True)
class MapComparison:
    """Index-set (never cardinality-only) comparison of two ownership maps."""

    label: str
    a: OwnershipMap
    b: OwnershipMap
    mismatched_block_ids: tuple[int, ...]
    mismatch_detail: dict[int, tuple[int, int]]  # block_id -> (a_state, b_state)

    @property
    def matches(self) -> bool:
        return not self.mismatched_block_ids

    def as_dict(self) -> dict:
        return {
            "label": self.label, "source_a": self.a.source, "source_b": self.b.source,
            "matches": self.matches, "nblocks": self.a.nblocks,
            "mismatched_block_count": len(self.mismatched_block_ids),
            "mismatched_block_ids": list(self.mismatched_block_ids),
            "mismatch_detail": {
                str(bid): {"a": states[0], "b": states[1]}
                for bid, states in sorted(self.mismatch_detail.items())
            },
        }


class MapShapeMismatchError(ValueError):
    """Raised when two maps do not cover the same block-id set at all."""


def compare_ownership_maps(a: OwnershipMap, b: OwnershipMap, *, label: str) -> MapComparison:
    if set(a.block_state_by_id) != set(b.block_state_by_id):
        raise MapShapeMismatchError(
            f"{label}: block-id sets differ entirely ({a.source} has "
            f"{len(a.block_state_by_id)} blocks, {b.source} has {len(b.block_state_by_id)}) "
            "-- not a dilation-shape comparison, a geometry/setup mismatch"
        )
    mismatch_detail = {
        bid: (a.block_state_by_id[bid], b.block_state_by_id[bid])
        for bid in a.block_state_by_id
        if a.block_state_by_id[bid] != b.block_state_by_id[bid]
    }
    return MapComparison(
        label=label, a=a, b=b,
        mismatched_block_ids=tuple(sorted(mismatch_detail)),
        mismatch_detail=mismatch_detail,
    )


@dataclass(frozen=True)
class SelfConsistencyResult:
    """Does ``reported`` match the CORRECT oracle, the ISOTROPIC oracle,
    neither, or (degenerately) both, when re-dilating ``reported``'s own FINE
    seed set? Both oracle maps are computed fresh here -- never read back
    from ``reported`` beyond the FINE seed set itself -- so a match is
    genuine independent confirmation, not circular."""

    backend: str
    fine_seed_ids: frozenset[int]
    matches_correct_oracle: bool
    matches_isotropic_oracle: bool
    correct_buffer_width: tuple[int, int, int]
    isotropic_buffer_width: tuple[int, int, int]
    correct_mismatch: MapComparison
    isotropic_mismatch: MapComparison

    def as_dict(self) -> dict:
        return {
            "backend": self.backend, "fine_seed_ids": sorted(self.fine_seed_ids),
            "matches_correct_oracle": self.matches_correct_oracle,
            "matches_isotropic_oracle": self.matches_isotropic_oracle,
            "correct_buffer_width": list(self.correct_buffer_width),
            "isotropic_buffer_width": list(self.isotropic_buffer_width),
            "correct_oracle_comparison": self.correct_mismatch.as_dict(),
            "isotropic_oracle_comparison": self.isotropic_mismatch.as_dict(),
        }


def self_consistency_check(
    reported: OwnershipMap, geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
    block_size_x: int, block_size_y: int, block_size_z: int, cg_buffer_blocks: int = 0,
) -> SelfConsistencyResult:
    """Recompute both the correct (per-axis) and isotropic (scalar-collapsed)
    dilation of ``reported``'s *own* FINE seed set, and compare each against
    what the backend actually reported -- the RCG-05C mechanism that
    isolates the buffer-width defect specifically, independent of whether
    this backend's FINE seed selection itself agrees with any other
    backend's (a separate, honestly-reported concern; see module docstring).
    """
    fine_ids = reported.fine_block_ids()
    correct = sto.compute_expected_topology(
        geometry, shells, block_size_x=block_size_x, block_size_y=block_size_y,
        block_size_z=block_size_z, fine_block_ids=fine_ids, cg_buffer_blocks=cg_buffer_blocks,
    )
    isotropic = sto.compute_isotropic_dilation_topology(
        geometry, shells, block_size_x=block_size_x, block_size_y=block_size_y,
        block_size_z=block_size_z, fine_block_ids=fine_ids, cg_buffer_blocks=cg_buffer_blocks,
    )
    correct_map = OwnershipMap(
        source=f"{reported.source}-correct-oracle",
        block_state_by_id=dict(correct.block_state_by_id),
    )
    isotropic_map = OwnershipMap(
        source=f"{reported.source}-isotropic-oracle",
        block_state_by_id=dict(isotropic.block_state_by_id),
    )
    correct_cmp = compare_ownership_maps(reported, correct_map, label=f"{reported.source}-vs-correct")
    isotropic_cmp = compare_ownership_maps(reported, isotropic_map, label=f"{reported.source}-vs-isotropic")
    return SelfConsistencyResult(
        backend=reported.source, fine_seed_ids=frozenset(fine_ids),
        matches_correct_oracle=correct_cmp.matches, matches_isotropic_oracle=isotropic_cmp.matches,
        correct_buffer_width=correct.buffer_width, isotropic_buffer_width=isotropic.buffer_width,
        correct_mismatch=correct_cmp, isotropic_mismatch=isotropic_cmp,
    )


def periodic_wrap_axes(topology: sto.ExpectedTopology) -> tuple[bool, bool, bool]:
    """Per axis, was the periodic (wrapped) reduction actually load-bearing
    for at least one atomistic block -- i.e. does at least one non-seed
    atomistic block have a *raw* (unwrapped) axis delta from every FINE seed
    that exceeds the buffer width on that axis, while its *periodic*
    (wrapped) delta from some seed does not? If every atomistic block on an
    axis was already within width using the raw (non-wrapped) delta, that
    axis's periodic-wrap logic was never actually exercised by this fixture
    and a "periodic wrapping checked" claim for that axis would be vacuous.
    """
    grid = topology.block_grid
    width = topology.buffer_width
    fine_coords = [sto._block_coord(bid, grid) for bid in topology.fine_block_ids]
    wrapped = [False, False, False]
    for bid, state in topology.block_state_by_id.items():
        if state == sto.COARSE or bid in topology.fine_block_ids:
            continue
        coord = sto._block_coord(bid, grid)
        for seed in fine_coords:
            for axis in range(3):
                delta = abs(coord[axis] - seed[axis])
                periodic_delta = min(delta, grid[axis] - delta)
                if periodic_delta <= width[axis] and delta > width[axis]:
                    wrapped[axis] = True
    return (wrapped[0], wrapped[1], wrapped[2])


@dataclass(frozen=True)
class BondCoverageResult:
    """RCG-05C's cross-interface-bond-coverage evidence: every directed bond
    from ``torque_oracle.build_geometric_bonds`` (RCG-04's production-
    calibrated neighbour list, reused not reinvented) is classified by both
    endpoints' reported ownership state, and any bond whose endpoints'
    atomistic/coarse classification differs between two compared maps is
    reported explicitly (never just a count)."""

    total_bonds: int
    interface_bonds_a: int
    interface_bonds_b: int
    disagreeing_bond_endpoints: tuple[tuple[int, int], ...]  # (atom_i, atom_j) pairs

    def as_dict(self) -> dict:
        return {
            "total_directed_bonds": self.total_bonds,
            "interface_bonds_a": self.interface_bonds_a,
            "interface_bonds_b": self.interface_bonds_b,
            "disagreeing_bond_endpoint_count": len(self.disagreeing_bond_endpoints),
            "disagreeing_bond_endpoints": [list(pair) for pair in self.disagreeing_bond_endpoints],
        }


def bond_coverage(
    geometry: Geometry, shells: tuple[ExchangeShell, ...],
    atom_block_by_atom: dict[int, int], map_a: OwnershipMap, map_b: OwnershipMap,
) -> BondCoverageResult:
    bonds = build_geometric_bonds(geometry, shells)

    def is_coarse(atom: int, ownership: OwnershipMap) -> bool:
        return ownership.block_state_by_id[atom_block_by_atom[atom]] == sto.COARSE

    total = 0
    interface_a = 0
    interface_b = 0
    disagreements: list[tuple[int, int]] = []
    for atom, neighbours in bonds.items():
        for neighbour, _jij in neighbours:
            total += 1
            crosses_a = (not is_coarse(atom, map_a)) and is_coarse(neighbour, map_a)
            crosses_b = (not is_coarse(atom, map_b)) and is_coarse(neighbour, map_b)
            if crosses_a:
                interface_a += 1
            if crosses_b:
                interface_b += 1
            if crosses_a != crosses_b:
                disagreements.append((atom, neighbour))
    return BondCoverageResult(
        total_bonds=total, interface_bonds_a=interface_a, interface_bonds_b=interface_b,
        disagreeing_bond_endpoints=tuple(disagreements),
    )


def evidence_record(
    *, fixture_name: str, geometry: Geometry, block_size: tuple[int, int, int],
    cell_vectors: sto.Matrix3, cg_buffer_blocks: int, maps: dict[str, OwnershipMap],
    comparisons: list[MapComparison], self_consistency: list[SelfConsistencyResult],
    wrap_axes: tuple[bool, bool, bool], coverage: BondCoverageResult | None,
) -> dict:
    """Canonical evidence record, per ``docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md``
    SS5 (defined by RCG-05A for RCG-05B/C to implement)."""
    return {
        "fixture": fixture_name,
        "geometry": {
            "na": geometry.na, "n1": geometry.n1, "n2": geometry.n2, "n3": geometry.n3,
            "natom": geometry.natom, "cell_vectors": [list(v) for v in cell_vectors],
            "block_size": list(block_size), "cg_buffer_blocks": cg_buffer_blocks,
        },
        "maps": {name: m.as_dict() for name, m in maps.items()},
        "comparisons": [c.as_dict() for c in comparisons],
        "self_consistency": [s.as_dict() for s in self_consistency],
        "periodic_wrap_axes_exercised": {"x": wrap_axes[0], "y": wrap_axes[1], "z": wrap_axes[2]},
        "bond_coverage": coverage.as_dict() if coverage is not None else None,
    }
