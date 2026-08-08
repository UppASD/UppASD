"""RCG-04B deterministic moving-state generators.

This module produces deterministic, byte-stable initial-state inputs for
later RCG-04 production e2e fixtures. It does not run the UppASD executable,
does not implement a trajectory parser or physical oracle, and does not make
any moving-dynamics or parity claim: see
``docs/RCG-04_MOVING_E2E_EVIDENCE.md`` (RCG-04B section) for what this slice
does and does not establish.

Atom ordering and Cartesian-position conventions below are taken directly
from the production Fortran source, not reinvented:

- the global atom index for basis site ``i0`` (1-based) in unit cell
  ``(i1, i2, i3)`` (0-based) is ``i0 + i1*na + i2*n1*na + i3*n2*n1*na``,
  exactly ``source/System/magnetizationinit.f90``'s ``Initmag`` loops
  (e.g. the ``initmag==8`` spin-spiral branch, lines 449-455) and
  ``source/System/geometry.f90``'s block-ownership loops;
- the Cartesian coordinate of that atom is
  ``coord = i1*C1 + i2*C2 + i3*C3 + Bas(i0)``
  (``source/System/geometry.f90:445``). This module only implements the
  identity-cell case (``cell 1 0 0 / 0 1 0 / 0 0 1``), which is every
  geometry file under ``tests/coarse_graining/e2e`` today; for that case
  ``posfile``'s fractional coordinates and Cartesian coordinates coincide,
  so ``basis`` below is taken directly from ``posfile``-style numbers
  without a separate direct/Cartesian conversion step.

Two output mechanisms are used, matching what AdaptiveCG production
currently accepts (``source/CoarseGraining/adaptivecgproduction.f90:695``
allows ``Initmag`` in ``{1, 2, 3, 5, 8}``; ``Initmag=4`` restart loading is
explicitly rejected at line 692-693, "unsupported until adaptive state
serialization is implemented"):

- ``conical_spiral_state``/``chiral_partner_pair`` produce an ordinary
  ``momfile`` (per-basis template, matching every existing fixture) plus the
  ``Initmag 8``/``initpropvec``/``initrotvec``/``initrotang`` records that
  drive UppASD's *existing* spin-spiral-modulation initializer
  (``magnetizationinit.f90:402-524``). This reuses the same mechanism the
  existing (flagged-degenerate, see RCG-04A evidence section 3.5)
  ``initmag_spin_spiral`` fixture already uses, generalized to a genuine
  cone rather than a planar spiral, and is directly consumable by
  AdaptiveCG production today.
- ``domain_wall_pair_state`` produces a full per-atom restart-format text
  block (matching ``source/System/restart.f90:read_mag_conf_std``'s exact
  expected format), because a domain wall cannot be expressed as a global
  propagation-vector rotation. This is consumable by *ordinary* UppASD
  (``Initmag 4``) today, but **not** by AdaptiveCG production, which
  rejects ``Initmag=4`` outright (cited above). This is a real, pre-existing
  production capability gap discovered while building this generator, not a
  defect introduced by it; RCG-04G (adaptive wall motion) depends on this
  generator and cannot consume it end-to-end until that gap is resolved by
  a separately reviewed change. See the RCG-04B evidence section for the
  full writeup.

Analytic justification for the conical-spiral degenerate-parameter
rejection (independently derived from the written Hamiltonian, not from any
production diagnostic): for a Bravais lattice with isotropic, centrosymmetric
pairwise exchange (every bond has an equal-coupling mirror partner, which
crystal-symmetry-generated neighbour lists such as every ``jfile`` under
``tests/coarse_graining/e2e`` satisfy), a conical-spiral ansatz
``m(x) = R_axis(q.x) . s0`` is an exact fixed point of the LLG effective
field if and only if ``s0`` is an eigenvector of
``M(q) = sum_delta J(delta) R_axis(q.delta)``. Decomposing ``R_axis`` in the
eigenbasis of rotations about ``axis`` (the axis direction itself, plus the
two perpendicular raising/lowering directions) makes ``M(q)`` diagonal with
eigenvalue ``J0 = sum_delta J(delta)`` along the axis and eigenvalue
``J(q) = sum_delta J(delta) cos(q.delta)`` (real, by centrosymmetry) in the
perpendicular plane. A cone angle ``theta`` strictly between 0 and 180
degrees, and not equal to 90 degrees, mixes the axis eigenvector
(``theta=0/180``) and the perpendicular-plane eigenvectors (``theta=90``)
with different eigenvalues whenever ``J0 != J(q)`` -- which is generically
true for any nonzero ``q`` in the long-wavelength regime, where
``J0 - J(q) ~= q^2 * sum_delta J(delta) delta^2 / 2`` (nonzero for any
nonzero exchange stiffness). ``theta`` near 0/180 (collinear/uniform) or 90
degrees (planar spiral -- the RCG-04A-flagged case) is therefore rejected;
``q = 0`` is rejected outright (trivially stationary for any ``theta``).
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math

GENERATOR_VERSION = "rcg04b-moving-state-generator-v1"

# Degenerate cone angles (collinear/uniform at 0 and 180, planar spiral at
# 90) are rejected within this angular tolerance, per the analytic
# derivation in the module docstring.
DEFAULT_DEGENERACY_TOLERANCE_DEG = 5.0


class DegenerateStateError(ValueError):
    """Raised when requested parameters would produce a stationary state."""


class MalformedGeneratorInputError(ValueError):
    """Raised for structurally invalid generator arguments."""


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Geometry:
    """A periodic supercell geometry, matching UppASD's ``posfile``/``ncell``.

    ``basis`` gives each basis atom's fractional-of-cell coordinate in the
    same order as ``posfile`` (1-based site index implied by list order).
    Only the identity ``cell`` (see module docstring) is supported.
    """

    na: int
    n1: int
    n2: int
    n3: int
    basis: tuple[tuple[float, float, float], ...]

    def __post_init__(self) -> None:
        if self.na <= 0 or self.n1 <= 0 or self.n2 <= 0 or self.n3 <= 0:
            raise MalformedGeneratorInputError(
                f"na/n1/n2/n3 must all be positive, got na={self.na} "
                f"n1={self.n1} n2={self.n2} n3={self.n3}"
            )
        if len(self.basis) != self.na:
            raise MalformedGeneratorInputError(
                f"geometry declares na={self.na} but {len(self.basis)} basis "
                "coordinates were supplied"
            )
        for index, coordinate in enumerate(self.basis, start=1):
            if len(coordinate) != 3:
                raise MalformedGeneratorInputError(
                    f"basis coordinate {index} must have 3 components, got {coordinate!r}"
                )

    @property
    def natom(self) -> int:
        return self.na * self.n1 * self.n2 * self.n3

    def atom_index(self, i0: int, i1: int, i2: int, i3: int) -> int:
        """1-based global atom index, matching magnetizationinit.f90."""
        if not (1 <= i0 <= self.na):
            raise MalformedGeneratorInputError(f"i0={i0} out of range [1,{self.na}]")
        if not (0 <= i1 < self.n1 and 0 <= i2 < self.n2 and 0 <= i3 < self.n3):
            raise MalformedGeneratorInputError(
                f"cell index ({i1},{i2},{i3}) out of range for n1={self.n1} "
                f"n2={self.n2} n3={self.n3}"
            )
        return i0 + i1 * self.na + i2 * self.n1 * self.na + i3 * self.n2 * self.n1 * self.na

    def iter_atoms(self):
        """Yield ``(atom_index, i0, i1, i2, i3)`` in production loop order.

        The nesting order (i3 outermost, then i2, i1, i0 innermost) matches
        every ``Initmag`` loop in ``magnetizationinit.f90`` exactly, so a
        generator/parser that walks atoms in this order needs no additional
        sort step to agree with production numbering.
        """
        for i3 in range(self.n3):
            for i2 in range(self.n2):
                for i1 in range(self.n1):
                    for i0 in range(1, self.na + 1):
                        yield self.atom_index(i0, i1, i2, i3), i0, i1, i2, i3

    def cartesian(self, i0: int, i1: int, i2: int, i3: int) -> tuple[float, float, float]:
        """Cartesian position in alat units (identity cell only)."""
        bx, by, bz = self.basis[i0 - 1]
        return (i1 + bx, i2 + by, i3 + bz)


# ---------------------------------------------------------------------------
# Small math helpers
# ---------------------------------------------------------------------------


def _normalize(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(sum(component * component for component in vector))
    if norm <= 0.0:
        raise MalformedGeneratorInputError(f"cannot normalize a zero-norm vector {vector!r}")
    return (vector[0] / norm, vector[1] / norm, vector[2] / norm)


def _rotate(vector: tuple[float, float, float], axis: tuple[float, float, float],
            angle_rad: float) -> tuple[float, float, float]:
    """Rodrigues' rotation formula.

    Algebraically identical to the rotation-matrix construction in
    ``magnetizationinit.f90:469-486`` (verified by hand in the module
    docstring's analytic-justification derivation); this closed form is
    used directly rather than hand-transcribing the 3x3 matrix.
    """
    k = _normalize(axis)
    cos_a, sin_a = math.cos(angle_rad), math.sin(angle_rad)
    k_dot_v = k[0] * vector[0] + k[1] * vector[1] + k[2] * vector[2]
    k_cross_v = (
        k[1] * vector[2] - k[2] * vector[1],
        k[2] * vector[0] - k[0] * vector[2],
        k[0] * vector[1] - k[1] * vector[0],
    )
    return tuple(
        vector[c] * cos_a + k_cross_v[c] * sin_a + k[c] * k_dot_v * (1.0 - cos_a)
        for c in range(3)
    )


def content_hash(text: str) -> str:
    """Stable content hash for provenance/regeneration checks."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def _manifest(generator_name: str, geometry: Geometry, parameters: dict,
              generated_text: str, extra: dict | None = None) -> dict:
    record = {
        "generator_name": generator_name,
        "generator_version": GENERATOR_VERSION,
        "geometry": {
            "na": geometry.na, "n1": geometry.n1, "n2": geometry.n2, "n3": geometry.n3,
            "natom": geometry.natom,
            "basis": [list(coordinate) for coordinate in geometry.basis],
        },
        "parameters": parameters,
        "content_sha256": content_hash(generated_text),
    }
    if extra:
        record["extra"] = extra
    return record


def manifest_json(manifest: dict) -> str:
    """Canonical, deterministic JSON serialization of a manifest dict."""
    return json.dumps(manifest, sort_keys=True, indent=2) + "\n"


# ---------------------------------------------------------------------------
# Conical spiral / chiral +-q partner (Initmag=8 mechanism)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MomfileSpiralState:
    """Output of `conical_spiral_state`: a momfile plus Initmag=8 records."""

    momfile_text: str
    inpsd_records: dict[str, str]
    direction_by_atom: dict[int, tuple[float, float, float]]
    initial_torque_scale_hint: float
    manifest: dict


def _require_nondegenerate_cone(cone_angle_deg: float, turns: float,
                                 tolerance_deg: float) -> None:
    if turns == 0:
        raise DegenerateStateError(
            "turns=0 gives a zero propagation vector; a q=0 conical spiral is "
            "an exact stationary uniform state under isotropic exchange for "
            "any cone angle (see module docstring derivation)"
        )
    theta = cone_angle_deg % 180.0
    for special, description in ((0.0, "collinear/uniform"), (90.0, "planar spiral"),
                                  (180.0, "collinear/uniform")):
        if abs(theta - special) < tolerance_deg:
            raise DegenerateStateError(
                f"cone_angle_deg={cone_angle_deg} is within {tolerance_deg} degrees "
                f"of {special} degrees ({description}), an exact stationary "
                "eigenconfiguration under isotropic exchange with centrosymmetric "
                "neighbour shells (see module docstring derivation). Choose an "
                "angle strictly between the degenerate values to guarantee "
                "nonzero torque."
            )


def conical_spiral_state(geometry: Geometry, *, cone_angle_deg: float, turns: float,
                          axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
                          phase_deg: float = 0.0, modulation_cell_axis: int = 0,
                          moment_magnitude: float = 2.23, landeg: float = 2.0,
                          degeneracy_tolerance_deg: float = DEFAULT_DEGENERACY_TOLERANCE_DEG,
                          ) -> MomfileSpiralState:
    """Deterministic conical spin spiral via UppASD's Initmag=8 mechanism.

    ``turns`` is the number of full 2*pi rotations across the periodic
    supercell along ``modulation_cell_axis`` (0/1/2 for n1/n2/n3), matching
    the existing ``initmag_spin_spiral`` fixture's convention
    (``initpropvec 1/n1 0 0`` for one turn along n1). ``cone_angle_deg`` is
    the polar angle from ``axis`` (90 degrees is the planar/degenerate case
    flagged in RCG-04A; 0/180 are the uniform/degenerate case).
    ``phase_deg`` becomes ``initrotang`` directly.

    Raises ``DegenerateStateError`` for a cone angle or wavevector that
    would give an exactly stationary state (see module docstring).
    """
    if not (0.0 <= cone_angle_deg <= 180.0):
        raise MalformedGeneratorInputError(f"cone_angle_deg must be in [0,180], got {cone_angle_deg}")
    if modulation_cell_axis not in (0, 1, 2):
        raise MalformedGeneratorInputError(f"modulation_cell_axis must be 0, 1, or 2, got {modulation_cell_axis}")
    if moment_magnitude <= 0.0:
        raise MalformedGeneratorInputError(f"moment_magnitude must be positive, got {moment_magnitude}")
    _require_nondegenerate_cone(cone_angle_deg, turns, degeneracy_tolerance_deg)

    axis_unit = _normalize(axis)
    theta = math.radians(cone_angle_deg)
    # Base cone vector lies in the (x,z) half-plane at azimuth 0 before the
    # per-atom rotation is applied; combined with `axis_unit` this exactly
    # reproduces s(r) = (sin(theta)cos(q.x+phi), sin(theta)sin(q.x+phi),
    # cos(theta)) when axis_unit = (0,0,1) (the module docstring's
    # recommended normalized form), and generalizes to an arbitrary
    # rotation axis via the same Rodrigues construction production uses.
    base_vector = (math.sin(theta), 0.0, math.cos(theta))
    if abs(axis_unit[2]) < 1.0 - 1e-12:
        # Re-express the base vector relative to an arbitrary axis by
        # building an orthonormal frame (axis_unit, e1) and placing the
        # cone vector at polar angle theta from axis_unit, azimuth 0 in e1.
        reference = (0.0, 0.0, 1.0) if abs(axis_unit[0]) > 0.9 else (1.0, 0.0, 0.0)
        e1 = _normalize((
            reference[1] * axis_unit[2] - reference[2] * axis_unit[1],
            reference[2] * axis_unit[0] - reference[0] * axis_unit[2],
            reference[0] * axis_unit[1] - reference[1] * axis_unit[0],
        ))
        base_vector = tuple(
            math.sin(theta) * e1[c] + math.cos(theta) * axis_unit[c] for c in range(3)
        )
    momfile_lines = []
    for i0 in range(1, geometry.na + 1):
        momfile_lines.append(
            f"{i0} 1 {moment_magnitude:.17g} {base_vector[0]:.17g} "
            f"{base_vector[1]:.17g} {base_vector[2]:.17g} {landeg:.17g}"
        )
    momfile_text = "\n".join(momfile_lines) + "\n"

    axis_length = (geometry.n1, geometry.n2, geometry.n3)[modulation_cell_axis]
    propvec_fraction = [0.0, 0.0, 0.0]
    propvec_fraction[modulation_cell_axis] = turns / axis_length
    initpropvec = " ".join(f"{value:.17g}" for value in propvec_fraction)
    initrotvec = " ".join(f"{value:.17g}" for value in axis_unit)

    direction_by_atom: dict[int, tuple[float, float, float]] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        x, y, z = geometry.cartesian(i0, i1, i2, i3)
        modulation_coord = (x, y, z)[modulation_cell_axis]
        rotang = 2.0 * math.pi * (turns / axis_length) * modulation_coord + math.radians(phase_deg)
        direction_by_atom[atom] = _rotate(base_vector, axis_unit, rotang)

    parameters = {
        "cone_angle_deg": cone_angle_deg, "turns": turns, "axis": list(axis),
        "phase_deg": phase_deg, "modulation_cell_axis": modulation_cell_axis,
        "moment_magnitude": moment_magnitude, "landeg": landeg,
    }
    manifest = _manifest("conical_spiral_state", geometry, parameters, momfile_text)
    # sin(theta)*cos(theta) is the leading-order torque-scale factor from
    # the eigenvalue-gap derivation in the module docstring; it vanishes
    # exactly (not approximately) only at the rejected degenerate angles.
    torque_hint = abs(math.sin(theta) * math.cos(theta))
    return MomfileSpiralState(
        momfile_text=momfile_text,
        inpsd_records={
            "Initmag": "8", "initpropvec": initpropvec, "initrotvec": initrotvec,
            "initrotang": f"{phase_deg:.17g}",
        },
        direction_by_atom=direction_by_atom,
        initial_torque_scale_hint=torque_hint,
        manifest=manifest,
    )


def chiral_partner_pair(geometry: Geometry, *, cone_angle_deg: float, turns: float,
                         axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
                         phase_deg: float = 0.0, modulation_cell_axis: int = 0,
                         moment_magnitude: float = 2.23, landeg: float = 2.0,
                         degeneracy_tolerance_deg: float = DEFAULT_DEGENERACY_TOLERANCE_DEG,
                         ) -> tuple[MomfileSpiralState, MomfileSpiralState]:
    """Deterministic +-q chiral partner states for DMI-sensitive fixtures.

    Both states share every parameter except the sign of ``turns`` (hence
    of the propagation vector). Neither this function nor its caller reads
    any DMI file: the sign convention under test is not used to choose
    which of the two partners is "expected" to win (RCG-04H's job, using
    the RCG-02-accepted handedness oracle). Choosing ``turns`` away from a
    material's analytic DMI-energy-minimizing wavevector is the caller's
    responsibility when this is used for a DMI fixture.
    """
    plus = conical_spiral_state(
        geometry, cone_angle_deg=cone_angle_deg, turns=turns, axis=axis,
        phase_deg=phase_deg, modulation_cell_axis=modulation_cell_axis,
        moment_magnitude=moment_magnitude, landeg=landeg,
        degeneracy_tolerance_deg=degeneracy_tolerance_deg,
    )
    minus = conical_spiral_state(
        geometry, cone_angle_deg=cone_angle_deg, turns=-turns, axis=axis,
        phase_deg=phase_deg, modulation_cell_axis=modulation_cell_axis,
        moment_magnitude=moment_magnitude, landeg=landeg,
        degeneracy_tolerance_deg=degeneracy_tolerance_deg,
    )
    plus.manifest["extra"] = {"chirality": "+q"}
    minus.manifest["extra"] = {"chirality": "-q"}
    return plus, minus


# ---------------------------------------------------------------------------
# Periodic domain-wall pair (restart-file mechanism)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DomainWallPairState:
    """Output of `domain_wall_pair_state`: a restart-format text block."""

    restart_text: str
    direction_by_atom: dict[int, tuple[float, float, float]]
    wall_centers_cells: tuple[float, float]
    manifest: dict


def _wall_polar_angle(u: float, centers: tuple[float, float], width_cells: float) -> float:
    """Two-kink periodic soliton profile: 0 -> pi -> 0 across the period.

    theta(u) = 2*atan(exp((u-c1)/w)) - 2*atan(exp((u-c2)/w)); see the
    module-level derivation in the RCG-04B evidence document for why this
    closes periodically (up to exponentially small corrections in
    exp(-separation/width)) whenever both walls are separated from each
    other and from the periodic boundary by several widths.
    """
    c1, c2 = centers
    return 2.0 * math.atan(math.exp((u - c1) / width_cells)) - \
        2.0 * math.atan(math.exp((u - c2) / width_cells))


def domain_wall_pair_state(geometry: Geometry, *, axis_cell_index: int,
                            wall_centers_cells: tuple[float, float], width_cells: float,
                            easy_axis: tuple[float, float, float] = (0.0, 0.0, 1.0),
                            wall_type: str = "NEEL", chirality: int = 1,
                            cant_deg: float = 2.0, moment_magnitude: float = 2.23,
                            simid: str = "rcg04bwall", separation_margin_widths: float = 4.0,
                            ) -> DomainWallPairState:
    """Deterministic periodic domain-wall pair (two walls; PBC zero winding).

    Output is a full per-atom restart-format text block, matching
    ``source/System/restart.f90``'s exact reader/writer format (7-line
    header, then ``iter ens iatom |M| Mx My Mz`` per atom), because a
    localized wall profile cannot be expressed through the global
    ``Initmag=8`` propagation-vector mechanism used by the spiral
    generators above. AdaptiveCG production currently rejects the
    corresponding ``Initmag=4`` restart-load path (see module docstring);
    this output is therefore consumable by ordinary UppASD today but not
    yet by AdaptiveCG-enabled production.
    """
    axis_length = (geometry.n1, geometry.n2, geometry.n3)[axis_cell_index]
    if axis_cell_index not in (0, 1, 2):
        raise MalformedGeneratorInputError(f"axis_cell_index must be 0, 1, or 2, got {axis_cell_index}")
    if width_cells <= 0.0:
        raise MalformedGeneratorInputError(f"width_cells must be positive, got {width_cells}")
    if wall_type not in ("NEEL", "BLOCH"):
        raise MalformedGeneratorInputError(f"wall_type must be 'NEEL' or 'BLOCH', got {wall_type!r}")
    if chirality not in (1, -1):
        raise MalformedGeneratorInputError(f"chirality must be 1 or -1, got {chirality}")
    c1, c2 = wall_centers_cells
    if not (0.0 <= c1 < c2 < axis_length):
        raise MalformedGeneratorInputError(
            f"wall_centers_cells must satisfy 0 <= c1 < c2 < {axis_length} (axis length), "
            f"got {wall_centers_cells}"
        )
    gaps = (c1, c2 - c1, axis_length - c2)
    min_gap = min(gaps)
    if min_gap < separation_margin_widths * width_cells:
        raise DegenerateStateError(
            f"wall separations {gaps} (from 0, between centres, and back to the "
            f"periodic boundary) must each exceed {separation_margin_widths} x "
            f"width_cells={width_cells} ({separation_margin_widths * width_cells}) "
            "for the periodic two-kink profile to close continuously; got a "
            f"minimum separation of {min_gap}"
        )

    axis_unit = _normalize(easy_axis)
    # Neel walls rotate in the plane spanned by the wall-propagation
    # direction and the easy axis; Bloch walls rotate in the perpendicular
    # plane. Build e_neel as the propagation direction's component
    # perpendicular to the easy axis (equal to the propagation direction
    # itself whenever the two are already orthogonal, the case used by
    # every current fixture geometry), and e_bloch to complete a
    # right-handed orthonormal frame with axis_unit.
    propagation_axis = tuple(1.0 if c == axis_cell_index else 0.0 for c in range(3))
    parallel_component = sum(propagation_axis[c] * axis_unit[c] for c in range(3))
    perpendicular = tuple(
        propagation_axis[c] - parallel_component * axis_unit[c] for c in range(3)
    )
    perpendicular_norm = math.sqrt(sum(c * c for c in perpendicular))
    if perpendicular_norm < 1.0e-8:
        raise MalformedGeneratorInputError(
            f"easy_axis={easy_axis} is parallel to the wall propagation axis "
            f"(axis_cell_index={axis_cell_index}); a domain wall requires the "
            "easy axis to be transverse to its propagation direction"
        )
    e_neel = tuple(c / perpendicular_norm for c in perpendicular)
    e_bloch = (
        axis_unit[1] * e_neel[2] - axis_unit[2] * e_neel[1],
        axis_unit[2] * e_neel[0] - axis_unit[0] * e_neel[2],
        axis_unit[0] * e_neel[1] - axis_unit[1] * e_neel[0],
    )
    in_plane = e_neel if wall_type == "NEEL" else e_bloch
    cant_rad = math.radians(cant_deg)

    direction_by_atom: dict[int, tuple[float, float, float]] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        position = geometry.cartesian(i0, i1, i2, i3)
        u = position[axis_cell_index]
        theta = _wall_polar_angle(u, (c1, c2), width_cells)
        sin_t, cos_t = math.sin(theta), math.cos(theta)
        vector = tuple(
            chirality * sin_t * in_plane[c] + cos_t * axis_unit[c] for c in range(3)
        )
        # Small deterministic cant, breaking exact wall-plane symmetry
        # pinning: a fixed rotation of the local frame about axis_unit by
        # cant_deg, scaled by sin(theta) so the cant vanishes (rather than
        # diverges) in the uniform domains away from either wall core.
        cant_axis = e_bloch if wall_type == "NEEL" else e_neel
        vector = tuple(
            vector[c] + sin_t * math.sin(cant_rad) * cant_axis[c] for c in range(3)
        )
        direction_by_atom[atom] = _normalize(vector)

    lines = [
        "#" * 80,
        "# File type: R",
        "# Simulation type: S",
        f"# Number of atoms: {geometry.natom}",
        "# Number of ensembles: 1",
        "#" * 80,
        "#iter    ens   iatom        |Mom|            M_x            M_y            M_z",
    ]
    for atom in range(1, geometry.natom + 1):
        mx, my, mz = direction_by_atom[atom]
        lines.append(
            f"{0:8d}{1:8d}{atom:8d}  {moment_magnitude:16.8e}{mx:16.8e}{my:16.8e}{mz:16.8e}"
        )
    restart_text = "\n".join(lines) + "\n"

    parameters = {
        "axis_cell_index": axis_cell_index, "wall_centers_cells": list(wall_centers_cells),
        "width_cells": width_cells, "easy_axis": list(easy_axis), "wall_type": wall_type,
        "chirality": chirality, "cant_deg": cant_deg, "moment_magnitude": moment_magnitude,
        "simid": simid,
    }
    manifest = _manifest("domain_wall_pair_state", geometry, parameters, restart_text)
    return DomainWallPairState(
        restart_text=restart_text,
        direction_by_atom=direction_by_atom,
        wall_centers_cells=wall_centers_cells,
        manifest=manifest,
    )
