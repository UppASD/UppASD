"""RCG-04D independent initial effective-field/torque oracle.

This module is deliberately narrow: it computes the discrete initial
Heisenberg-exchange effective field and LLG-driving torque for a periodic
lattice **independently of any UppASD production diagnostic**, purely from
(a) atom positions derived from ``moving_state_generator.Geometry`` (the same
position/atom-index formulas as production, but not production's own
neighbour-list code), (b) an exchange ``jfile``'s declared shell distances and
``Jij`` values, and (c) a caller-supplied per-atom initial direction map (for
RCG-04D, ``moving_state_generator.conical_spiral_state``'s own
independently-computed ``direction_by_atom``, never a value read back from a
production run). It is used only as a pre-acceptance nontriviality gate for a
moving fixture (RCG-04D governing rule: "compute the discrete initial
effective field and torque independently of the production diagnostic being
tested"), not as an oracle for the trajectory itself -- the trajectory
equivalence claim is proved by comparing the feature-off and all-fine
production runs against each other (see ``run_moving_off_fine.py``), not by
comparing either to this oracle.

Hamiltonian convention and its empirical calibration
------------------------------------------------------
Source reading (``source/Hamiltonian/hamiltonianactions.f90:204``,
``heisenberg_field``) shows the per-atom exchange effective field is a
direct, unweighted sum over the neighbour list, ``B_i = sum_j J_ij M_j``, no
extra sign flip. ``source/Hamiltonian/energy.f90:558`` (``update_ene``) then
computes ``-factor * (M_i . B_i)`` with ``factor=0.5`` for exchange
(``energy.f90:213``). Reading the source alone leaves two open questions
that matter for an independently *correct* oracle and are not safe to guess:
(1) whether ``M_i`` here still carries the ``mmom`` scalar or has effectively
been re-normalized upstream, and (2) whether a shell displacement whose
periodic image coincides with another shell displacement's image (a real
effect for this repository's small ``ncell 6 2 2`` / ``block_size 1 2 2``
host geometry -- e.g. the ``(1,0,0)`` shell's symmetry-generated ``(0,1,0)``
and ``(0,-1,0)`` partners both land on the *same* physical neighbour atom,
since ``n2=2`` makes a "+1 cell" and a "-1 cell" step in ``y`` identical)
should contribute once or twice.

Both questions were resolved by **running the real production binary** on
two disposable, non-tracked synthetic ``momfile``s built from this
repository's own shared ``tests/coarse_graining/e2e/{posfile,jfile}`` (not
by inspecting more source), and comparing this module's independently
computed neighbour list and energy formula against the printed
``localenergy.<simid>.out``/``totenergy.<simid>.out`` per-atom exchange
energy (``plotenergy 2``, ``Initmag 3``, ``Nstep 1``, no CG):

1. **Uniform aligned state** (every atom ``(1,0,0)``, ``mmom=2.23``):
   production prints ``Exc = -12.7251832`` per atom, for every one of the 48
   atoms identically. This module's ``build_geometric_bonds`` (single-count
   Cartesian-distance matching against the jfile shells, one bond per
   physical neighbour atom -- **not** a naive symmetry-orbit count, which
   would double-count the coinciding-image case above) gives each atom
   exactly 25 neighbours (shell-by-shell: 8 at distance 0.866, 4 at 1.0, 5 at
   1.414, 8 at 1.658). ``-sum(Jij for the 25 neighbours) = -12.72518322837978``,
   matching production to all 9 printed significant figures.
2. **Re-running with every ``mmom`` changed from 2.23 to 1.0** (uniform
   state unchanged otherwise) reproduces the *identical* ``-12.7251832``,
   confirming the printed exchange energy here does not scale with ``mmom``
   -- i.e. the convention that reconciles with production is the plain
   **unit-vector** form ``E_i = -sum_j J_ij (S_i . S_j)`` (factor 1, not
   0.5; source-level ``mmom``/factor-of-two details evidently cancel
   upstream of this diagnostic and are not re-derived here since the
   externally observable convention is what an independent oracle needs).
3. **Sublattice flip** (every basis-site-1 atom set to ``(0,0,1)``, every
   basis-site-2 atom left at ``(1,0,0)`` -- basis-site 1/2 are exactly the
   two interpenetrating simple-cubic sublattices of this jfile's BCC
   structure, so this flips only the *inter*-sublattice bonds' alignment):
   production prints ``Exc = -2.72937118`` per atom, exactly matching this
   module's prediction restricted to each atom's same-sublattice (intra)
   neighbours only (``-(4*J[1.0] + 5*J[1.414]) = -2.7293711796337``, using
   the identical 25-neighbour list from point 1). This independently
   confirms the neighbour *topology* (not just its size) atom-class by
   atom-class, not merely the aggregate count.

Both checks match to 9 significant figures (not an approximate/near match),
across two independent, structurally different configurations, which is
sufficient confidence to adopt ``E_i = -sum_j J_ij (S_i . S_j)``,
``B_i = sum_j J_ij S_j`` (unit vectors, no ``mmom``/Landeg factor) as this
module's calibrated convention. This calibration used only the existing
*uniform, stationary* fixture inputs (the same ground state
``feature_off``/``static_all_fine`` already use) -- never the RCG-04D moving
state itself -- so the moving-state torque computed below remains
independent of the diagnostic it gates, per the governing rule.

Because the torque this module reports (``S_i x B_i``) omits only a single,
atom-independent, strictly positive physical prefactor (the LLG gyromagnetic
factor, and here also the omitted ``mmom`` scale, which is uniform across
every atom in every RCG-04D fixture), neither omission can affect the
nonzero-torque nontriviality claim or the *relative* max/RMS comparison
across atoms -- the only things this module is used to gate.

Neighbour-list construction: this module does **not** read or trust
production's ``ham%nlist``/``ham%ncoup`` (that would make the oracle
circular). It builds bonds purely from Cartesian distance: for every ordered
atom pair, the minimum-image displacement (periodic box lengths
``n1*alat, n2*alat, n3*alat``) is compared against each declared jfile shell
distance, and a bond is added when they match within tolerance -- one bond
per physical neighbour atom, exactly the convention calibration point 1
above validated. This is only safe when every shell magnitude projects to at
most half of the corresponding box length in each Cartesian direction
(``validate_minimum_image_safe`` checks this explicitly and raises
otherwise -- see that function's docstring).
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from moving_state_generator import Geometry

Vector = tuple[float, float, float]


class MalformedJfileError(ValueError):
    """Raised when jfile text does not match the expected 6-column shell format."""


class UnsafeMinimumImageError(ValueError):
    """Raised when a shell distance is too large for plain minimum-image PBC."""


@dataclass(frozen=True)
class ExchangeShell:
    """One exchange shell: a bond vector magnitude and its coupling."""

    dx: float
    dy: float
    dz: float
    jij: float

    @property
    def distance(self) -> float:
        return math.sqrt(self.dx * self.dx + self.dy * self.dy + self.dz * self.dz)


def parse_jfile_shells(text: str) -> tuple[ExchangeShell, ...]:
    """Parse a UppASD ``jfile`` (``isite jsite rx ry rz Jij`` per line).

    Only the geometric displacement and ``Jij`` columns are used; ``isite``/
    ``jsite`` are not interpreted (every fixture jfile under
    ``tests/coarse_graining/e2e`` uses a single chemical type).
    """
    shells: list[ExchangeShell] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != 6:
            raise MalformedJfileError(
                f"jfile line {line_number}: expected 6 fields, got {len(fields)}: {line!r}"
            )
        try:
            int(fields[0])  # isite: not interpreted, validated as integer only
            int(fields[1])  # jsite: not interpreted, validated as integer only
            dx, dy, dz, jij = (float(value) for value in fields[2:6])
        except ValueError as exc:
            raise MalformedJfileError(f"jfile line {line_number}: {exc}") from exc
        shells.append(ExchangeShell(dx=dx, dy=dy, dz=dz, jij=jij))
    if not shells:
        raise MalformedJfileError("jfile text contained no shell records")
    return tuple(shells)


def validate_minimum_image_safe(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                                 alat: float = 1.0) -> None:
    """Raise unless every shell is safely inside plain minimum-image PBC.

    Plain minimum-image PBC (wrap each Cartesian component to the nearest
    periodic replica) enumerates each physical bond exactly once, without
    ambiguity or double-counting, only when every shell's per-axis
    displacement magnitude is at most half that axis's box length. Every
    shell under ``tests/coarse_graining/e2e`` satisfies this (checked here,
    not assumed) for the standard ``ncell 6 2 2`` host geometry; see the
    module docstring's calibration section for the case (a shell component
    landing exactly on the half-box boundary) that makes this the correct
    single-count convention rather than a symmetry-orbit count.
    """
    box = (geometry.n1 * alat, geometry.n2 * alat, geometry.n3 * alat)
    for shell in shells:
        for component, length in zip((shell.dx, shell.dy, shell.dz), box):
            if abs(component) > length / 2.0 + 1.0e-9:
                raise UnsafeMinimumImageError(
                    f"shell component {component} exceeds half the box length "
                    f"{length / 2.0} -- plain minimum-image PBC is not safe for this "
                    "geometry/jfile combination; a full periodic-image sum would be "
                    "required instead"
                )


def _minimum_image(displacement: Vector, box: Vector) -> Vector:
    return tuple(
        d - round(d / length) * length for d, length in zip(displacement, box)
    )


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


def build_geometric_bonds(geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
                           alat: float = 1.0, distance_tolerance: float = 1.0e-6,
                           ) -> dict[int, list[tuple[int, float]]]:
    """Independent geometric neighbour list: ``{atom: [(neighbour, Jij), ...]}``.

    Built purely from Cartesian distance matching against ``shells``, never
    from production's own neighbour-list construction: one bond per physical
    neighbour atom (see module docstring calibration point 1 for why this,
    rather than a symmetry-orbit multiplicity count, is the convention that
    matches real production output for this host geometry). Raises
    ``UnsafeMinimumImageError`` first (via ``validate_minimum_image_safe``)
    if the geometry/jfile combination is not safe for this method.
    """
    validate_minimum_image_safe(geometry, shells, alat=alat)
    box = (geometry.n1 * alat, geometry.n2 * alat, geometry.n3 * alat)
    positions: dict[int, Vector] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        cx, cy, cz = geometry.cartesian(i0, i1, i2, i3)
        positions[atom] = (cx * alat, cy * alat, cz * alat)

    shell_by_distance = [(shell.distance, shell.jij) for shell in shells]
    bonds: dict[int, list[tuple[int, float]]] = {atom: [] for atom in positions}
    natom = geometry.natom
    for i in range(1, natom + 1):
        pos_i = positions[i]
        for j in range(1, natom + 1):
            if i == j:
                continue
            raw = tuple(positions[j][c] - pos_i[c] for c in range(3))
            wrapped = _minimum_image(raw, box)
            distance = _norm(wrapped)
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


@dataclass(frozen=True)
class TorqueReport:
    """Per-atom exchange field/torque, plus summary nontriviality metrics."""

    field_by_atom: dict[int, Vector]
    torque_by_atom: dict[int, Vector]
    max_torque: float
    rms_torque: float
    max_field_misalignment_deg: float

    def as_dict(self) -> dict:
        return {
            "max_torque": self.max_torque,
            "rms_torque": self.rms_torque,
            "max_field_misalignment_deg": self.max_field_misalignment_deg,
            "natom": len(self.torque_by_atom),
        }


def initial_torque_report(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                           direction_by_atom: dict[int, Vector], *,
                           alat: float = 1.0,
                           distance_tolerance: float = 1.0e-6) -> TorqueReport:
    """Independent initial exchange field/torque for every atom.

    ``B_i = sum_{j in neighbours(i)} J_ij * S_j`` and ``torque_i = S_i x B_i``
    (unit vectors throughout -- see module docstring's calibration section
    for why this, not a moment-weighted ``mmom_i * mmom_j`` form, is the
    convention that reproduces real production exchange-energy output).
    Reported up to the atom-independent, strictly positive LLG gyromagnetic
    prefactor, which is omitted throughout: it scales every atom's torque
    identically and so cannot affect either the nonzero-torque nontriviality
    claim or the *relative* max/RMS comparison.
    """
    bonds = build_geometric_bonds(
        geometry, shells, alat=alat, distance_tolerance=distance_tolerance,
    )
    field_by_atom: dict[int, Vector] = {}
    torque_by_atom: dict[int, Vector] = {}
    misalignment_degrees: list[float] = []
    for atom, neighbours in bonds.items():
        field = [0.0, 0.0, 0.0]
        for neighbour, jij in neighbours:
            neighbour_direction = direction_by_atom[neighbour]
            for component in range(3):
                field[component] += jij * neighbour_direction[component]
        field_vector = (field[0], field[1], field[2])
        field_by_atom[atom] = field_vector
        torque_by_atom[atom] = _cross(direction_by_atom[atom], field_vector)
        field_norm = _norm(field_vector)
        if field_norm > 0.0:
            cos_angle = max(-1.0, min(1.0,
                _dot(direction_by_atom[atom], field_vector) / field_norm))
            misalignment_degrees.append(math.degrees(math.acos(cos_angle)))
    torque_magnitudes = [_norm(torque) for torque in torque_by_atom.values()]
    max_torque = max(torque_magnitudes) if torque_magnitudes else 0.0
    rms_torque = (
        math.sqrt(sum(t * t for t in torque_magnitudes) / len(torque_magnitudes))
        if torque_magnitudes else 0.0
    )
    return TorqueReport(
        field_by_atom=field_by_atom,
        torque_by_atom=torque_by_atom,
        max_torque=max_torque,
        rms_torque=rms_torque,
        max_field_misalignment_deg=max(misalignment_degrees) if misalignment_degrees else 0.0,
    )


def exchange_energy_per_atom(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                              direction_by_atom: dict[int, Vector], *,
                              alat: float = 1.0,
                              distance_tolerance: float = 1.0e-6) -> dict[int, float]:
    """Per-atom exchange energy ``E_i = -sum_j J_ij (S_i . S_j)`` (calibrated form).

    Used only for calibrating/self-testing this module against real
    production ``localenergy.<simid>.out`` output (see module docstring);
    not used by the RCG-04D moving-state acceptance gate itself, which only
    needs the torque, not the energy.
    """
    bonds = build_geometric_bonds(
        geometry, shells, alat=alat, distance_tolerance=distance_tolerance,
    )
    energies: dict[int, float] = {}
    for atom, neighbours in bonds.items():
        total = 0.0
        for neighbour, jij in neighbours:
            total += jij * _dot(direction_by_atom[atom], direction_by_atom[neighbour])
        energies[atom] = -total
    return energies


def coordination_numbers(geometry: Geometry, shells: tuple[ExchangeShell, ...], *,
                          alat: float = 1.0, distance_tolerance: float = 1.0e-6,
                          ) -> dict[int, int]:
    """Neighbour count per atom -- an independent geometric sanity check."""
    bonds = build_geometric_bonds(
        geometry, shells, alat=alat, distance_tolerance=distance_tolerance,
    )
    return {atom: len(neighbours) for atom, neighbours in bonds.items()}


# ---------------------------------------------------------------------------
# RCG-04G: uniaxial anisotropy field/torque extension
# ---------------------------------------------------------------------------
#
# Formula read directly from the production code path this fixture actually
# exercises -- ``source/CoarseGraining/adaptivecgproduction.f90:
# evaluate_atomistic_anisotropy`` (the atomistic-onsite term used by
# AdaptiveCG's own energy evaluator, not the legacy Hamiltonian this module's
# exchange convention above was calibrated against) -- rather than guessed:
#
#   c = S_i . e_axis
#   E_i = k1*c**2 + k2*(2*c**2 - c**4)
#   field_i = -(2*c*(k1 + 2*k2*(1-c**2))) * e_axis / (mub * moment_mub_i)
#
# For every RCG-04G fixture, k2=0 and moment_mub_i is uniform across atoms
# (moment_magnitude passed to domain_wall_pair_state), so the reported
# torque below omits the common positive prefactor 1/(mub*moment_mub) --
# exactly the same "uniform, positive, atom-independent" omission this
# module's docstring already justifies for the exchange field's gyromagnetic
# prefactor: it cannot change a zero-vs-nonzero torque conclusion or the
# relative comparison across atoms, only its absolute scale.


def anisotropy_field(direction: Vector, axis: Vector, k1: float, k2: float = 0.0) -> Vector:
    """Uniaxial anisotropy field for one atom, in the calibrated production convention."""
    axis_norm = _norm(axis)
    axis_unit = tuple(component / axis_norm for component in axis)
    c = _dot(direction, axis_unit)
    derivative = 2.0 * c * (k1 + 2.0 * k2 * (1.0 - c * c))
    return tuple(-derivative * component for component in axis_unit)


def combined_torque_report(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                            direction_by_atom: dict[int, Vector], *,
                            anisotropy_axis: Vector, anisotropy_k1: dict[int, float],
                            anisotropy_k2: dict[int, float] | None = None,
                            alat: float = 1.0,
                            distance_tolerance: float = 1.0e-6) -> TorqueReport:
    """Exchange-plus-uniaxial-anisotropy initial torque, independent of production.

    ``anisotropy_k1``/``anisotropy_k2`` are keyed by basis site id (1-based,
    matching ``Geometry.iter_atoms``'s ``i0``), matching how ``kfile_cg``
    assigns one (k1, k2) pair per basis site. Used as the RCG-04G
    pre-acceptance nontriviality gate for the domain-wall-pair fixture, the
    same role RCG-04D's exchange-only ``initial_torque_report`` plays for
    the conical-spiral fixture -- extended here because a domain wall's
    existence depends on the exchange/anisotropy competition, not exchange
    alone.
    """
    k2_by_site = anisotropy_k2 or {}
    bonds = build_geometric_bonds(
        geometry, shells, alat=alat, distance_tolerance=distance_tolerance,
    )
    field_by_atom: dict[int, Vector] = {}
    torque_by_atom: dict[int, Vector] = {}
    misalignment_degrees: list[float] = []
    atom_basis: dict[int, int] = {atom: i0 for atom, i0, _, _, _ in geometry.iter_atoms()}
    for atom, neighbours in bonds.items():
        field = [0.0, 0.0, 0.0]
        for neighbour, jij in neighbours:
            neighbour_direction = direction_by_atom[neighbour]
            for component in range(3):
                field[component] += jij * neighbour_direction[component]
        basis = atom_basis[atom]
        aniso = anisotropy_field(
            direction_by_atom[atom], anisotropy_axis,
            anisotropy_k1[basis], k2_by_site.get(basis, 0.0),
        )
        field_vector = (field[0] + aniso[0], field[1] + aniso[1], field[2] + aniso[2])
        field_by_atom[atom] = field_vector
        torque_by_atom[atom] = _cross(direction_by_atom[atom], field_vector)
        field_norm = _norm(field_vector)
        if field_norm > 0.0:
            cos_angle = max(-1.0, min(1.0,
                _dot(direction_by_atom[atom], field_vector) / field_norm))
            misalignment_degrees.append(math.degrees(math.acos(cos_angle)))
    torque_magnitudes = [_norm(torque) for torque in torque_by_atom.values()]
    max_torque = max(torque_magnitudes) if torque_magnitudes else 0.0
    rms_torque = (
        math.sqrt(sum(t * t for t in torque_magnitudes) / len(torque_magnitudes))
        if torque_magnitudes else 0.0
    )
    return TorqueReport(
        field_by_atom=field_by_atom,
        torque_by_atom=torque_by_atom,
        max_torque=max_torque,
        rms_torque=rms_torque,
        max_field_misalignment_deg=max(misalignment_degrees) if misalignment_degrees else 0.0,
    )


# ---------------------------------------------------------------------------
# RCG-04H: independent DMI field/energy oracle, anchored to the RCG-02
# accepted handedness convention
# ---------------------------------------------------------------------------
#
# Convention, restated from docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md (the
# governing, already-accepted source): for a directed neighbour list
# containing *both* orientations of every physical bond (``D_ji = -D_ij``,
# supplied explicitly by the dmfile, never inferred by this module),
#
#   E_D = (mu_B/2) * sum_i sum_{j in N_D(i)} D_ij . (M_i x M_j)
#   B_i = -(1/mu_B) dE_D/dM_i = sum_{j in N_D(i)} D_ij x M_j
#
# For the accepted dimer M_1=+x, M_2=+y, D_12=+D*z, D_21=-D*z: E_D=mu_B*D,
# B_1=-D*x, B_2=-D*y (RCG-02's own worked example -- reproduced exactly by
# this module's functions below and cross-checked against it in
# test_torque_oracle.py, independent of this module's DMI file under test).
# Positive D_zx (a D-vector along +z on the +x-displacement bond) therefore
# raises the right-handed "+q" planar chain m=(cos(qx),sin(qx),0) and favours
# its left-handed "-q" partner -- RCG-02's own stated result, reproduced
# here as ``dmi_energy_reduced`` on ``moving_state_generator.chiral_partner_pair``
# output, never derived from whichever dmfile a caller is testing.
#
# Following this module's existing calibrated convention (unit-vector
# ``S``, no ``mu_B``/``mmom`` prefactor -- see the module docstring), the
# reported quantities below omit exactly the same uniform, strictly
# positive factor as ``initial_torque_report``/``exchange_energy_per_atom``:
# ``field = sum_j D_ij x S_j`` and ``E_D_reduced = sum_i sum_j D_ij.(S_i x S_j) / 2``.
# This cannot affect a sign, an energy *ordering*, or a zero-vs-nonzero
# conclusion -- only the two things this module is used for.


@dataclass(frozen=True)
class DmiBond:
    """One directed DMI shell entry: displacement vector and its D vector."""

    dx: float
    dy: float
    dz: float
    d: Vector


def parse_dmfile_bonds(text: str) -> tuple[DmiBond, ...]:
    """Parse a UppASD ``dmfile`` (``isite jsite rx ry rz Dx Dy Dz`` per line).

    Only the geometric displacement and ``D`` columns are used; ``isite``/
    ``jsite`` are not interpreted, matching ``parse_jfile_shells``. Unlike
    exchange shells, both orientations of a physical bond must appear as
    separate lines with their own (generally different) ``D`` vector -- this
    parser does not synthesize a reciprocal entry, matching the accepted
    RCG-02 convention that the directed list is supplied complete.
    """
    bonds: list[DmiBond] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) != 8:
            raise MalformedJfileError(
                f"dmfile line {line_number}: expected 8 fields, got {len(fields)}: {line!r}"
            )
        try:
            int(fields[0])
            int(fields[1])
            dx, dy, dz, ddx, ddy, ddz = (float(value) for value in fields[2:8])
        except ValueError as exc:
            raise MalformedJfileError(f"dmfile line {line_number}: {exc}") from exc
        bonds.append(DmiBond(dx=dx, dy=dy, dz=dz, d=(ddx, ddy, ddz)))
    if not bonds:
        raise MalformedJfileError("dmfile text contained no shell records")
    return tuple(bonds)


def negate_dmi_bonds(dmi_bonds: tuple[DmiBond, ...]) -> tuple[DmiBond, ...]:
    """Return the sign-reversed DMI operator: every ``D`` vector negated.

    Used only to build the RCG-04H negative control's disposable
    sign-reversed oracle input from the accepted ``dmfile_chiral`` text --
    never to regenerate an "expected" chirality from whatever sign happens
    to be under test (the fixed oracle in ``run_moving_dmi_chiral.py`` is
    derived once, from the accepted-sign convention, and reused unchanged
    against this reversed operator).
    """
    return tuple(DmiBond(dx=b.dx, dy=b.dy, dz=b.dz, d=(-b.d[0], -b.d[1], -b.d[2])) for b in dmi_bonds)


def build_directed_dmi_bonds(geometry: Geometry, dmi_bonds: tuple[DmiBond, ...], *,
                              alat: float = 1.0, distance_tolerance: float = 1.0e-6,
                              ) -> dict[int, list[tuple[int, Vector]]]:
    """Independent directed DMI neighbour list: ``{atom: [(neighbour, D_ij), ...]}``.

    Unlike ``build_geometric_bonds`` (exchange), matching is by the full
    displacement *vector*, not just its magnitude, because ``D`` is
    direction-dependent (the accepted convention gives the ``+x`` and
    ``-x`` bonds of the same pair opposite ``D``, both of which must be
    supplied as separate ``dmi_bonds`` entries and matched to the correct
    physical direction).
    """
    box = (geometry.n1 * alat, geometry.n2 * alat, geometry.n3 * alat)
    positions: dict[int, Vector] = {}
    for atom, i0, i1, i2, i3 in geometry.iter_atoms():
        cx, cy, cz = geometry.cartesian(i0, i1, i2, i3)
        positions[atom] = (cx * alat, cy * alat, cz * alat)

    bonds: dict[int, list[tuple[int, Vector]]] = {atom: [] for atom in positions}
    natom = geometry.natom
    for i in range(1, natom + 1):
        pos_i = positions[i]
        for j in range(1, natom + 1):
            if i == j:
                continue
            raw = tuple(positions[j][c] - pos_i[c] for c in range(3))
            wrapped = _minimum_image(raw, box)
            for bond in dmi_bonds:
                delta = (
                    wrapped[0] - bond.dx * alat,
                    wrapped[1] - bond.dy * alat,
                    wrapped[2] - bond.dz * alat,
                )
                if _norm(delta) <= distance_tolerance:
                    bonds[i].append((j, bond.d))
    return bonds


def dmi_field(geometry: Geometry, dmi_bonds: tuple[DmiBond, ...],
              direction_by_atom: dict[int, Vector], *, alat: float = 1.0,
              distance_tolerance: float = 1.0e-6) -> dict[int, Vector]:
    """Independent initial DMI effective field, ``B_i = sum_j D_ij x S_j``.

    Accepted RCG-02 convention (see module section docstring above),
    calibrated-unit-vector form (no ``mu_B``/``mmom`` prefactor).
    """
    directed = build_directed_dmi_bonds(
        geometry, dmi_bonds, alat=alat, distance_tolerance=distance_tolerance,
    )
    field_by_atom: dict[int, Vector] = {}
    for atom, neighbours in directed.items():
        field = [0.0, 0.0, 0.0]
        for neighbour, d_ij in neighbours:
            cross = _cross(d_ij, direction_by_atom[neighbour])
            for component in range(3):
                field[component] += cross[component]
        field_by_atom[atom] = (field[0], field[1], field[2])
    return field_by_atom


def dmi_energy_reduced(geometry: Geometry, dmi_bonds: tuple[DmiBond, ...],
                        direction_by_atom: dict[int, Vector], *, alat: float = 1.0,
                        distance_tolerance: float = 1.0e-6) -> float:
    """Independent total DMI energy, ``E_D = sum_i sum_j D_ij.(S_i x S_j) / 2``.

    Accepted RCG-02 sign convention, calibrated-unit-vector ("reduced")
    form -- see module section docstring above for the omitted, uniform,
    strictly positive ``mu_B*mmom_i*mmom_j`` factor and why it cannot
    affect a sign or an energy *ordering*. The ``1/2`` matches RCG-02's own
    stated factor, which removes the duplicate from the directed
    (both-orientations) list.
    """
    directed = build_directed_dmi_bonds(
        geometry, dmi_bonds, alat=alat, distance_tolerance=distance_tolerance,
    )
    total = 0.0
    for atom, neighbours in directed.items():
        for neighbour, d_ij in neighbours:
            total += _dot(d_ij, _cross(direction_by_atom[atom], direction_by_atom[neighbour]))
    return total / 2.0


def anisotropy_energy_reduced(geometry: Geometry, direction_by_atom: dict[int, Vector], *,
                               axis: Vector, k1_by_site: dict[int, float],
                               k2_by_site: dict[int, float] | None = None) -> float:
    """Independent total uniaxial anisotropy energy, calibrated-reduced form.

    ``E_i = k1*c**2 + k2*(2*c**2 - c**4)``, ``c = S_i . e_axis`` -- the same
    per-atom formula (and the same omitted ``mu_B*moment_mub**2`` prefactor)
    ``anisotropy_field`` above already uses; see that function's docstring
    for the source (``adaptivecgproduction.f90:evaluate_atomistic_anisotropy``).
    """
    k2 = k2_by_site or {}
    axis_norm = _norm(axis)
    axis_unit = tuple(component / axis_norm for component in axis)
    atom_basis: dict[int, int] = {atom: i0 for atom, i0, _, _, _ in geometry.iter_atoms()}
    total = 0.0
    for atom, direction in direction_by_atom.items():
        basis = atom_basis[atom]
        c = _dot(direction, axis_unit)
        total += k1_by_site[basis] * c * c + k2.get(basis, 0.0) * (2.0 * c * c - c ** 4)
    return total


def dmi_anisotropy_torque_report(geometry: Geometry, shells: tuple[ExchangeShell, ...],
                                  dmi_bonds: tuple[DmiBond, ...],
                                  direction_by_atom: dict[int, Vector], *,
                                  anisotropy_axis: Vector, anisotropy_k1: dict[int, float],
                                  anisotropy_k2: dict[int, float] | None = None,
                                  alat: float = 1.0,
                                  distance_tolerance: float = 1.0e-6) -> TorqueReport:
    """Exchange + DMI + uniaxial-anisotropy initial torque, independent of production.

    Same role as ``combined_torque_report`` (RCG-04G), extended with the DMI
    contribution above -- the RCG-04H pre-acceptance nontriviality gate.
    """
    exchange = combined_torque_report(
        geometry, shells, direction_by_atom, anisotropy_axis=anisotropy_axis,
        anisotropy_k1=anisotropy_k1, anisotropy_k2=anisotropy_k2, alat=alat,
        distance_tolerance=distance_tolerance,
    )
    dmi = dmi_field(geometry, dmi_bonds, direction_by_atom, alat=alat,
                     distance_tolerance=distance_tolerance)
    field_by_atom: dict[int, Vector] = {}
    torque_by_atom: dict[int, Vector] = {}
    misalignment_degrees: list[float] = []
    for atom, direction in direction_by_atom.items():
        exchange_field = exchange.field_by_atom[atom]
        dmi_field_atom = dmi.get(atom, (0.0, 0.0, 0.0))
        field_vector = tuple(exchange_field[c] + dmi_field_atom[c] for c in range(3))
        field_by_atom[atom] = field_vector
        torque_by_atom[atom] = _cross(direction, field_vector)
        field_norm = _norm(field_vector)
        if field_norm > 0.0:
            cos_angle = max(-1.0, min(1.0, _dot(direction, field_vector) / field_norm))
            misalignment_degrees.append(math.degrees(math.acos(cos_angle)))
    torque_magnitudes = [_norm(torque) for torque in torque_by_atom.values()]
    max_torque = max(torque_magnitudes) if torque_magnitudes else 0.0
    rms_torque = (
        math.sqrt(sum(t * t for t in torque_magnitudes) / len(torque_magnitudes))
        if torque_magnitudes else 0.0
    )
    return TorqueReport(
        field_by_atom=field_by_atom,
        torque_by_atom=torque_by_atom,
        max_torque=max_torque,
        rms_torque=rms_torque,
        max_field_misalignment_deg=max(misalignment_degrees) if misalignment_degrees else 0.0,
    )
