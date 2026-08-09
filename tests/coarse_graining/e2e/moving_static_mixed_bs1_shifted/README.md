# RCG-04F: E2E-MOVING-STATIC, block_size_x=1, interface shifted by half a period

Identical to `../moving_static_mixed_bs1` (same `momfile`, `block_size_x=1`,
`ncell 24 2 2`, physical/integration parameters) except the FINE seed is
blocks 13-18 instead of 1-6 -- a shift of 12 cells, exactly half the
`turns=1`/`n1=24` conical-spiral period (`initpropvec_x=1/24`), so the
interface now sits at the opposite phase of the modulation.

**Why this is expected to make (almost) no difference:** the conical-spiral
ansatz `s(x) = R_axis(q.x).s0` has, by construction (RCG-04B's own
degeneracy analysis), a local torque/field-misalignment magnitude that does
not depend on `x` for an exact uniform-pitch spiral on a translationally
invariant lattice -- the whole reason `MIN_MAX_TORQUE`/`MIN_RMS_TORQUE` are
gated with a *single* per-run floor rather than a per-atom one. A 12-cell
shift of the FINE/COARSE mask is, for this Hamiltonian and this state,
equivalent (up to a global `pi` rotation about the spiral axis, which
leaves every scalar/angular observable used in this slice invariant) to the
unshifted case. `run_moving_static_mixed.py` measures, rather than assumes,
that the two fixtures' spatial error profiles agree to within ordinary
run-to-run floating-point/integration reproducibility -- a genuine
regression check on the mask/geometry indexing (a bug there would *not*
respect this symmetry), not a vacuous re-run of the same case.
