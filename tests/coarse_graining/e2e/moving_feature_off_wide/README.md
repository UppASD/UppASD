# RCG-04E: wide-geometry feature-off half (parity re-check before trusting the oracle)

**Purpose:** re-verify, at the `ncell 24 2 2` (192-atom) geometry RCG-04E
actually uses for its all-coarse block-size sweep, the same
feature-off/all-fine parity RCG-04D already established at `ncell 6 2 2`
(48 atoms). RCG-04D's parity finding exercises the exact same per-atom code
path (`adaptive_cg_cpu_step`'s `atomistic_atom` branch in
`source/CoarseGraining/adaptivecgproduction.f90`, using `gama*Landeg(atom)`
after RCG-04D's fix) regardless of atom count, so this is not expected to
reveal anything new -- but RCG-04E's all-coarse claims all depend on
`../moving_all_fine_wide` being a trustworthy atomistic oracle, and the
governing rules require that trust to be demonstrated, not assumed to
generalize from a different geometry size. Driven together with
`../moving_all_fine_wide` by `run_moving_all_coarse.py`'s
`verify_wide_reference_is_accepted()`.

**Construction:** identical mechanism to `../moving_feature_off`
(`docs/RCG-04_MOVING_E2E_EVIDENCE.md`, RCG-04D section), same
`conical_spiral_state(cone_angle_deg=40, turns=1, axis=(0,0,1),
modulation_cell_axis=0, moment_magnitude=2.23, landeg=2.0)` call, only
`ncell`/`initpropvec` differ (`24 2 2` / `turns/n1 = 1/24` instead of `6 2 2`
/ `1/6`). `momfile` is **byte-identical** to `../moving_feature_off/momfile`
and `../moving_all_fine/momfile` (verified by content hash in
`run_moving_all_coarse.py`) -- `conical_spiral_state`'s momfile output is a
per-basis-site template (`na=2` lines) and does not depend on `n1`/`n2`/`n3`
at all; only the `Initmag=8` `initpropvec` record (computed from `turns/n1`)
differs, exactly as expected.

Same Hamiltonian (exchange-only, shared `../jfile`), `damping 0.05`,
`timestep 1.0e-16`, `Nstep 50`, `do_tottraj Y`/`tottraj_step 5` (11 sampled
steps) as RCG-04D, so the RCG-04D `MAX_ANGULAR_ERROR_RAD`/
`MAX_COMPONENT_ERROR`/displacement-floor budgets apply unchanged.
`do_adaptive_cg` is unset (ordinary production path).
