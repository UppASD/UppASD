# RCG-04D E2E-MOVING-OFF-FINE: AdaptiveCG all-fine half

**Purpose:** the AdaptiveCG-enabled (`cg_operator TENSOR`, `cg_mask_mode
STATIC`, no mask file -> every block fine) half of the RCG-04D moving
feature-off/all-fine parity pair. Paired fixture: `../moving_feature_off`
(see that fixture's `README.md` for the shared initial-state construction,
generator parameters, integration parameters, and the production defect
found and fixed while building this pair). Driven together by
`tests/coarse_graining/run_moving_off_fine.py`. Full analysis and evidence:
`docs/RCG-04_MOVING_E2E_EVIDENCE.md` (RCG-04D section).

**Why not reuse `../static_all_fine`:** same reasoning as
`../moving_feature_off` vs `../feature_off` -- `static_all_fine` is a
uniform, zero-torque state and is retained unmodified as a smoke test; it
cannot support a moving-dynamics parity claim.

**`momfile` is byte-identical to `../moving_feature_off/momfile`** (verified
by `run_moving_off_fine.py` via content hash, not merely by convention: both
files are the direct, unmodified output of the same
`conical_spiral_state(...)` call recorded in `GENERATOR_MANIFEST.json`).

**`inpsd.dat` differs from `../moving_feature_off/inpsd.dat` in exactly the
intended AdaptiveCG block** (`simid`, `do_adaptive_cg`, `block_size_x/y/z`,
`cg_operator`, `cg_mask_mode`, `cg_diagnostics`) -- every other physically or
numerically relevant key (geometry, Hamiltonian inputs, `SDEalgh`,
`Initmag`/`initpropvec`/`initrotvec`/`initrotang`, `damping`, `timestep`,
`Nstep`, `do_tottraj`/`tottraj_step`) is identical. `run_moving_off_fine.py`
verifies this with a normalized key-by-key comparison of both `inpsd.dat`
files, not a visual/textual diff.

**Block topology:** `block_size_x/y/z = 1 2 2` on the shared `ncell 6 2 2`
host geometry gives 6 spatial blocks (not 12 -- `blocks_x = n1/block_size_x =
6`, `blocks_y = n2/block_size_y = 1`, `blocks_z = n3/block_size_z = 1`; the
production banner confirms `AdaptiveCG: atoms=48 blocks=6`), all fine
(`initial_fine=6`, `initial_coarse=0`, `active_atoms=48`). This is a
correction to an inventory-slice miscount in the RCG-04A evidence document
(section 3, which stated "12 blocks (4 atoms/block)" for this same host
geometry/block-size combination); noted here since this fixture's own
production banner is the direct evidence for the correct count, not fixed
retroactively in that earlier document per the governing rules' preference
against unrelated churn.

**The production defect found while building this fixture, and its fix,**
are documented in full in `../moving_feature_off/README.md` and the RCG-04D
evidence document: `adaptive_cg_cpu_step`
(`source/CoarseGraining/adaptivecgproduction.f90`) was missing the physical
gyromagnetic-ratio constant `gama` for the atomistic fine-region LLG update,
which silently froze this exact fixture (displacement ~1e-13 rad after 50
steps versus the paired feature-off run's ~0.124 rad) despite every setup
diagnostic (`capability accepted`, `active_atom_updates`, nonzero
`last_field_checksums_t`) looking normal. After the fix, this fixture's
trajectory agrees with `../moving_feature_off`'s to a maximum angular error
of ~5.8e-5 rad (~0.0033 degrees) out of ~0.124 rad (~7.1 degrees) of real
motion -- the residual is consistent with two independently implemented
Heun-family integrators (production's ordinary `evolve_first` vs. this
AdaptiveCG path's own predictor-corrector) applied to the same physics, not
a remaining defect.
