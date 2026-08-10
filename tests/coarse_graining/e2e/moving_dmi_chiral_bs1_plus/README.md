# RCG-04H: mixed fine/interface/coarse, +q, accepted DMI sign

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** the same fine/interface/coarse partition as
`../moving_static_mixed_bs1` (blocks 1-6 FINE of 24 `block_size_x=1`
blocks -> 48 fine / 32 interface / 112 coarse atoms), applied to the `+q`
chiral partner state with DMI (`../dmfile_chiral`, `D_zx=+0.02`) and
uniaxial anisotropy (`../kfile_cg_x`, easy axis `(1,0,0)`) both enabled.
Compared against `../moving_dmi_chiral_all_fine_plus` for complete-trajectory
accuracy, and against `../moving_dmi_chiral_bs1_minus` for the `+q`/`-q`
DMI energy ordering and signed-chirality claims.
