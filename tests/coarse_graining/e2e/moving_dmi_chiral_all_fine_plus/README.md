# RCG-04H: all-fine +q reference (DMI + uniaxial anisotropy)

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** all-fine (`cg_operator TENSOR`, STATIC, no mask -- every
block FINE) atomistic reference for the `+q` chiral partner, DMI
(`../dmfile_chiral`, `D_zx=+0.02`) and uniaxial anisotropy
(`../kfile_cg_x`) enabled. Used as the accuracy reference the mixed
fine/interface/coarse cases (`../moving_dmi_chiral_bs1_plus`,
`../moving_dmi_chiral_bs2_plus`) are compared against, exactly as
`moving_all_fine_wide` serves RCG-04E/F.
