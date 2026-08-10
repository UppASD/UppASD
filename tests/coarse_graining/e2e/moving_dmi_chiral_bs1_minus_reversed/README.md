# RCG-04H: negative control, -q, REVERSED DMI sign

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** identical to `../moving_dmi_chiral_bs1_minus` except it
reads `../dmfile_chiral_reversed` instead of `../dmfile_chiral`; see
`../moving_dmi_chiral_bs1_plus_reversed/README.md`.
