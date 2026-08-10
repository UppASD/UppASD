# RCG-04H: negative control, +q, REVERSED DMI sign

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** identical to `../moving_dmi_chiral_bs1_plus` except it
reads `../dmfile_chiral_reversed` (`D_zx=-0.02`, the exact `D -> -D`
negation of `../dmfile_chiral`, verified by
`run_moving_dmi_chiral.py:verify_reversed_dmfile_is_exact_negation`)
instead of `../dmfile_chiral`. Used only for the RCG-04H sign
negative control: the accepted-sign oracle (derived once from `D_zx>0`)
must fail against this run, evaluated as a genuinely different, always-
tracked interaction file rather than a temporary source mutation.
