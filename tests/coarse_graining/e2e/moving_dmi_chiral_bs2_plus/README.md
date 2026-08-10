# RCG-04H: mixed fine/interface/coarse, +q, block_size_x=2 (refinement point)

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** the same physical 48-fine/32-interface/112-coarse-atom
partition as `../moving_dmi_chiral_bs1_plus`, discretised at
`block_size_x=2` (blocks 1-3 FINE of 12 total) instead of 1 -- a genuine
spatial-refinement pair at fixed physical evolution, matching
`../moving_static_mixed_bs2`'s role for RCG-04F.
