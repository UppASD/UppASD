# RCG-04H: mixed fine/interface/coarse, -q, accepted DMI sign

Fixture family for RCG-04H (DMI-HYBRID-CROSSING): a deterministic
mixed-resolution (fine/interface/coarse) DMI + uniaxial-anisotropy chirality
and dynamics claim, anchored to the accepted RCG-02 handedness convention.
See tests/coarse_graining/run_moving_dmi_chiral.py and
docs/RCG-04_MOVING_E2E_EVIDENCE.md (RCG-04H section) for the full
construction, oracle, and results.

**This fixture:** identical to `../moving_dmi_chiral_bs1_plus` (same
partition, same DMI file, same anisotropy) except the initial state is the
`-q` chiral partner (`inpsd.dat`'s `initpropvec` sign only; `momfile` is
byte-identical to every other RCG-04H fixture). For `D_zx=+0.02`, this
state is expected (RCG-02 convention) to have *lower* DMI energy and
*negative* signed chirality throughout, versus `../moving_dmi_chiral_bs1_plus`'s
higher energy and positive chirality.
