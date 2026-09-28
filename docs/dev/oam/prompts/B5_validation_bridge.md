# B5 — CI, pyswatter cross-check, bridge test  ·  B5.1–B5.3 Sonnet, B5.4 Opus  ·  gate G3

Blueprint B section **B5**. Paste the shared preamble first. Requires `[OAM-B4]`. B5.4 also requires `[OAM-A3]` and the GA decision.

## B5.1 — CI (Sonnet, commit `[OAM-B5a]`)

Register `run_traj_checks.py` in CTest with the label `oam-traj`, next to `oam-lswt`.

## B5.2 / B5.3 — pyswatter (Sonnet prepares, maintainer runs)

Prepare three case directories (`ell+1`, `shift_b`, `w_area`) and a shell script containing the exact `pyswatter-animate spin-oam-balance` commands from the blueprint, with the origin matched. Include the B5.3 sign case: counter-clockwise `m_x + i m_y` winding must give positive λ.

You may not have pyswatter. Don't emulate it; leave the script for the maintainer.

## B5.4 — Bridge test (Opus, commit `[OAM-B5b]`)

Before coding, write the derivation into `OAM_QUESTIONS.md`: how the Holstein–Primakoff amplitudes in `T_n(k0)` map onto site `m_x + i m_y` in the C1 frame, including the hole components, the `sqrt(2S)`, and the Fourier sign of `setup_ektij` (`exp(-i k·r)`). The mapping must reproduce the LSWT frequency when propagated. Check that first: the precession frequency of the packet has to match `E_n(k0)/ħ` to 1%.

Then follow the steps in Blueprint B5.4. Tolerances:
- `λ_centroid` within 5% of `F_n(k0)/ħ`, stable over the run;
- the l = 1 shift is +1 within 0.05.

A mismatch by a sign is the C13 question. Report it; don't fix it.

**Commit, report, stop for gate G3.**
