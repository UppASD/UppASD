"""Measurement profiles (WP-03 section E).

A profile turns into ordinary, allow-listed `inpsd.dat` overrides -- never a
change to production source. Both profiles reuse the measurement-cadence
keys `cases.GLOBALLY_SAFE_OVERRIDE_KEYS` already exposes
(`avrg_step`/`cumu_step`/`tottraj_step`/`ene_step`); a case only supports a
profile if it declares the relevant keys in its own
`allowed_input_overrides` (`cases.build_overrides` raises otherwise).

DYNAMICS_ONLY pushes every measurement cadence past the run's own `nstep`,
so the optional measurement/output work such a cadence gates never triggers
within the measured run -- the run still performs every timestep the chosen
dynamics and Hamiltonian inherently require.

PRODUCTION_LIGHT uses one fixed, explicitly documented cadence, representative
of a real production run that wants periodic monitoring rather than none.
"""

from __future__ import annotations

DYNAMICS_ONLY = "DYNAMICS_ONLY"
PRODUCTION_LIGHT = "PRODUCTION_LIGHT"

PROFILES = (DYNAMICS_ONLY, PRODUCTION_LIGHT)

# The representative measurement cadence PRODUCTION_LIGHT samples use.
# Documented here rather than left implicit, per blueprint section 7: "Use a
# fixed, documented, representative measurement cadence."
PRODUCTION_LIGHT_CADENCE = {
    "avrg_step": 100,
    "cumu_step": 100,
    "tottraj_step": 100,
    "ene_step": 100,
}

_CADENCE_KEYS = tuple(PRODUCTION_LIGHT_CADENCE)


class MeasurementProfileError(ValueError):
    """An unknown measurement profile was requested."""


def profile_overrides(profile, nstep):
    """Return the `inpsd.dat` overrides implementing ``profile`` at ``nstep``.

    These are handed to ``cases.generate_run_directory`` as
    ``extra_overrides``; that call raises ``CaseManifestError`` if the case
    does not declare the relevant cadence keys in
    ``allowed_input_overrides``, rather than silently skipping them.
    """
    if profile == DYNAMICS_ONLY:
        never_within_run = nstep + 1
        return {key: never_within_run for key in _CADENCE_KEYS}
    if profile == PRODUCTION_LIGHT:
        return dict(PRODUCTION_LIGHT_CADENCE)
    raise MeasurementProfileError(
        f"unknown measurement_profile {profile!r}; expected one of {PROFILES}"
    )
