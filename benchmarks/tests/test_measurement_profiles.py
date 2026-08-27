"""Unit tests for the WP-03 measurement-profile overrides (section E)."""

import pytest

from harness import cases
from harness import measurement_profiles as mp


@pytest.fixture
def infra_case():
    import pathlib

    path = pathlib.Path(__file__).resolve().parent.parent / "cases" / "INFRA_test_only" / "case.yaml"
    return cases.load_case_manifest(path)


def test_dynamics_only_pushes_every_cadence_past_nstep():
    overrides = mp.profile_overrides(mp.DYNAMICS_ONLY, nstep=500)
    assert overrides == {"avrg_step": 501, "cumu_step": 501, "tottraj_step": 501, "ene_step": 501}
    for value in overrides.values():
        assert value > 500


def test_production_light_uses_the_documented_fixed_cadence():
    overrides = mp.profile_overrides(mp.PRODUCTION_LIGHT, nstep=1_000_000)
    assert overrides == mp.PRODUCTION_LIGHT_CADENCE
    # Fixed regardless of nstep -- this is the point of a *representative*
    # cadence, not one scaled to make the run always sample once.
    assert mp.profile_overrides(mp.PRODUCTION_LIGHT, nstep=10) == overrides


def test_unknown_profile_is_rejected():
    with pytest.raises(mp.MeasurementProfileError):
        mp.profile_overrides("SOMETHING_ELSE", nstep=10)


def test_dynamics_only_overrides_are_legal_for_the_infra_case(infra_case):
    overrides = mp.profile_overrides(mp.DYNAMICS_ONLY, nstep=10)
    combined = cases.build_overrides(infra_case, "T0", "2x2x2", extra_overrides=overrides)
    for key, value in overrides.items():
        assert combined[key] == value


def test_production_light_overrides_are_legal_for_the_infra_case(infra_case):
    overrides = mp.profile_overrides(mp.PRODUCTION_LIGHT, nstep=10)
    combined = cases.build_overrides(infra_case, "T0", "2x2x2", extra_overrides=overrides)
    for key, value in overrides.items():
        assert combined[key] == value
