"""Validation tests for the UppASD benchmark result contract (WP-01).

These are infrastructure fixtures only: no UppASD executable is run, no
physics is exercised and no timing is measured. The fixtures exist to prove
that the contract can represent the records later work packages must write,
and that it rejects records that would silently corrupt a campaign.
"""

import json
import pathlib
import re

import jsonschema
import pytest

from harness import schema_validate

FIXTURES = pathlib.Path(__file__).resolve().parent / "fixtures"
VALID_DIR = FIXTURES / "valid"
INVALID_DIR = FIXTURES / "invalid"
EXPECTED_FAILURES = INVALID_DIR / "expected_failures.json"

VALID_RECORDS = sorted(VALID_DIR.glob("*.json"))
INVALID_RECORDS = sorted(
    path for path in INVALID_DIR.glob("*.json") if path != EXPECTED_FAILURES
)


def _load(path):
    with open(path) as handle:
        return json.load(handle)


def _expected_failures():
    return _load(EXPECTED_FAILURES)


# ---------------------------------------------------------------------------
# The schemas themselves
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("record_kind", sorted(schema_validate.SCHEMA_FILES))
def test_schema_is_a_valid_json_schema(record_kind):
    schema = schema_validate.load_schema(record_kind)
    jsonschema.Draft202012Validator.check_schema(schema)


@pytest.mark.parametrize("record_kind", sorted(schema_validate.SCHEMA_FILES))
def test_schema_rejects_unknown_fields(record_kind):
    schema = schema_validate.load_schema(record_kind)
    assert schema["additionalProperties"] is False


def test_unknown_record_kind_is_rejected():
    with pytest.raises(schema_validate.RecordInvalid):
        schema_validate.load_schema("timing_summary")


# ---------------------------------------------------------------------------
# Positive fixtures
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", VALID_RECORDS, ids=lambda p: p.stem)
def test_valid_fixture_validates(path):
    assert schema_validate.validate_file(path) is True


def test_required_positive_fixtures_are_present():
    """The contract must demonstrably represent every workload class in scope."""
    stems = {path.stem for path in VALID_RECORDS}
    assert "cpu_bccfe_dynamics_only" in stems
    assert "cuda_bccfe_single_precision" in stems
    assert "cpu_dhcpnd_long_range" in stems
    assert "cpu_dipole_fft_grid_only" in stems
    assert "hip_unmeasured_backend" in stems


# ---------------------------------------------------------------------------
# Negative fixtures
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", INVALID_RECORDS, ids=lambda p: p.stem)
def test_invalid_fixture_is_rejected(path):
    expected = _expected_failures()
    assert path.name in expected, (
        "%s has no entry in expected_failures.json; every negative fixture "
        "must document why it must fail" % path.name
    )
    with pytest.raises(schema_validate.RecordInvalid) as excinfo:
        schema_validate.validate_file(path)
    assert expected[path.name]["message_contains"] in str(excinfo.value)


def test_expected_failure_manifest_matches_fixtures():
    listed = set(_expected_failures())
    present = {path.name for path in INVALID_RECORDS}
    assert listed == present


# ---------------------------------------------------------------------------
# Contract invariants that later work packages depend on
# ---------------------------------------------------------------------------

CONTRACT_FIELDS = [
    # identity
    "schema_version", "run_id", "campaign_id", "case_id", "variant_id",
    "size_id", "sample_index",
    # physics/workload
    "natom", "directed_interactions", "mean_neighbors", "max_neighbors",
    "fft_grid", "fft_grid_points",
    # execution
    "backend", "gpu_backend", "omp_threads", "omp_places", "omp_proc_bind",
    "omp_dynamic",
    # precision
    "requested_precision", "effective_cpu_precision", "effective_gpu_precision",
    "comparison_precision_class",
    # source/build
    "git_commit", "git_dirty", "compiler", "compiler_version", "compile_flags",
    "cmake_options", "binary_checksum", "build_id",
    # hardware
    "cpu_model", "cpu_sockets", "cpu_physical_cores", "cpu_threads",
    "numa_nodes", "gpu_model", "gpu_id", "gpu_driver", "gpu_runtime",
    # simulation
    "temperature", "timestep", "nstep", "measurement_profile",
    # timing
    "timing_method", "process_wall_seconds", "setup_seconds",
    "steady_step_seconds",
    # quality
    "numerical_valid", "environment_valid", "quality_flags",
]


@pytest.mark.parametrize("field", CONTRACT_FIELDS)
def test_contract_field_is_required(field):
    schema = schema_validate.load_schema("raw_sample")
    assert field in schema["properties"]
    assert field in schema["required"]


def test_mpi_is_not_a_benchmark_dimension():
    """MPI is explicitly out of scope; no MPI field may enter the contract."""
    mpi_field = re.compile(r"(^|_)mpi(_|$)")
    for record_kind in schema_validate.SCHEMA_FILES:
        schema = schema_validate.load_schema(record_kind)
        assert not [name for name in schema["properties"] if mpi_field.search(name)]


def test_hip_backend_is_representable():
    schema = schema_validate.load_schema("raw_sample")
    assert "HIP" in schema["properties"]["gpu_backend"]["enum"]


def test_mixed_precision_has_an_unsupported_state():
    schema = schema_validate.load_schema("raw_sample")
    assert "MIXED" in schema["properties"]["requested_precision"]["enum"]
    assert "unsupported" in schema["properties"]["precision_support_state"]["enum"]


def test_requested_and_effective_precision_are_separate_fields():
    schema = schema_validate.load_schema("raw_sample")
    properties = schema["properties"]
    assert properties["requested_precision"] is not properties["effective_cpu_precision"]
    assert properties["effective_cpu_precision"] != properties["requested_precision"]


def test_numerical_and_environment_validity_are_separate_fields():
    schema = schema_validate.load_schema("raw_sample")
    assert schema["properties"]["numerical_valid"]["type"] == "boolean"
    assert schema["properties"]["environment_valid"]["type"] == "boolean"


def test_raw_records_carry_no_statistics():
    """Individual raw records must remain individual."""
    schema = schema_validate.load_schema("raw_sample")
    for statistic in ("sample_count", "median", "mad", "minimum", "maximum"):
        assert statistic not in schema["properties"]


def test_aggregates_carry_statistics_and_their_provenance():
    schema = schema_validate.load_schema("aggregate")
    for statistic in ("sample_count", "median", "mad", "minimum", "maximum"):
        assert statistic in schema["required"]
    assert "run_ids" in schema["required"]


def test_numerically_valid_but_environmentally_unsuitable_is_representable():
    """A correct run on a busy GPU is valid physics but not authoritative timing."""
    record = _load(VALID_DIR / "cuda_bccfe_single_precision.json")
    record["environment_valid"] = False
    record["quality_flags"] = ["gpu_busy"]
    assert schema_validate.validate_record(record) is True
    assert record["numerical_valid"] is True
