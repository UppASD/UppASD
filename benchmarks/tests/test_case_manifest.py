"""Tests for the WP-02 case-manifest infrastructure.

These are infrastructure fixtures only, exercising manifest parsing,
immutable-template handling, allow-listed input overrides, dimensional
scaling, template fingerprints and the workload-metadata hooks. No UppASD
executable is run.
"""

import json
import pathlib
import shutil

import pytest
import yaml

from harness import cases
from harness import workload_metadata

BENCHMARKS_DIR = pathlib.Path(__file__).resolve().parent.parent
CASES_DIR = BENCHMARKS_DIR / "cases"
INFRA_CASE_PATH = CASES_DIR / "INFRA_test_only" / "case.yaml"


def _template_dir():
    return CASES_DIR / "INFRA_test_only" / "template"


def _base_manifest():
    return yaml.safe_load(INFRA_CASE_PATH.read_text())


def _write_manifest(tmp_path, manifest, name="case.yaml"):
    dst_template = tmp_path / "template"
    if not dst_template.exists():
        shutil.copytree(_template_dir(), dst_template)
    manifest = dict(manifest)
    manifest["template_directory"] = "template"
    path = tmp_path / name
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return path


# ---------------------------------------------------------------------------
# The case-manifest schema itself
# ---------------------------------------------------------------------------

def test_case_manifest_schema_is_a_valid_json_schema():
    import jsonschema

    jsonschema.Draft202012Validator.check_schema(cases._load_case_schema())


def test_case_manifest_schema_rejects_unknown_fields():
    assert cases._load_case_schema()["additionalProperties"] is False


# ---------------------------------------------------------------------------
# Manifest loading and validation
# ---------------------------------------------------------------------------

def test_infra_case_is_a_valid_manifest():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    assert case.id == "INFRA_test_only"
    assert case.infrastructure_test_only is True
    assert case.template_dir.is_dir()


def test_invalid_manifest_missing_required_field(tmp_path):
    manifest = _base_manifest()
    del manifest["workload_class"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


def test_invalid_manifest_bad_enum(tmp_path):
    manifest = _base_manifest()
    manifest["dimensionality"] = "4D"
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


def test_invalid_manifest_missing_template_directory(tmp_path):
    manifest = _base_manifest()
    manifest["template_directory"] = "does_not_exist"
    path = tmp_path / "case.yaml"
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


# ---------------------------------------------------------------------------
# Allow-listed input overrides
# ---------------------------------------------------------------------------

def test_case_level_override_outside_global_safe_set_is_rejected(tmp_path):
    manifest = _base_manifest()
    manifest["allowed_input_overrides"].append("exchange")
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


def test_variant_override_outside_case_allow_list_is_rejected(tmp_path):
    manifest = _base_manifest()
    manifest["variants"][0]["overrides"]["exchange"] = "./other_jfile"
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


def test_runtime_override_outside_allow_list_is_rejected():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    with pytest.raises(cases.CaseManifestError):
        cases.build_overrides(case, "T0", "2x2x2", extra_overrides={"exchange": "./other_jfile"})


def test_runtime_override_in_allow_list_is_applied():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    overrides = cases.build_overrides(case, "T0", "2x2x2", extra_overrides={"Nstep": 50})
    assert overrides["Nstep"] == 50
    assert overrides["temp"] == 0.0
    assert overrides["ncell"] == [2, 2, 2]


def test_manifest_without_ncell_in_allow_list_is_rejected(tmp_path):
    manifest = _base_manifest()
    manifest["allowed_input_overrides"] = [
        k for k in manifest["allowed_input_overrides"] if k != "ncell"
    ]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


# ---------------------------------------------------------------------------
# Dimensional scaling
# ---------------------------------------------------------------------------

def test_legal_3d_replication_folds_into_ncell():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    size = case.resolve_size("4x4x4")
    assert cases.replication_to_ncell(case, size) == [4, 4, 4]


def test_legal_2d_replication_uses_fixed_thickness(tmp_path):
    manifest = _base_manifest()
    manifest["dimensionality"] = "2D"
    manifest["scalable_dimensions"] = ["Nx", "Ny"]
    manifest["thickness_cells"] = 1
    manifest["sizes"] = [{"id": "8x8", "replication": [8, 8]}]
    path = _write_manifest(tmp_path, manifest)
    case = cases.load_case_manifest(path)
    size = case.resolve_size("8x8")
    assert cases.replication_to_ncell(case, size) == [8, 8, 1]


def test_2d_size_with_wrong_tuple_length_is_rejected(tmp_path):
    manifest = _base_manifest()
    manifest["dimensionality"] = "2D"
    manifest["scalable_dimensions"] = ["Nx", "Ny"]
    manifest["thickness_cells"] = 1
    manifest["sizes"] = [{"id": "bad", "replication": [8, 8, 8]}]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(cases.CaseManifestError):
        cases.load_case_manifest(path)


def test_case_specific_replication_is_used_verbatim():
    manifest = _base_manifest()
    manifest["dimensionality"] = "case_specific"
    manifest["scalable_dimensions"] = ["hex_ring"]
    manifest["sizes"] = [{"id": "ring3", "replication": [3, 3, 3]}]
    case = cases.Case(manifest=manifest, manifest_path=INFRA_CASE_PATH, template_dir=_template_dir())
    size = case.resolve_size("ring3")
    assert cases.replication_to_ncell(case, size) == [3, 3, 3]


# ---------------------------------------------------------------------------
# Immutable templates, generated run directories, fingerprints
# ---------------------------------------------------------------------------

def test_source_template_is_never_modified(tmp_path):
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    before = cases.compute_template_fingerprint(case.template_dir)
    before_inpsd = (case.template_dir / "inpsd.dat").read_text()

    cases.generate_run_directory(case, "finiteT", "4x4x4", tmp_path)

    after = cases.compute_template_fingerprint(case.template_dir)
    after_inpsd = (case.template_dir / "inpsd.dat").read_text()
    assert before == after
    assert before_inpsd == after_inpsd


def test_generated_run_directory_is_separate_from_template(tmp_path):
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    run = cases.generate_run_directory(case, "T0", "2x2x2", tmp_path)
    assert run.path.resolve() != case.template_dir.resolve()
    inpsd_text = (run.path / "inpsd.dat").read_text()
    assert "ncell 2 2 2" in inpsd_text
    assert "temp 0.0" in inpsd_text


def test_generate_run_directory_refuses_the_template_directory_itself():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    with pytest.raises(cases.CaseManifestError):
        cases.generate_run_directory(case, "T0", "2x2x2", case.template_dir.parent, run_id="template")


def test_generate_run_directory_refuses_to_reuse_an_existing_directory(tmp_path):
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    cases.generate_run_directory(case, "T0", "2x2x2", tmp_path, run_id="dup")
    with pytest.raises(cases.CaseManifestError):
        cases.generate_run_directory(case, "T0", "2x2x2", tmp_path, run_id="dup")


def test_template_fingerprint_is_deterministic():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    first = cases.compute_template_fingerprint(case.template_dir)
    second = cases.compute_template_fingerprint(case.template_dir)
    assert first == second
    assert len(first) == 64  # sha256 hex digest


def test_template_fingerprint_changes_with_content(tmp_path):
    src = tmp_path / "tpl"
    shutil.copytree(_template_dir(), src)
    original = cases.compute_template_fingerprint(src)
    (src / "inpsd.dat").write_text((src / "inpsd.dat").read_text() + "\nNstep 99\n")
    changed = cases.compute_template_fingerprint(src)
    assert original != changed


def test_generated_run_directory_records_the_fingerprint_used(tmp_path):
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    run = cases.generate_run_directory(case, "T0", "2x2x2", tmp_path)
    recorded = json.loads((run.path / "_benchmark_run.json").read_text())
    assert recorded["template_fingerprint"] == run.template_fingerprint
    assert recorded["template_fingerprint"] == cases.compute_template_fingerprint(case.template_dir)


# ---------------------------------------------------------------------------
# Infrastructure-only cases cannot enter authoritative campaigns
# ---------------------------------------------------------------------------

def test_infra_case_is_excluded_from_authoritative_campaigns():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    filtered = cases.filter_cases_for_campaign({case.id: case}, authoritative=True)
    assert case.id not in filtered


def test_infra_case_is_kept_for_non_authoritative_campaigns():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    filtered = cases.filter_cases_for_campaign({case.id: case}, authoritative=False)
    assert case.id in filtered


def test_require_not_infrastructure_only_raises_for_infra_case():
    case = cases.load_case_manifest(INFRA_CASE_PATH)
    with pytest.raises(cases.CaseManifestError):
        cases.require_not_infrastructure_only(case)


# ---------------------------------------------------------------------------
# discover_cases
# ---------------------------------------------------------------------------

def test_discover_cases_finds_the_infra_fixture():
    found = cases.discover_cases(CASES_DIR)
    assert "INFRA_test_only" in found


def test_discover_cases_rejects_duplicate_ids(tmp_path):
    manifest = _base_manifest()
    (tmp_path / "case_a").mkdir()
    (tmp_path / "case_b").mkdir()
    shutil.copytree(_template_dir(), tmp_path / "case_a" / "template")
    shutil.copytree(_template_dir(), tmp_path / "case_b" / "template")
    for sub in ("case_a", "case_b"):
        m = dict(manifest)
        m["template_directory"] = "template"
        (tmp_path / sub / "case.yaml").write_text(yaml.safe_dump(m, sort_keys=False))
    with pytest.raises(cases.CaseManifestError):
        cases.discover_cases(tmp_path)


# ---------------------------------------------------------------------------
# Workload metadata hooks (section E)
# ---------------------------------------------------------------------------

STRUCT_FIXTURE = """\
#######################################################
# Number of atoms:       4
# Maximum num of neighbours:        2
#######################################################
#  iatom jatom  itype  jtype        r_{ij}^x        r_{ij}^y        r_{ij}^z          J_{ij}         |r_{ij}|
       1        2      1      1    5.000000E-01    0.000000E+00    0.000000E+00    1.000000E+00    5.000000E-01
       1        3      1      1    0.000000E+00    5.000000E-01    0.000000E+00    1.000000E+00    5.000000E-01
       2        1      1      1   -5.000000E-01    0.000000E+00    0.000000E+00    1.000000E+00    5.000000E-01
       2        4      1      1    0.000000E+00    5.000000E-01    0.000000E+00    1.000000E+00    5.000000E-01
       3        1      1      1    0.000000E+00   -5.000000E-01    0.000000E+00    1.000000E+00    5.000000E-01
       3        4      1      1    5.000000E-01    0.000000E+00    0.000000E+00    1.000000E+00    5.000000E-01
       4        2      1      1    0.000000E+00   -5.000000E-01    0.000000E+00    1.000000E+00    5.000000E-01
       4        3      1      1   -5.000000E-01    0.000000E+00    0.000000E+00    1.000000E+00    5.000000E-01
"""


def test_neighbor_list_parser_reads_real_struct_output(tmp_path):
    (tmp_path / "inpsd.dat").write_text("simid demo\nncell 1 1 1\n")
    (tmp_path / "struct.demo.out").write_text(STRUCT_FIXTURE)
    metadata = workload_metadata.neighbor_list_from_struct_output(None, None, tmp_path)
    assert metadata["natom"] == 4
    assert metadata["max_neighbors"] == 2
    assert metadata["directed_interactions"] == 8
    assert metadata["mean_neighbors"] == pytest.approx(2.0)
    assert metadata["mean_neighbors"] <= metadata["max_neighbors"]


def test_neighbor_list_parser_requires_struct_output(tmp_path):
    (tmp_path / "inpsd.dat").write_text("simid demo\n")
    with pytest.raises(workload_metadata.WorkloadMetadataError):
        workload_metadata.neighbor_list_from_struct_output(None, None, tmp_path)


def test_fft_grid_parser_derives_padded_grid_from_replication(tmp_path):
    manifest = {
        "id": "FFT_synthetic",
        "dimensionality": "3D",
        "scalable_dimensions": ["Nx", "Ny", "Nz"],
    }
    case = cases.Case(manifest=manifest, manifest_path=tmp_path, template_dir=tmp_path)
    size = {"id": "3x3x2", "replication": [3, 3, 2]}
    (tmp_path / "inpsd.dat").write_text("simid demo\nposfile ./posfile\n")
    (tmp_path / "posfile").write_text("1 1 0 0 0\n2 1 0.5 0.5 0.5\n")

    metadata = workload_metadata.fft_grid_from_replication(case, size, tmp_path)

    assert metadata["fft_grid"] == [3, 3, 2]
    assert metadata["fft_grid_padded"] == [5, 5, 3]
    assert metadata["fft_grid_points"] == 5 * 5 * 3
    assert metadata["natom"] == 2 * 3 * 3 * 2


def test_compute_workload_metadata_dispatches_by_method_name(tmp_path):
    (tmp_path / "inpsd.dat").write_text("simid demo\nncell 1 1 1\n")
    (tmp_path / "struct.demo.out").write_text(STRUCT_FIXTURE)
    manifest = {"id": "demo", "workload_metadata_method": "neighbor_list_from_struct_output"}
    case = cases.Case(manifest=manifest, manifest_path=tmp_path, template_dir=tmp_path)
    metadata = workload_metadata.compute_workload_metadata(case, None, tmp_path)
    assert metadata["natom"] == 4


def test_compute_workload_metadata_rejects_unknown_method(tmp_path):
    manifest = {"id": "demo", "workload_metadata_method": "made_up"}
    case = cases.Case(manifest=manifest, manifest_path=tmp_path, template_dir=tmp_path)
    with pytest.raises(workload_metadata.WorkloadMetadataError):
        workload_metadata.compute_workload_metadata(case, None, tmp_path)
