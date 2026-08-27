"""Tests for the WP-09 campaign manifest schema, loader and cell resolution.

Nothing here executes the UppASD binary -- these are exactly the "schema
tests; manifest tests; runner parsing" categories blueprint section 16/WP-09
section F allows on ordinary shared CI, with no fixed wall-time threshold
anywhere in this file.
"""

from __future__ import annotations

import pathlib

import jsonschema
import pytest
import yaml

from harness import campaigns
from harness import cases as cases_mod
from harness import omp_sweep
from harness import omp_topology

BENCHMARKS_DIR = pathlib.Path(__file__).resolve().parent.parent
CAMPAIGNS_DIR = BENCHMARKS_DIR / "campaigns"
CASES_DIR = BENCHMARKS_DIR / "cases"

TRACKED_MANIFESTS = ["lean.yaml", "full_cpu.yaml", "full_cuda.yaml", "full_crossover.yaml", "ci_infra.yaml"]


def _base_manifest(**overrides):
    manifest = {
        "schema_version": "1.0.0",
        "id": "test_campaign",
        "description": "A synthetic campaign for schema testing.",
        "tier": "LEAN",
        "backend_scope": "CPU_ONLY",
        "cases": [
            {"case_id": "B01_bccFe", "variant_ids": ["bcc_fe_t0"], "size_ids": ["13x13x13"]},
        ],
        "cpu": {"thread_policy": "LEAN_SAMPLE", "proc_bind": "close"},
        "sample_count": 2,
        "measurement_profile": "DYNAMICS_ONLY",
        "environment_quality_mode": "DEVELOPER",
        "report_banner": campaigns.LEAN_REPORT_BANNER,
    }
    manifest.update(overrides)
    return manifest


def _write_manifest(tmp_path, manifest, name="campaign.yaml"):
    path = tmp_path / name
    path.write_text(yaml.safe_dump(manifest, sort_keys=False, allow_unicode=True))
    return path


# ---------------------------------------------------------------------------
# The campaign schema itself
# ---------------------------------------------------------------------------

def test_campaign_schema_is_a_valid_json_schema():
    jsonschema.Draft202012Validator.check_schema(campaigns._load_campaign_schema())


def test_campaign_schema_rejects_unknown_fields():
    assert campaigns._load_campaign_schema()["additionalProperties"] is False


# ---------------------------------------------------------------------------
# Manifest loading and validation
# ---------------------------------------------------------------------------

def test_valid_lean_manifest_loads(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest())
    campaign = campaigns.load_campaign_manifest(path)
    assert campaign.id == "test_campaign"
    assert campaign.tier == campaigns.LEAN
    assert campaign.authoritative is False
    assert campaign.has_cpu is True
    assert campaign.has_gpu is False
    assert campaign.report_banner == campaigns.LEAN_REPORT_BANNER
    assert campaign.require_clean_environment is False


def test_full_manifest_is_authoritative_and_strict(tmp_path):
    manifest = _base_manifest(
        tier="FULL", sample_count=5, environment_quality_mode="STRICT", report_banner=None,
    )
    del manifest["report_banner"]
    path = _write_manifest(tmp_path, manifest)
    campaign = campaigns.load_campaign_manifest(path)
    assert campaign.authoritative is True
    assert campaign.require_clean_environment is True


def test_lean_tier_requires_report_banner(tmp_path):
    manifest = _base_manifest()
    del manifest["report_banner"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_lean_tier_rejects_wrong_banner_text(tmp_path):
    manifest = _base_manifest(report_banner="not the real banner")
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_full_tier_forbids_report_banner(tmp_path):
    manifest = _base_manifest(
        tier="FULL", sample_count=5, environment_quality_mode="STRICT",
        report_banner=campaigns.LEAN_REPORT_BANNER,
    )
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_full_tier_requires_at_least_five_samples(tmp_path):
    manifest = _base_manifest(
        tier="FULL", sample_count=3, environment_quality_mode="STRICT", report_banner=None,
    )
    del manifest["report_banner"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_lean_tier_forbids_strict_environment(tmp_path):
    manifest = _base_manifest(environment_quality_mode="STRICT")
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_full_tier_forbids_developer_environment(tmp_path):
    manifest = _base_manifest(
        tier="FULL", sample_count=5, environment_quality_mode="DEVELOPER", report_banner=None,
    )
    del manifest["report_banner"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_cpu_only_requires_cpu_block(tmp_path):
    manifest = _base_manifest()
    del manifest["cpu"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_cpu_only_forbids_gpu_block(tmp_path):
    manifest = _base_manifest(gpu={"backend": "CUDA", "precisions": ["SINGLE"]})
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_cuda_only_requires_gpu_block(tmp_path):
    manifest = _base_manifest(backend_scope="CUDA_ONLY")
    del manifest["cpu"]
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_cuda_only_manifest_loads(tmp_path):
    manifest = _base_manifest(backend_scope="CUDA_ONLY", gpu={"backend": "CUDA", "precisions": ["SINGLE", "DOUBLE"]})
    del manifest["cpu"]
    path = _write_manifest(tmp_path, manifest)
    campaign = campaigns.load_campaign_manifest(path)
    assert campaign.has_cpu is False
    assert campaign.has_gpu is True
    assert campaign.gpu_precisions == ["SINGLE", "DOUBLE"]
    assert campaign.gpu_index == 0


def test_cpu_and_cuda_requires_single_campaign_execution_flag(tmp_path):
    manifest = _base_manifest(
        backend_scope="CPU_AND_CUDA", gpu={"backend": "CUDA", "precisions": ["SINGLE"]},
    )
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.load_campaign_manifest(path)


def test_cpu_and_cuda_with_flag_loads(tmp_path):
    manifest = _base_manifest(
        backend_scope="CPU_AND_CUDA", gpu={"backend": "CUDA", "precisions": ["SINGLE"]},
        single_campaign_execution=True,
    )
    path = _write_manifest(tmp_path, manifest)
    campaign = campaigns.load_campaign_manifest(path)
    assert campaign.has_cpu is True
    assert campaign.has_gpu is True


# ---------------------------------------------------------------------------
# The five tracked manifests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("filename", TRACKED_MANIFESTS)
def test_tracked_manifest_loads(filename):
    campaigns.load_campaign_manifest(CAMPAIGNS_DIR / filename)


def test_lean_yaml_is_lean_and_covers_two_cases():
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "lean.yaml")
    assert campaign.tier == campaigns.LEAN
    assert campaign.authoritative is False
    assert {c["case_id"] for c in campaign.cases} == {"B01_bccFe", "B02_skyrmion2D"}
    assert campaign.backend_scope == campaigns.CPU_AND_CUDA


def test_full_cpu_yaml_is_authoritative_cpu_only():
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "full_cpu.yaml")
    assert campaign.authoritative is True
    assert campaign.backend_scope == campaigns.CPU_ONLY
    assert campaign.has_gpu is False


def test_full_cuda_yaml_is_authoritative_cuda_only():
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "full_cuda.yaml")
    assert campaign.authoritative is True
    assert campaign.backend_scope == campaigns.CUDA_ONLY
    assert set(campaign.gpu_precisions) == {"SINGLE", "DOUBLE"}


def test_full_crossover_yaml_combines_both_backends_in_one_campaign():
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "full_crossover.yaml")
    assert campaign.authoritative is True
    assert campaign.backend_scope == campaigns.CPU_AND_CUDA
    assert campaign.single_campaign_execution is True


def test_all_five_families_covered_across_full_manifests():
    families = {"B01_bccFe", "B02_skyrmion2D", "B03_skyrmion3D", "B04_dhcpNd", "B05_dipoleFFT"}
    for filename in ("full_cpu.yaml", "full_cuda.yaml", "full_crossover.yaml"):
        campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / filename)
        assert {c["case_id"] for c in campaign.cases} == families


def test_ci_infra_yaml_is_lean_and_infrastructure_only():
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "ci_infra.yaml")
    assert campaign.tier == campaigns.LEAN
    assert [c["case_id"] for c in campaign.cases] == ["INFRA_test_only"]


# ---------------------------------------------------------------------------
# resolve_cells
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def real_cases():
    return cases_mod.discover_cases(CASES_DIR)


def test_resolve_cells_expands_all_wildcard(real_cases):
    manifest = _base_manifest(
        cases=[{"case_id": "B01_bccFe", "variant_ids": "ALL", "size_ids": "ALL"}],
    )
    campaign = campaigns.Campaign(manifest)
    cells = campaigns.resolve_cells(campaign, real_cases)
    case = real_cases["B01_bccFe"]
    expected_variants = {v["id"] for v in case.manifest["variants"]}
    expected_sizes = {s["id"] for s in case.manifest["sizes"]}
    assert {c["variant_id"] for c in cells} == expected_variants
    assert {c["size_id"] for c in cells} == expected_sizes
    assert all(c["backend"] == "CPU" for c in cells)
    assert len(cells) == len(expected_variants) * len(expected_sizes)


def test_resolve_cells_gpu_cells_carry_precision(real_cases):
    manifest = _base_manifest(
        backend_scope="CPU_AND_CUDA",
        cases=[{"case_id": "B01_bccFe", "variant_ids": ["bcc_fe_t0"], "size_ids": ["13x13x13"]}],
        gpu={"backend": "CUDA", "precisions": ["SINGLE", "DOUBLE"]},
        single_campaign_execution=True,
    )
    campaign = campaigns.Campaign(manifest)
    cells = campaigns.resolve_cells(campaign, real_cases)
    backends = sorted((c["backend"], c["requested_precision"]) for c in cells)
    assert backends == [("CPU", "DOUBLE"), ("GPU", "DOUBLE"), ("GPU", "SINGLE")]


def test_resolve_cells_rejects_unknown_case(real_cases):
    manifest = _base_manifest(cases=[{"case_id": "NOT_A_CASE", "variant_ids": "ALL", "size_ids": "ALL"}])
    campaign = campaigns.Campaign(manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.resolve_cells(campaign, real_cases)


def test_resolve_cells_rejects_unknown_variant(real_cases):
    manifest = _base_manifest(
        cases=[{"case_id": "B01_bccFe", "variant_ids": ["not_a_variant"], "size_ids": "ALL"}]
    )
    campaign = campaigns.Campaign(manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.resolve_cells(campaign, real_cases)


def test_resolve_cells_rejects_infrastructure_only_case_for_authoritative_campaign(real_cases):
    manifest = _base_manifest(
        tier="FULL", sample_count=5, environment_quality_mode="STRICT", report_banner=None,
        cases=[{"case_id": "INFRA_test_only", "variant_ids": "ALL", "size_ids": "ALL"}],
    )
    del manifest["report_banner"]
    campaign = campaigns.Campaign(manifest)
    with pytest.raises(cases_mod.CaseManifestError):
        campaigns.resolve_cells(campaign, real_cases)


def test_resolve_cells_allows_infrastructure_only_case_for_lean_campaign(real_cases):
    campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / "ci_infra.yaml")
    cells = campaigns.resolve_cells(campaign, real_cases)
    assert len(cells) > 0
    assert all(c["case_id"] == "INFRA_test_only" for c in cells)


def test_full_manifests_never_resolve_infrastructure_only_case(real_cases):
    for filename in ("full_cpu.yaml", "full_cuda.yaml", "full_crossover.yaml"):
        campaign = campaigns.load_campaign_manifest(CAMPAIGNS_DIR / filename)
        # Must not raise: none of these manifests name INFRA_test_only.
        campaigns.resolve_cells(campaign, real_cases)


# ---------------------------------------------------------------------------
# Machine hardware policy: thread-policy resolution against real topology
# ---------------------------------------------------------------------------

def _synthetic_topology(physical_cores):
    return omp_topology.HostTopology(
        cpu_sockets=1,
        physical_cores=physical_cores,
        logical_threads=physical_cores,
        physical_cpu_ids=list(range(physical_cores)),
        numa_node_cpus={0: set(range(physical_cores))},
        metadata_incomplete=False,
    )


def test_lean_sweep_manifest_includes_one_moderate_and_full():
    topology = _synthetic_topology(16)
    manifest = campaigns.build_lean_sweep_manifest(topology)
    assert manifest["thread_counts"] == [1, 8, 16]
    assert manifest["smt_mode"] == "physical_only"


def test_lean_sweep_manifest_deduplicates_on_tiny_hosts():
    topology = _synthetic_topology(1)
    manifest = campaigns.build_lean_sweep_manifest(topology)
    assert manifest["thread_counts"] == [1]


def test_full_physical_sweep_manifest_covers_every_physical_count():
    topology = _synthetic_topology(4)
    manifest = campaigns.build_full_physical_sweep_manifest(topology)
    assert manifest["thread_counts"] == [1, 2, 3, 4]


def test_resolve_cpu_sweep_returns_hardware_validated_omp_sweep():
    topology = _synthetic_topology(4)
    manifest = _base_manifest(cpu={"thread_policy": "LEAN_SAMPLE", "proc_bind": "close"})
    campaign = campaigns.Campaign(manifest)
    sweep = campaigns.resolve_cpu_sweep(campaign, topology)
    configurations = omp_sweep.resolve_run_configurations(sweep, topology)
    assert [c.threads for c in configurations] == [1, 2, 4]


def test_resolve_cpu_sweep_rejects_gpu_only_campaign():
    manifest = _base_manifest(backend_scope="CUDA_ONLY", gpu={"backend": "CUDA", "precisions": ["SINGLE"]})
    del manifest["cpu"]
    campaign = campaigns.Campaign(manifest)
    with pytest.raises(campaigns.CampaignManifestError):
        campaigns.resolve_cpu_sweep(campaign, _synthetic_topology(4))
