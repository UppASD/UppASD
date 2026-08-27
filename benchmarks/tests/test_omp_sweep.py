"""Tests for the WP-05 declarative OpenMP sweep manifest and run-configuration
binding.
"""

from __future__ import annotations

import jsonschema
import pytest
import yaml

from harness import omp_sweep
from harness import omp_topology


def _base_manifest(**overrides):
    manifest = {
        "schema_version": "1.0.0",
        "id": "standard_physical",
        "description": "Standard physical-cores-first thread sweep.",
        "thread_counts": [1, 2, 4],
        "proc_bind": "close",
        "places": "cores",
        "dynamic": False,
        "smt_mode": "physical_only",
    }
    manifest.update(overrides)
    return manifest


def _write_manifest(tmp_path, manifest, name="sweep.yaml"):
    path = tmp_path / name
    path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    return path


# ---------------------------------------------------------------------------
# The sweep schema itself
# ---------------------------------------------------------------------------

def test_omp_sweep_schema_is_a_valid_json_schema():
    jsonschema.Draft202012Validator.check_schema(omp_sweep._load_sweep_schema())


def test_omp_sweep_schema_rejects_unknown_fields():
    assert omp_sweep._load_sweep_schema()["additionalProperties"] is False


# ---------------------------------------------------------------------------
# Manifest loading and validation
# ---------------------------------------------------------------------------

def test_valid_manifest_loads(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest())
    sweep = omp_sweep.load_omp_sweep_manifest(path)
    assert sweep.id == "standard_physical"
    assert sweep.thread_counts == [1, 2, 4]
    assert sweep.proc_bind == "close"
    assert sweep.places == "cores"
    assert sweep.smt_mode == "physical_only"
    assert sweep.allow_smt is False


def test_smt_extension_manifest_allow_smt_is_true(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest(id="smt_extension", smt_mode="smt_extension"))
    sweep = omp_sweep.load_omp_sweep_manifest(path)
    assert sweep.allow_smt is True


def test_manifest_must_include_thread_count_one(tmp_path):
    manifest = _base_manifest(thread_counts=[2, 4, 8])
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


def test_manifest_rejects_duplicate_thread_counts(tmp_path):
    manifest = _base_manifest(thread_counts=[1, 2, 2])
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


def test_manifest_rejects_sweeping_both_proc_binds_at_once(tmp_path):
    # proc_bind is a single enum value, not an array -- one campaign-level choice.
    manifest = _base_manifest(proc_bind=["close", "spread"])
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


def test_manifest_rejects_places_other_than_cores(tmp_path):
    manifest = _base_manifest(places="threads")
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


def test_manifest_rejects_dynamic_true(tmp_path):
    manifest = _base_manifest(dynamic=True)
    path = _write_manifest(tmp_path, manifest)
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


def test_manifest_rejects_unsupported_extension(tmp_path):
    path = tmp_path / "sweep.txt"
    path.write_text("not a manifest")
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.load_omp_sweep_manifest(path)


# ---------------------------------------------------------------------------
# build_omp_env
# ---------------------------------------------------------------------------

def test_build_omp_env_sets_every_binding_variable_explicitly():
    env = omp_sweep.build_omp_env(4, proc_bind="close", base_env={"PATH": "/usr/bin"})
    assert env["OMP_NUM_THREADS"] == "4"
    assert env["OMP_PLACES"] == "cores"
    assert env["OMP_PROC_BIND"] == "close"
    assert env["OMP_DYNAMIC"] == "FALSE"
    assert env["PATH"] == "/usr/bin"  # base env passed through untouched


def test_build_omp_env_supports_spread():
    env = omp_sweep.build_omp_env(2, proc_bind="spread")
    assert env["OMP_PROC_BIND"] == "spread"


def test_build_omp_env_rejects_unknown_proc_bind():
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.build_omp_env(2, proc_bind="master")


def test_build_omp_env_rejects_non_cores_places():
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.build_omp_env(2, proc_bind="close", places="sockets")


def test_build_omp_env_rejects_non_positive_threads():
    with pytest.raises(omp_sweep.OmpSweepError):
        omp_sweep.build_omp_env(0, proc_bind="close")


# ---------------------------------------------------------------------------
# resolve_run_configurations -- end to end against a synthetic topology
# ---------------------------------------------------------------------------

def _synthetic_topology():
    return omp_topology.HostTopology(
        cpu_sockets=2,
        physical_cores=8,
        logical_threads=16,
        physical_cpu_ids=[0, 2, 4, 6, 8, 10, 12, 14],
        numa_node_cpus={0: set(range(0, 8)), 1: set(range(8, 16))},
        metadata_incomplete=False,
    )


def test_resolve_run_configurations_produces_one_config_per_thread_count(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest(thread_counts=[1, 2, 8]))
    sweep = omp_sweep.load_omp_sweep_manifest(path)
    topology = _synthetic_topology()

    configs = omp_sweep.resolve_run_configurations(sweep, topology, base_env={})
    assert [c.threads for c in configs] == [1, 2, 8]
    assert configs[0].numa_placement == omp_topology.NUMA_SINGLE_NODE
    assert configs[-1].numa_placement == omp_topology.NUMA_FULL_NODE
    for config in configs:
        assert config.env["OMP_NUM_THREADS"] == str(config.threads)
        assert config.env["OMP_PROC_BIND"] == "close"
        assert config.env["OMP_DYNAMIC"] == "FALSE"


def test_resolve_run_configurations_rejects_oversubscribing_physical_cores(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest(thread_counts=[1, 16]))
    sweep = omp_sweep.load_omp_sweep_manifest(path)
    topology = _synthetic_topology()

    with pytest.raises(omp_topology.OmpTopologyError):
        omp_sweep.resolve_run_configurations(sweep, topology)


def test_resolve_run_configurations_smt_extension_allows_logical_threads(tmp_path):
    path = _write_manifest(tmp_path, _base_manifest(thread_counts=[1, 16], smt_mode="smt_extension"))
    sweep = omp_sweep.load_omp_sweep_manifest(path)
    topology = _synthetic_topology()

    configs = omp_sweep.resolve_run_configurations(sweep, topology, base_env={})
    assert [c.threads for c in configs] == [1, 16]
