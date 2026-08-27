"""Tests for the WP-05 host topology / OpenMP thread-count validation module.

Real hardware is not assumed available or a fixed shape in CI, so every test
builds a synthetic `/proc/cpuinfo` + `/sys/devices/system/node` tree (the
same style as `test_provenance.py`'s `cpu_provenance` tests) rather than
depending on the machine this happens to run on.
"""

from __future__ import annotations

import pytest

from harness import omp_topology


def _write_cpuinfo(tmp_path, *, sockets, cores_per_socket, threads_per_core):
    """A synthetic multi-socket, optionally-SMT `/proc/cpuinfo`."""
    blocks = []
    processor = 0
    for socket in range(sockets):
        for core in range(cores_per_socket):
            for _thread in range(threads_per_core):
                blocks.append(
                    f"processor\t: {processor}\n"
                    f"model name\t: Fake CPU\n"
                    f"physical id\t: {socket}\n"
                    f"core id\t: {core}\n"
                )
                processor += 1
    path = tmp_path / "cpuinfo"
    path.write_text("\n".join(blocks))
    return path


def _write_numa_tree(tmp_path, node_cpu_lists):
    """`node_cpu_lists` maps node id -> cpulist string, e.g. {0: "0-3", 1: "4-7"}."""
    node_root = tmp_path / "node"
    node_root.mkdir()
    for node_id, cpulist in node_cpu_lists.items():
        node_dir = node_root / f"node{node_id}"
        node_dir.mkdir()
        (node_dir / "cpulist").write_text(cpulist + "\n")
    return node_root


# ---------------------------------------------------------------------------
# parse_cpu_range_list
# ---------------------------------------------------------------------------

def test_parse_cpu_range_list():
    assert omp_topology.parse_cpu_range_list("0-3,5,7-8") == {0, 1, 2, 3, 5, 7, 8}
    assert omp_topology.parse_cpu_range_list("") == set()
    assert omp_topology.parse_cpu_range_list("4") == {4}


# ---------------------------------------------------------------------------
# capture_host_topology
# ---------------------------------------------------------------------------

def test_capture_host_topology_single_socket_no_smt(tmp_path):
    cpuinfo = _write_cpuinfo(tmp_path, sockets=1, cores_per_socket=4, threads_per_core=1)
    node_root = _write_numa_tree(tmp_path, {0: "0-3"})

    topology = omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)
    assert topology.cpu_sockets == 1
    assert topology.physical_cores == 4
    assert topology.logical_threads == 4
    assert topology.physical_cpu_ids == [0, 1, 2, 3]
    assert topology.smt_enabled is False
    assert topology.numa_nodes == 1
    assert topology.metadata_incomplete is False


def test_capture_host_topology_two_socket_smt(tmp_path):
    # 2 sockets x 4 cores x 2 SMT threads = 16 logical cpus, 8 physical cores.
    cpuinfo = _write_cpuinfo(tmp_path, sockets=2, cores_per_socket=4, threads_per_core=2)
    node_root = _write_numa_tree(tmp_path, {0: "0-7", 1: "8-15"})

    topology = omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)
    assert topology.cpu_sockets == 2
    assert topology.physical_cores == 8
    assert topology.logical_threads == 16
    assert topology.smt_enabled is True
    assert topology.numa_nodes == 2
    # One representative (lowest-numbered sibling) logical id per physical core.
    assert topology.physical_cpu_ids == [0, 2, 4, 6, 8, 10, 12, 14]


def test_capture_host_topology_missing_paths_is_honestly_incomplete(tmp_path):
    topology = omp_topology.capture_host_topology(
        cpuinfo_path=tmp_path / "no-such-file", node_root=tmp_path / "no-such-dir"
    )
    assert topology.metadata_incomplete is True
    assert topology.numa_nodes == 1


# ---------------------------------------------------------------------------
# validate_thread_count / validate_thread_sweep
# ---------------------------------------------------------------------------

@pytest.fixture
def four_core_topology(tmp_path):
    cpuinfo = _write_cpuinfo(tmp_path, sockets=1, cores_per_socket=4, threads_per_core=2)
    node_root = _write_numa_tree(tmp_path, {0: "0-7"})
    return omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)


def test_validate_thread_count_within_physical_cores(four_core_topology):
    assert omp_topology.validate_thread_count(four_core_topology, 4) == 4


def test_validate_thread_count_rejects_oversubscribing_physical_cores(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError, match="oversubscribe"):
        omp_topology.validate_thread_count(four_core_topology, 8, allow_smt=False)


def test_validate_thread_count_allows_smt_extension_up_to_logical_threads(four_core_topology):
    assert omp_topology.validate_thread_count(four_core_topology, 8, allow_smt=True) == 8
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.validate_thread_count(four_core_topology, 9, allow_smt=True)


def test_validate_thread_count_rejects_zero_or_negative(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.validate_thread_count(four_core_topology, 0)


def test_validate_thread_sweep_sorts_and_dedupes_valid_input(four_core_topology):
    assert omp_topology.validate_thread_sweep([4, 1, 2], four_core_topology) == [1, 2, 4]


def test_validate_thread_sweep_rejects_duplicates(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError, match="unique"):
        omp_topology.validate_thread_sweep([1, 2, 2], four_core_topology)


def test_validate_thread_sweep_rejects_empty(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.validate_thread_sweep([], four_core_topology)


def test_validate_thread_sweep_does_not_oversubscribe_the_node(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.validate_thread_sweep([1, 2, 4, 8], four_core_topology, allow_smt=False)


def test_physical_core_thread_counts(four_core_topology):
    assert omp_topology.physical_core_thread_counts(four_core_topology) == [1, 2, 3, 4]
    assert omp_topology.physical_core_thread_counts(four_core_topology, max_threads=2) == [1, 2]


# ---------------------------------------------------------------------------
# expected_cpu_ids_for_threads
# ---------------------------------------------------------------------------

def test_expected_cpu_ids_close_packs_first_n_physical_cores(four_core_topology):
    assert omp_topology.expected_cpu_ids_for_threads(four_core_topology, 2, proc_bind="close") == [0, 2]


def test_expected_cpu_ids_spread_distributes_across_places(four_core_topology):
    ids = omp_topology.expected_cpu_ids_for_threads(four_core_topology, 2, proc_bind="spread")
    assert len(ids) == 2
    assert set(ids).issubset(set(four_core_topology.physical_cpu_ids))


def test_expected_cpu_ids_rejects_unknown_proc_bind(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.expected_cpu_ids_for_threads(four_core_topology, 1, proc_bind="master")


def test_expected_cpu_ids_rejects_more_threads_than_physical_cores(four_core_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.expected_cpu_ids_for_threads(four_core_topology, 5, proc_bind="close")


# ---------------------------------------------------------------------------
# NUMA classification
# ---------------------------------------------------------------------------

@pytest.fixture
def two_socket_topology(tmp_path):
    # 2 sockets x 4 cores, no SMT -- physical_cpu_ids = [0..7], node0={0-3}, node1={4-7}.
    cpuinfo = _write_cpuinfo(tmp_path, sockets=2, cores_per_socket=4, threads_per_core=1)
    node_root = _write_numa_tree(tmp_path, {0: "0-3", 1: "4-7"})
    return omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)


def test_classify_numa_placement_single_node(two_socket_topology):
    assert omp_topology.classify_numa_placement(two_socket_topology, [0, 1]) == omp_topology.NUMA_SINGLE_NODE


def test_classify_numa_placement_full_node(two_socket_topology):
    assert omp_topology.classify_numa_placement(two_socket_topology, [0, 4]) == omp_topology.NUMA_FULL_NODE


def test_classify_numa_placement_multi_node_on_larger_host(tmp_path):
    cpuinfo = _write_cpuinfo(tmp_path, sockets=3, cores_per_socket=4, threads_per_core=1)
    node_root = _write_numa_tree(tmp_path, {0: "0-3", 1: "4-7", 2: "8-11"})
    topology = omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)
    assert omp_topology.classify_numa_placement(topology, [0, 4]) == omp_topology.NUMA_MULTI_NODE
    assert omp_topology.classify_numa_placement(topology, [0, 4, 8]) == omp_topology.NUMA_FULL_NODE


def test_classify_numa_placement_single_socket_host_is_always_single_node(tmp_path):
    cpuinfo = _write_cpuinfo(tmp_path, sockets=1, cores_per_socket=4, threads_per_core=1)
    node_root = _write_numa_tree(tmp_path, {0: "0-3"})
    topology = omp_topology.capture_host_topology(cpuinfo_path=cpuinfo, node_root=node_root)
    assert omp_topology.classify_numa_placement(topology, [0, 1, 2, 3]) == omp_topology.NUMA_SINGLE_NODE


def test_classify_numa_placement_rejects_empty_cpu_set(two_socket_topology):
    with pytest.raises(omp_topology.OmpTopologyError):
        omp_topology.classify_numa_placement(two_socket_topology, [])


def test_classify_thread_count_numa_close_binding_fills_one_socket_first(two_socket_topology):
    # threads=2 close-packs cpu ids [0, 1] -> both within node0.
    assert (
        omp_topology.classify_thread_count_numa(two_socket_topology, 2, proc_bind="close")
        == omp_topology.NUMA_SINGLE_NODE
    )
    # threads=8 uses every physical core -> spans both sockets.
    assert (
        omp_topology.classify_thread_count_numa(two_socket_topology, 8, proc_bind="close")
        == omp_topology.NUMA_FULL_NODE
    )


def test_classify_observed_affinity_uses_real_affinity_string(two_socket_topology):
    assert (
        omp_topology.classify_observed_affinity(two_socket_topology, "0-3")
        == omp_topology.NUMA_SINGLE_NODE
    )
    assert (
        omp_topology.classify_observed_affinity(two_socket_topology, "0,4")
        == omp_topology.NUMA_FULL_NODE
    )


def test_classify_observed_affinity_none_when_no_affinity_available(two_socket_topology):
    assert omp_topology.classify_observed_affinity(two_socket_topology, None) is None
    assert omp_topology.classify_observed_affinity(two_socket_topology, "") is None
