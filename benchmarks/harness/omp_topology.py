"""Host CPU/NUMA topology and OpenMP thread-count validation (WP-05 sections A-D).

Everything here reads real evidence -- `/proc/cpuinfo` and
`/sys/devices/system/node/nodeN/cpulist` -- the same honest-about-what-it-
could-not-determine style as `cpu_provenance.py`, whose `capture_cpu_topology`
this module builds on rather than duplicates. What `cpu_provenance` does not
provide is per-physical-core identity (which logical cpu ids are real cores
versus SMT siblings) and per-NUMA-node cpu membership; both are needed here to:

* validate a requested OpenMP thread count against real hardware without
  oversubscribing the node (section A);
* prioritize physical cores in the standard campaign and keep an SMT sweep
  clearly distinguished (section B);
* record whether a thread count's placement spans one NUMA node/socket,
  several, or the full node (section D).

Every filesystem root is a parameter, exactly like `cpu_provenance`, so tests
run against a synthetic tree instead of the real machine.
"""

from __future__ import annotations

import pathlib
import re

from harness import cpu_provenance

NUMA_SINGLE_NODE = "SINGLE_NODE"
NUMA_MULTI_NODE = "MULTI_NODE"
NUMA_FULL_NODE = "FULL_NODE"

_CPUINFO_BLOCK_RE = re.compile(r"\n\s*\n")


class OmpTopologyError(ValueError):
    """A requested OpenMP thread configuration is not legal on this host."""


class HostTopology:
    """A resolved snapshot of this host's CPU/NUMA layout.

    ``physical_cpu_ids`` is one representative logical cpu id per real
    physical core (the lowest-numbered sibling of each distinct
    ``(physical id, core id)`` pair), sorted ascending -- this is what
    ``OMP_PLACES=cores`` enumerates. ``numa_node_cpus`` maps NUMA node id to
    the *complete* set of logical cpu ids (physical cores and their SMT
    siblings alike) that belong to it, since a NUMA node owns whole cores,
    SMT siblings included.
    """

    def __init__(self, *, cpu_sockets, physical_cores, logical_threads,
                 physical_cpu_ids, numa_node_cpus, metadata_incomplete):
        self.cpu_sockets = cpu_sockets
        self.physical_cores = physical_cores
        self.logical_threads = logical_threads
        self.physical_cpu_ids = list(physical_cpu_ids)
        self.numa_node_cpus = dict(numa_node_cpus)
        self.metadata_incomplete = metadata_incomplete

    @property
    def smt_enabled(self):
        """Whether this host exposes more logical threads than physical cores."""
        return self.logical_threads > self.physical_cores

    @property
    def numa_nodes(self):
        return len(self.numa_node_cpus) or 1

    @property
    def logical_cpu_ids(self):
        """Every logical cpu id this host reports (physical cores and their
        SMT siblings alike), ascending. The SMT-extension place pool."""
        all_ids = set()
        for members in self.numa_node_cpus.values():
            all_ids |= members
        return sorted(all_ids)

    def __repr__(self):
        return (
            f"HostTopology(cpu_sockets={self.cpu_sockets!r}, "
            f"physical_cores={self.physical_cores!r}, "
            f"logical_threads={self.logical_threads!r}, "
            f"numa_nodes={self.numa_nodes!r}, smt_enabled={self.smt_enabled!r})"
        )


def parse_cpu_range_list(text):
    """Parse a Linux-style cpu range list (``"0-3,5,7-8"``) into a set of ints.

    This is the sysfs-side inverse of `cpu_provenance._format_affinity_ranges`;
    both `nodeN/cpulist` and `/proc/self/status`'s `Cpus_allowed_list` use this
    format. Blank input parses to an empty set rather than raising.
    """
    text = text.strip()
    if not text:
        return set()
    ids = set()
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            start, _, end = chunk.partition("-")
            ids.update(range(int(start), int(end) + 1))
        else:
            ids.add(int(chunk))
    return ids


def _parse_physical_core_ids(cpuinfo_text):
    """One representative logical cpu id per distinct `(physical id, core id)`.

    On hosts without those keys (some virtualized/ARM environments), falls
    back to treating every logical processor as its own physical core -- the
    same honest conservative fallback `cpu_provenance.parse_cpuinfo` uses.
    """
    representative_by_core = {}
    for block in _CPUINFO_BLOCK_RE.split(cpuinfo_text.strip()):
        fields = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
        if "processor" not in fields:
            continue
        processor = int(fields["processor"])
        physical_id = fields.get("physical id", "0")
        core_id = fields.get("core id", fields["processor"])
        core_key = (physical_id, core_id)
        representative_by_core[core_key] = min(processor, representative_by_core.get(core_key, processor))
    return sorted(representative_by_core.values())


def capture_host_topology(*, cpuinfo_path="/proc/cpuinfo", node_root="/sys/devices/system/node"):
    """Return a :class:`HostTopology` built from real `/proc/cpuinfo` and
    `/sys/devices/system/node` evidence.

    Falls back to `os.cpu_count()`-derived, single-node/single-core-per-thread
    assumptions when either path is unreadable, with `metadata_incomplete`
    set so callers know the fallback was used rather than measured.
    """
    cpuinfo_path = pathlib.Path(cpuinfo_path)
    topology = cpu_provenance.capture_cpu_topology(cpuinfo_path=cpuinfo_path, node_root=node_root)

    if cpuinfo_path.is_file():
        physical_cpu_ids = _parse_physical_core_ids(cpuinfo_path.read_text())
    else:
        physical_cpu_ids = list(range(topology["cpu_physical_cores"]))

    node_root = pathlib.Path(node_root)
    numa_node_cpus = {}
    if node_root.is_dir():
        for entry in sorted(node_root.iterdir()):
            match = re.fullmatch(r"node(\d+)", entry.name)
            if not match:
                continue
            cpulist_path = entry / "cpulist"
            if cpulist_path.is_file():
                numa_node_cpus[int(match.group(1))] = parse_cpu_range_list(cpulist_path.read_text())
    if not numa_node_cpus:
        # Single (unknown) NUMA node spanning every logical cpu this host
        # reports -- the same conservative fallback as cpu_provenance.
        all_logical = set(range(topology["cpu_threads"]))
        numa_node_cpus = {0: all_logical}

    return HostTopology(
        cpu_sockets=topology["cpu_sockets"],
        physical_cores=topology["cpu_physical_cores"],
        logical_threads=topology["cpu_threads"],
        physical_cpu_ids=physical_cpu_ids,
        numa_node_cpus=numa_node_cpus,
        metadata_incomplete=topology["metadata_incomplete"],
    )


def validate_thread_count(topology, threads, *, allow_smt=False):
    """Raise :class:`OmpTopologyError` if `threads` would oversubscribe the node.

    `allow_smt=False` (the standard, physical-cores-first campaign) caps at
    `topology.physical_cores`; `allow_smt=True` (the separately-enabled SMT
    extension) caps at `topology.logical_threads`. Never silently clamps --
    an illegal request is refused, not corrected.
    """
    threads = int(threads)
    if threads < 1:
        raise OmpTopologyError(f"thread count must be >= 1, got {threads}")
    ceiling = topology.logical_threads if allow_smt else topology.physical_cores
    if threads > ceiling:
        kind = "logical threads" if allow_smt else "physical cores"
        hint = (
            "reduce the requested thread count"
            if allow_smt or threads > topology.logical_threads
            else "enable the SMT extension explicitly (smt_mode: smt_extension) to use hardware threads"
        )
        raise OmpTopologyError(
            f"requested {threads} OpenMP threads would oversubscribe this host's "
            f"{ceiling} available {kind}; {hint}"
        )
    return threads


def validate_thread_sweep(thread_counts, topology, *, allow_smt=False):
    """Validate a whole thread-count sweep against real hardware.

    Rejects an empty sweep, duplicate thread counts and any count that would
    oversubscribe the node (see `validate_thread_count`). Returns the sweep
    as a sorted list of ints; never reorders silently past validation, never
    drops an illegal entry instead of raising.
    """
    thread_counts = list(thread_counts)
    if not thread_counts:
        raise OmpTopologyError("thread sweep must declare at least one thread count")
    as_ints = [int(t) for t in thread_counts]
    unique_sorted = sorted(set(as_ints))
    if len(unique_sorted) != len(as_ints):
        raise OmpTopologyError(f"thread_counts must be unique, got {thread_counts}")
    for threads in unique_sorted:
        validate_thread_count(topology, threads, allow_smt=allow_smt)
    return unique_sorted


def physical_core_thread_counts(topology, *, max_threads=None):
    """Thread counts (1..physical_cores, capped at `max_threads`) legal for
    the standard physical-cores-first campaign on this host."""
    ceiling = topology.physical_cores if max_threads is None else min(max_threads, topology.physical_cores)
    return list(range(1, ceiling + 1))


def expected_cpu_ids_for_threads(topology, threads, *, proc_bind="close", allow_smt=False):
    """Model which logical cpu ids `OMP_PLACES=cores` +
    `OMP_PROC_BIND=<proc_bind>` will occupy for `threads` threads.

    This is a documented *placement model*, not a runtime measurement: it
    assumes the process is not further constrained by an external affinity
    mask below the full host, and it orders places ascending by logical cpu
    id. The standard, physical-cores-first campaign (`allow_smt=False`, the
    default) draws places only from `topology.physical_cpu_ids` -- one
    representative per physical core, matching `OMP_PLACES=cores`'s own
    intent to avoid SMT siblings. `allow_smt=True` (the separately-enabled
    SMT extension) draws from every logical cpu id
    (`topology.logical_cpu_ids`), since an SMT-extension sweep may request
    more threads than there are physical cores. `close` packs the first
    `threads` places in order; `spread` distributes `threads` places as
    evenly as the ordered place list allows. Use `classify_observed_affinity`
    instead when a real, measured process affinity mask is available (e.g.
    from `cpu_provenance.capture_process_affinity`).
    """
    threads = int(threads)
    places = topology.logical_cpu_ids if allow_smt else topology.physical_cpu_ids
    if threads > len(places):
        kind = "logical cpu" if allow_smt else "physical-core"
        raise OmpTopologyError(
            f"requested {threads} OpenMP threads exceeds {len(places)} known {kind} places"
        )
    if proc_bind == "close":
        return list(places[:threads])
    if proc_bind == "spread":
        step = len(places) / threads
        chosen = sorted({places[int(index * step)] for index in range(threads)})
        return chosen
    raise OmpTopologyError(f"unsupported OMP_PROC_BIND {proc_bind!r}; expected 'close' or 'spread'")


def classify_numa_placement(topology, cpu_ids):
    """Classify a set of occupied logical cpu ids as one NUMA node/socket,
    several, or the full node (blueprint section 9/D).

    A single-socket host always classifies as `NUMA_SINGLE_NODE`: "full node"
    is meaningful only once more than one NUMA node exists to span.
    """
    cpu_ids = frozenset(int(c) for c in cpu_ids)
    if not cpu_ids:
        raise OmpTopologyError("cannot classify NUMA placement of an empty cpu set")
    spanned_nodes = {
        node for node, members in topology.numa_node_cpus.items() if cpu_ids & members
    }
    if not spanned_nodes:
        raise OmpTopologyError(
            f"cpu ids {sorted(cpu_ids)} do not overlap any known NUMA node "
            f"({topology.numa_node_cpus!r})"
        )
    total_nodes = len(topology.numa_node_cpus)
    if len(spanned_nodes) == 1:
        return NUMA_SINGLE_NODE
    if len(spanned_nodes) >= total_nodes:
        return NUMA_FULL_NODE
    return NUMA_MULTI_NODE


def classify_thread_count_numa(topology, threads, *, proc_bind="close", allow_smt=False):
    """`classify_numa_placement` for the modelled placement of `threads`
    threads under `proc_bind` (see `expected_cpu_ids_for_threads`)."""
    cpu_ids = expected_cpu_ids_for_threads(topology, threads, proc_bind=proc_bind, allow_smt=allow_smt)
    return classify_numa_placement(topology, cpu_ids)


def classify_observed_affinity(topology, affinity_range_string):
    """`classify_numa_placement` for a real, measured process affinity mask
    (`cpu_provenance.capture_process_affinity`'s compact range-string format),
    rather than the `expected_cpu_ids_for_threads` model. Returns `None` when
    no affinity string is available, rather than guessing.
    """
    if not affinity_range_string:
        return None
    cpu_ids = parse_cpu_range_list(affinity_range_string)
    if not cpu_ids:
        return None
    return classify_numa_placement(topology, cpu_ids)
