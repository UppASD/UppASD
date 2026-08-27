"""CPU topology, OpenMP environment and frequency provenance (WP-04 section B).

Everything here reads real evidence -- `/proc/cpuinfo`, `/sys/devices/system/`,
`os.sched_getaffinity`, the process environment -- and is honest about what it
could not determine rather than guessing: unreadable topology comes back as
`None` fields, and the caller is expected to raise `cpu_affinity_unknown` /
`cpu_frequency_unstable` / `metadata_incomplete` accordingly (see
`provenance.py`). Every filesystem root is a parameter so tests can point at a
synthetic tree instead of the real machine (WP-04 section G).
"""

from __future__ import annotations

import os
import pathlib
import re


def parse_cpuinfo(text):
    """Parse `/proc/cpuinfo` text into `(model_name, sockets, physical_cores, logical_threads)`.

    Physical cores are counted as the number of distinct `(physical id, core
    id)` pairs actually reported; on machines without those keys (some
    virtualized/ARM environments) every logical processor is counted as its
    own physical core, which is the honest conservative fallback.
    """
    blocks = re.split(r"\n\s*\n", text.strip())
    model_name = None
    sockets = set()
    cores = set()
    logical = 0
    for block in blocks:
        fields = {}
        for line in block.splitlines():
            if ":" not in line:
                continue
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
        if "processor" not in fields:
            continue
        logical += 1
        if model_name is None and "model name" in fields:
            model_name = fields["model name"]
        physical_id = fields.get("physical id", "0")
        core_id = fields.get("core id", fields["processor"])
        sockets.add(physical_id)
        cores.add((physical_id, core_id))
    return model_name, (len(sockets) or 1), (len(cores) or logical or 1), logical


def capture_cpu_topology(*, cpuinfo_path="/proc/cpuinfo", node_root="/sys/devices/system/node"):
    """Return `{cpu_model, cpu_sockets, cpu_physical_cores, cpu_threads, numa_nodes}`.

    Any field that cannot be read from real evidence comes back `None`
    (`cpu_model`) or falls back to a documented conservative default of `1`
    (the count fields) rather than raising -- callers decide whether that
    warrants `metadata_incomplete`.
    """
    cpuinfo_path = pathlib.Path(cpuinfo_path)
    if cpuinfo_path.is_file():
        model, sockets, physical_cores, threads = parse_cpuinfo(cpuinfo_path.read_text())
        incomplete = model is None
    else:
        model, sockets, physical_cores, threads = None, 1, os.cpu_count() or 1, os.cpu_count() or 1
        incomplete = True

    node_root = pathlib.Path(node_root)
    if node_root.is_dir():
        numa_nodes = sum(1 for entry in node_root.iterdir() if re.fullmatch(r"node\d+", entry.name)) or 1
    else:
        numa_nodes = 1
        incomplete = True

    return {
        "cpu_model": model or "unaudited",
        "cpu_sockets": sockets,
        "cpu_physical_cores": physical_cores,
        "cpu_threads": threads,
        "numa_nodes": numa_nodes,
        "metadata_incomplete": incomplete,
    }


_TRUE_TOKENS = {"true", "1", "yes"}
_FALSE_TOKENS = {"false", "0", "no"}


def capture_omp_environment(env=None, *, default_threads=1):
    """Return `{omp_threads, omp_places, omp_proc_bind, omp_dynamic}` from `env`.

    `omp_threads` falls back to `default_threads` when `OMP_NUM_THREADS` is
    unset -- blueprint section 9 requires campaigns to set it explicitly, so
    this default only ever backstops an under-specified developer run.
    `omp_places`/`omp_proc_bind` are `None` when unset (not `"unset"`, so the
    schema's null-means-not-configured distinction holds); `omp_dynamic` is
    `None` unless `OMP_DYNAMIC` parses unambiguously as a boolean.
    """
    env = os.environ if env is None else env
    threads_raw = env.get("OMP_NUM_THREADS")
    try:
        omp_threads = int(threads_raw) if threads_raw else default_threads
    except ValueError:
        omp_threads = default_threads

    dynamic_raw = env.get("OMP_DYNAMIC")
    if dynamic_raw is None:
        omp_dynamic = None
    else:
        token = dynamic_raw.strip().lower()
        if token in _TRUE_TOKENS:
            omp_dynamic = True
        elif token in _FALSE_TOKENS:
            omp_dynamic = False
        else:
            omp_dynamic = None

    return {
        "omp_threads": omp_threads,
        "omp_places": env.get("OMP_PLACES"),
        "omp_proc_bind": env.get("OMP_PROC_BIND"),
        "omp_dynamic": omp_dynamic,
    }


def _format_affinity_ranges(cpu_ids):
    ids = sorted(cpu_ids)
    if not ids:
        return ""
    ranges = []
    start = prev = ids[0]
    for cpu_id in ids[1:]:
        if cpu_id == prev + 1:
            prev = cpu_id
            continue
        ranges.append(f"{start}" if start == prev else f"{start}-{prev}")
        start = prev = cpu_id
    ranges.append(f"{start}" if start == prev else f"{start}-{prev}")
    return ",".join(ranges)


def capture_process_affinity():
    """Observed CPU affinity of the current process as a compact range string.

    Returns `None` (not raising) when the platform has no affinity concept
    (`os.sched_getaffinity` is Linux-only) -- schema semantics: `None` here
    is what raises `cpu_affinity_unknown` upstream.
    """
    getter = getattr(os, "sched_getaffinity", None)
    if getter is None:
        return None
    try:
        return _format_affinity_ranges(getter(0))
    except OSError:
        return None


def capture_cpu_frequency_state(cpu_ids, *, sysfs_root="/sys/devices/system/cpu"):
    """Best-effort per-core `{governor, current_freq_khz}` for `cpu_ids`.

    Missing `cpufreq` (containers, some VMs, non-Linux) is reported as an
    empty result rather than raising; a caller comparing two snapshots
    (before/after a run) is how `cpu_frequency_unstable` gets raised, not
    this function alone.
    """
    sysfs_root = pathlib.Path(sysfs_root)
    governors = {}
    frequencies_khz = {}
    for cpu_id in cpu_ids:
        cpufreq_dir = sysfs_root / f"cpu{cpu_id}" / "cpufreq"
        governor_path = cpufreq_dir / "scaling_governor"
        freq_path = cpufreq_dir / "scaling_cur_freq"
        try:
            governors[cpu_id] = governor_path.read_text().strip()
        except OSError:
            pass
        try:
            frequencies_khz[cpu_id] = int(freq_path.read_text().strip())
        except (OSError, ValueError):
            pass
    return {"governors": governors, "current_freq_khz": frequencies_khz}


def capture_background_load():
    """1-minute load average and machine logical CPU count, or `(None, None)`.

    `os.getloadavg()` is POSIX-only; callers normalize by CPU count to decide
    `background_load_high` (see `provenance.classify_background_load`).
    """
    try:
        load1, _load5, _load15 = os.getloadavg()
    except OSError:
        return None, None
    return load1, (os.cpu_count() or 1)
