"""Declarative OpenMP thread-sweep manifests and run-configuration binding
(WP-05 sections A-D).

A sweep manifest (`schema/omp_sweep.v1.schema.json`) declares *intent*: which
thread counts to sweep, one binding policy for the whole sweep, and whether
it is the standard physical-cores-first campaign or a separately-enabled SMT
extension. It says nothing about the host it will run on. This module turns
that declaration into concrete, hardware-validated run configurations: one
`RunConfiguration` per thread count, each carrying the explicit OpenMP
environment (`OMP_NUM_THREADS`/`OMP_PLACES`/`OMP_PROC_BIND`/`OMP_DYNAMIC`,
never left implicit) and a NUMA-placement label, ready to hand to
`runner.run_sample` as its `env` argument.
"""

from __future__ import annotations

import json
import os
import pathlib

import jsonschema
import yaml

from harness import omp_topology
from harness import schema_validate

SWEEP_SCHEMA_PATH = schema_validate.SCHEMA_DIR / "omp_sweep.v1.schema.json"

PHYSICAL_ONLY = "physical_only"
SMT_EXTENSION = "smt_extension"


class OmpSweepError(ValueError):
    """A sweep manifest, or an operation on a resolved sweep, is invalid."""


def _load_sweep_schema():
    with open(SWEEP_SCHEMA_PATH) as handle:
        return json.load(handle)


class OmpSweep:
    """A loaded, schema-validated OpenMP sweep manifest."""

    def __init__(self, manifest):
        self.manifest = manifest

    @property
    def id(self):
        return self.manifest["id"]

    @property
    def thread_counts(self):
        return list(self.manifest["thread_counts"])

    @property
    def proc_bind(self):
        return self.manifest["proc_bind"]

    @property
    def places(self):
        return self.manifest["places"]

    @property
    def smt_mode(self):
        return self.manifest["smt_mode"]

    @property
    def allow_smt(self):
        return self.smt_mode == SMT_EXTENSION


def load_omp_sweep_manifest(path):
    """Load and validate an OpenMP sweep manifest. Returns an :class:`OmpSweep`."""
    path = pathlib.Path(path)
    with open(path) as handle:
        if path.suffix in (".yaml", ".yml"):
            manifest = yaml.safe_load(handle)
        elif path.suffix == ".json":
            manifest = json.load(handle)
        else:
            raise OmpSweepError(f"{path}: unsupported manifest extension {path.suffix!r}")

    if not isinstance(manifest, dict):
        raise OmpSweepError(f"{path}: manifest must be a mapping")

    schema = _load_sweep_schema()
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(manifest), key=lambda e: list(map(str, e.absolute_path)))
    if errors:
        messages = "; ".join(
            f"{'/'.join(str(p) for p in e.absolute_path) or '<root>'}: {e.message}" for e in errors
        )
        raise OmpSweepError(f"{path}: {messages}")

    return OmpSweep(manifest)


def build_omp_env(threads, *, proc_bind, places="cores", base_env=None):
    """Return an environment mapping with every OpenMP binding variable set
    explicitly (blueprint section 9/C): `OMP_NUM_THREADS`, `OMP_PLACES=cores`,
    `OMP_PROC_BIND` (`close` or `spread`) and `OMP_DYNAMIC=FALSE`.

    Never leaves a binding variable to inherit an ambient, unrecorded value:
    `base_env` (default `os.environ`) supplies everything *else* the child
    process needs (`LD_LIBRARY_PATH`, `PATH`, ...), but the four OpenMP
    variables here always come from this call's own arguments.
    """
    if proc_bind not in ("close", "spread"):
        raise OmpSweepError(f"unsupported OMP_PROC_BIND {proc_bind!r}; expected 'close' or 'spread'")
    if places != "cores":
        raise OmpSweepError(f"unsupported OMP_PLACES {places!r}; the standard campaign requires 'cores'")
    threads = int(threads)
    if threads < 1:
        raise OmpSweepError(f"thread count must be >= 1, got {threads}")

    env = dict(os.environ if base_env is None else base_env)
    env["OMP_NUM_THREADS"] = str(threads)
    env["OMP_PLACES"] = places
    env["OMP_PROC_BIND"] = proc_bind
    env["OMP_DYNAMIC"] = "FALSE"
    return env


class RunConfiguration:
    """One hardware-validated `(threads, env, numa placement)` cell of a sweep."""

    def __init__(self, *, threads, env, proc_bind, places, expected_cpu_ids, numa_placement):
        self.threads = threads
        self.env = env
        self.proc_bind = proc_bind
        self.places = places
        self.expected_cpu_ids = expected_cpu_ids
        self.numa_placement = numa_placement

    def __repr__(self):
        return (
            f"RunConfiguration(threads={self.threads!r}, proc_bind={self.proc_bind!r}, "
            f"numa_placement={self.numa_placement!r})"
        )


def resolve_run_configurations(sweep, topology, *, base_env=None):
    """Validate `sweep` against real `topology` and return one
    :class:`RunConfiguration` per thread count.

    Raises :class:`omp_topology.OmpTopologyError` if any declared thread
    count would oversubscribe the node given `sweep.smt_mode` -- a sweep
    manifest never gets to silently clamp itself to fit whatever host it
    lands on.
    """
    validated_counts = omp_topology.validate_thread_sweep(
        sweep.thread_counts, topology, allow_smt=sweep.allow_smt
    )
    configurations = []
    for threads in validated_counts:
        env = build_omp_env(threads, proc_bind=sweep.proc_bind, places=sweep.places, base_env=base_env)
        expected_cpu_ids = omp_topology.expected_cpu_ids_for_threads(
            topology, threads, proc_bind=sweep.proc_bind, allow_smt=sweep.allow_smt
        )
        numa_placement = omp_topology.classify_numa_placement(topology, expected_cpu_ids)
        configurations.append(
            RunConfiguration(
                threads=threads,
                env=env,
                proc_bind=sweep.proc_bind,
                places=sweep.places,
                expected_cpu_ids=expected_cpu_ids,
                numa_placement=numa_placement,
            )
        )
    return configurations
