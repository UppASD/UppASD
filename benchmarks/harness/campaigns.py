"""Campaign manifest loading and cell resolution (WP-09).

A campaign manifest (`schema/campaign_manifest.v1.schema.json`) declares
*what* to run -- which cases/variants/sizes, which backends and precisions,
how many samples, under which environment-quality policy -- and leaves
*how* to run it to the harness already built by WP-01 through WP-07:
`cases.py` (case/variant/size resolution), `omp_sweep.py`/`omp_topology.py`
(OpenMP run configurations), `runner.py`/`provenance.py` (execution and
strict-mode gating), `gpu_campaign.py`/`gpu_memory.py` (GPU precision and
memory admission). This module never executes anything itself.

Two ideas this module keeps deliberately separate, per the blueprint's
"keep machine hardware policy separate from physical case manifests" (WP-09
prompt section E):

* the campaign manifest names a `cpu.thread_policy` -- `LEAN_SAMPLE` or
  `FULL_PHYSICAL_SWEEP` -- never a concrete thread-count list. It is
  portable across any host and carries no assumption about how many cores
  a given machine has.
* `resolve_cpu_sweep` turns that policy into a real `omp_sweep.OmpSweep`
  only once a real `omp_topology.HostTopology` is available, using exactly
  the same live-queried-never-hardcoded discipline `omp_topology.py` already
  uses throughout. No thread-count list is ever checked into a campaign
  manifest.

`resolve_cells` is the "runner parsing" CI can exercise without executing
anything: given a campaign and the case manifests it names, it expands
`variant_ids`/`size_ids` (including the `"ALL"` wildcard) and the CPU/GPU
backend axis into the flat list of `(case, variant_id, size_id, backend,
...)` cells a real campaign driver (WP-10) would iterate over.
"""

from __future__ import annotations

import json
import pathlib

import jsonschema
import yaml

from harness import cases as cases_mod
from harness import omp_sweep
from harness import omp_topology
from harness import schema_validate

CAMPAIGN_SCHEMA_PATH = schema_validate.SCHEMA_DIR / "campaign_manifest.v1.schema.json"

LEAN = "LEAN"
FULL = "FULL"

CPU_ONLY = "CPU_ONLY"
CUDA_ONLY = "CUDA_ONLY"
CPU_AND_CUDA = "CPU_AND_CUDA"

LEAN_SAMPLE = "LEAN_SAMPLE"
FULL_PHYSICAL_SWEEP = "FULL_PHYSICAL_SWEEP"

DEVELOPER = "DEVELOPER"
STRICT = "STRICT"

LEAN_REPORT_BANNER = "LEAN CAMPAIGN — NOT AUTHORITATIVE FOR CROSSOVER CLAIMS"

_ALL = "ALL"


class CampaignManifestError(ValueError):
    """A campaign manifest, or an operation on a resolved campaign, is invalid."""


class Campaign:
    """A loaded, schema-validated campaign manifest."""

    def __init__(self, manifest):
        self.manifest = manifest

    @property
    def id(self):
        return self.manifest["id"]

    @property
    def tier(self):
        return self.manifest["tier"]

    @property
    def authoritative(self):
        """FULL campaigns are the blueprint's authoritative tier; LEAN never is."""
        return self.tier == FULL

    @property
    def backend_scope(self):
        return self.manifest["backend_scope"]

    @property
    def has_cpu(self):
        return self.backend_scope in (CPU_ONLY, CPU_AND_CUDA)

    @property
    def has_gpu(self):
        return self.backend_scope in (CUDA_ONLY, CPU_AND_CUDA)

    @property
    def cases(self):
        return list(self.manifest["cases"])

    @property
    def sample_count(self):
        return self.manifest["sample_count"]

    @property
    def measurement_profile(self):
        return self.manifest["measurement_profile"]

    @property
    def environment_quality_mode(self):
        return self.manifest["environment_quality_mode"]

    @property
    def require_clean_environment(self):
        """The `runner.run_sample`/`provenance.gate_sample` strict-mode flag
        this campaign's `environment_quality_mode` implies."""
        return self.environment_quality_mode == STRICT

    @property
    def report_banner(self):
        return self.manifest.get("report_banner")

    @property
    def cpu_thread_policy(self):
        return self.manifest["cpu"]["thread_policy"] if self.has_cpu else None

    @property
    def cpu_proc_bind(self):
        return self.manifest["cpu"]["proc_bind"] if self.has_cpu else None

    @property
    def gpu_precisions(self):
        return list(self.manifest["gpu"]["precisions"]) if self.has_gpu else []

    @property
    def gpu_index(self):
        return self.manifest["gpu"].get("gpu_index", 0) if self.has_gpu else None

    @property
    def single_campaign_execution(self):
        """True only for a `CPU_AND_CUDA` campaign -- the schema-enforced
        acknowledgement that this manifest's CPU and CUDA cells run together
        under one `campaign_id`/build (blueprint section 10)."""
        return bool(self.manifest.get("single_campaign_execution", False))


def _load_campaign_schema():
    with open(CAMPAIGN_SCHEMA_PATH) as handle:
        return json.load(handle)


def load_campaign_manifest(path):
    """Load and schema-validate a campaign manifest. Returns a :class:`Campaign`.

    Accepts YAML (`.yaml`/`.yml`) or JSON (`.json`). Raises
    :class:`CampaignManifestError` on any schema violation. Does not resolve
    case references against `benchmarks/cases/` -- that is `resolve_cells`'s
    job, since it needs the caller's already-discovered case set.
    """
    path = pathlib.Path(path)
    with open(path) as handle:
        if path.suffix in (".yaml", ".yml"):
            manifest = yaml.safe_load(handle)
        elif path.suffix == ".json":
            manifest = json.load(handle)
        else:
            raise CampaignManifestError(f"{path}: unsupported manifest extension {path.suffix!r}")

    if not isinstance(manifest, dict):
        raise CampaignManifestError(f"{path}: manifest must be a mapping")

    schema = _load_campaign_schema()
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(manifest), key=lambda e: list(map(str, e.absolute_path)))
    if errors:
        messages = "; ".join(
            f"{'/'.join(str(p) for p in e.absolute_path) or '<root>'}: {e.message}" for e in errors
        )
        raise CampaignManifestError(f"{path}: {messages}")

    return Campaign(manifest)


def _resolve_id_list_or_all(selector, legal_ids, *, kind, case_id):
    if selector == _ALL:
        return list(legal_ids)
    unknown = [i for i in selector if i not in legal_ids]
    if unknown:
        raise CampaignManifestError(
            f"case {case_id!r}: unknown {kind} id(s) {unknown!r}; legal {kind} ids are {sorted(legal_ids)}"
        )
    return list(selector)


def resolve_cells(campaign, cases_by_id):
    """Expand `campaign` into the flat list of cells a driver would run.

    `cases_by_id` is normally `cases.discover_cases(...)`'s return value. Each
    cell is a dict: `case_id`, `variant_id`, `size_id`, plus `backend` ("CPU"
    or "GPU") and, for a GPU cell, `requested_precision`. CPU cells carry no
    concrete `omp_threads` here -- that is `resolve_cpu_sweep`'s job, since it
    requires a real host topology this function does not take.

    Raises :class:`CampaignManifestError` if a named case/variant/size does
    not exist, or if an authoritative (tier FULL) campaign names an
    `infrastructure_test_only` case (`cases.require_not_infrastructure_only`).
    """
    cells = []
    for selection in campaign.cases:
        case_id = selection["case_id"]
        if case_id not in cases_by_id:
            raise CampaignManifestError(f"campaign {campaign.id!r} names unknown case_id {case_id!r}")
        case = cases_by_id[case_id]
        if campaign.authoritative:
            cases_mod.require_not_infrastructure_only(case)

        legal_variants = {v["id"] for v in case.manifest["variants"]}
        legal_sizes = {s["id"] for s in case.manifest["sizes"]}
        variant_ids = _resolve_id_list_or_all(
            selection["variant_ids"], legal_variants, kind="variant", case_id=case_id
        )
        size_ids = _resolve_id_list_or_all(
            selection["size_ids"], legal_sizes, kind="size", case_id=case_id
        )

        for variant_id in variant_ids:
            for size_id in size_ids:
                if campaign.has_cpu:
                    cells.append({
                        "case_id": case_id, "variant_id": variant_id, "size_id": size_id,
                        "backend": "CPU", "requested_precision": "DOUBLE",
                    })
                if campaign.has_gpu:
                    for precision in campaign.gpu_precisions:
                        cells.append({
                            "case_id": case_id, "variant_id": variant_id, "size_id": size_id,
                            "backend": "GPU", "requested_precision": precision,
                        })
    return cells


def build_lean_sweep_manifest(topology, *, sweep_id="lean_sample", proc_bind="close"):
    """The `LEAN_SAMPLE` thread policy resolved against `topology`: 1 thread,
    one moderate count, and the host's full physical-core count (blueprint
    section 15 "OMP 1 + one moderate + one large thread count"), deduplicated
    for hosts with very few physical cores. Returns an `omp_sweep`-schema-
    shaped dict, not an `OmpSweep` -- see `resolve_cpu_sweep`.
    """
    physical = topology.physical_cores
    moderate = max(1, physical // 2)
    thread_counts = sorted({1, moderate, physical})
    return {
        "schema_version": "1.0.0",
        "id": sweep_id,
        "description": (
            f"LEAN_SAMPLE thread policy resolved against this host: {thread_counts} "
            f"(physical_cores={physical})."
        ),
        "thread_counts": thread_counts,
        "proc_bind": proc_bind,
        "places": "cores",
        "dynamic": False,
        "smt_mode": omp_sweep.PHYSICAL_ONLY,
    }


def build_full_physical_sweep_manifest(topology, *, sweep_id="full_physical_sweep", proc_bind="close"):
    """The `FULL_PHYSICAL_SWEEP` thread policy resolved against `topology`:
    every thread count from 1 to the host's physical-core count (blueprint
    section 15 "broad OpenMP scaling"). Returns an `omp_sweep`-schema-shaped
    dict, not an `OmpSweep` -- see `resolve_cpu_sweep`.
    """
    thread_counts = omp_topology.physical_core_thread_counts(topology)
    return {
        "schema_version": "1.0.0",
        "id": sweep_id,
        "description": (
            f"FULL_PHYSICAL_SWEEP thread policy resolved against this host: "
            f"1..{topology.physical_cores}."
        ),
        "thread_counts": thread_counts,
        "proc_bind": proc_bind,
        "places": "cores",
        "dynamic": False,
        "smt_mode": omp_sweep.PHYSICAL_ONLY,
    }


_SWEEP_BUILDERS = {
    LEAN_SAMPLE: build_lean_sweep_manifest,
    FULL_PHYSICAL_SWEEP: build_full_physical_sweep_manifest,
}


def resolve_cpu_sweep(campaign, topology):
    """`campaign.cpu_thread_policy` resolved against real `topology` as an
    `omp_sweep.OmpSweep`, ready for `omp_sweep.resolve_run_configurations`.

    Raises :class:`CampaignManifestError` if `campaign` has no `cpu` block.
    """
    if not campaign.has_cpu:
        raise CampaignManifestError(f"campaign {campaign.id!r} declares no cpu block")
    builder = _SWEEP_BUILDERS[campaign.cpu_thread_policy]
    manifest = builder(topology, proc_bind=campaign.cpu_proc_bind)
    return omp_sweep.OmpSweep(manifest)
