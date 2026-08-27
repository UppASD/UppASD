"""Case manifest loading and immutable-template handling (WP-02).

A *case* is a physical input/model family: an immutable production input
template plus a manifest declaring its legal variants, size ladder,
allow-listed input overrides and workload-metadata hook (see
benchmarks/cases/README.md and TERMINOLOGY.md).

This module never invents case content. It loads and validates manifests
written elsewhere, and turns (case, variant, size) into a generated work
directory without ever touching the tracked template.
"""

from __future__ import annotations

import hashlib
import json
import pathlib
import shutil

import jsonschema
import yaml

from harness import schema_validate

CASE_SCHEMA_PATH = schema_validate.SCHEMA_DIR / "case_manifest.v1.schema.json"

# The harness-wide vocabulary of input keywords a case manifest is ever
# permitted to declare in `allowed_input_overrides`. A case's own list only
# narrows this set further; it can never widen it (enforced in
# _validate_manifest). Deliberately absent: exchange/DMI parameters,
# interaction cutoffs, lattice geometry (cell/BC/Sym), posfile/momfile/
# exchange file paths, and moment magnitudes -- the framework must never
# silently alter any of those.
#
# `do_prnstruct` is included even though the blueprint's illustrative list
# (supercell replication, Nstep, seed, temperature, measurement cadence)
# does not name it: it requests UppASD's own struct.<simid>.out diagnostic
# dump (source/Hamiltonian/printhamiltonian.f90::prn_exchange), which is the
# real production output the E. workload-metadata neighbour parser reads.
# It changes no Hamiltonian, lattice or moment parameter.
#
# `gpu_mode` is included for the same reason `do_prnstruct` is: it is a
# runtime backend-dispatch switch, not a Hamiltonian/lattice/moment
# parameter. source/Input/inputhandler.f90's `gpu_mode` case sets `do_gpu`;
# source/uppasd.f90's `run_initial_phase`/`run_measurement_phase` dispatch
# to the GPU driver (`sd_iphaseGPU`/`sd_mphaseGPU`) only when `do_gpu=='Y'`,
# else the ordinary Fortran CPU path -- and on a build without CUDA/HIP
# compiled in, that GPU driver is `source/Tools/nocuda.f90`'s stub module,
# whose subroutines silently `return` with no work done. A case template
# shared between CPU and GPU backends (as the same case_id/variant_id/
# size_id is, across backend records) therefore cannot bake in one fixed
# `gpu_mode`: baking in 1 would make a CPU-only build silently skip the
# simulation it was asked to run, and baking in 0 would make a GPU build
# never dispatch to the device. It must be a per-run override instead --
# a WP-06-style GPU campaign passes `extra_overrides={"gpu_mode": 1}` for
# its GPU samples, leaving the template's own baseline (0) for CPU samples.
GLOBALLY_SAFE_OVERRIDE_KEYS = frozenset(
    {
        "ncell",  # supercell replication -- the sanctioned size axis
        "Nstep",  # measurement-phase step count
        "tseed",  # RNG seed
        "temp",  # production-phase temperature, when a case variant defines it
        "avrg_step",
        "cumu_step",
        "tottraj_step",
        "ene_step",  # measurement cadence
        "do_prnstruct",  # request struct.<simid>.out for neighbour metadata
        "gpu_mode",  # runtime CPU/GPU backend dispatch, not a physics parameter
    }
)

_INPUT_FILENAME = "inpsd.dat"

# Case-level workload_class (blueprint section 3/5) -> the per-record
# workload_class enum in benchmark_record.v1.schema.json. See
# benchmarks/cases/README.md's mapping table; kept here because it is the
# one place both the case manifest and the WP-03 runner (which writes
# records) both depend on.
RECORD_WORKLOAD_CLASS = {
    "short_range_neighbor": "NEIGHBOR_LIST",
    "medium_range_neighbor": "NEIGHBOR_LIST",
    "long_range_neighbor": "NEIGHBOR_LIST",
    "fft_dipole": "FFT_DIPOLE",
    "neighbor_plus_dipole": "NEIGHBOR_LIST_PLUS_FFT_DIPOLE",
}


class CaseManifestError(ValueError):
    """A case manifest, or an operation on a resolved case, is invalid."""


class Case:
    """A loaded, validated case manifest bound to its immutable template."""

    def __init__(self, manifest, manifest_path, template_dir):
        self.manifest = manifest
        self.manifest_path = pathlib.Path(manifest_path)
        self.template_dir = pathlib.Path(template_dir)

    @property
    def id(self):
        return self.manifest["id"]

    @property
    def infrastructure_test_only(self):
        return bool(self.manifest.get("infrastructure_test_only", False))

    @property
    def record_workload_class(self):
        """This case's `workload_class`, mapped onto the per-record enum."""
        return RECORD_WORKLOAD_CLASS[self.manifest["workload_class"]]

    def resolve_variant(self, variant_id):
        for variant in self.manifest["variants"]:
            if variant["id"] == variant_id:
                return variant
        raise CaseManifestError(f"case {self.id!r} has no variant {variant_id!r}")

    def resolve_size(self, size_id):
        for size in self.manifest["sizes"]:
            if size["id"] == size_id:
                return size
        raise CaseManifestError(f"case {self.id!r} has no size {size_id!r}")


class RunDirectory:
    """A generated, writable per-run work directory produced from a case."""

    def __init__(self, path, case, variant_id, size_id, template_fingerprint, overrides):
        self.path = pathlib.Path(path)
        self.case = case
        self.variant_id = variant_id
        self.size_id = size_id
        self.template_fingerprint = template_fingerprint
        self.overrides = overrides


def _load_case_schema():
    with open(CASE_SCHEMA_PATH) as handle:
        return json.load(handle)


def load_case_manifest(path):
    """Load and validate a case manifest. Returns a :class:`Case`.

    Accepts YAML (`.yaml`/`.yml`) or JSON (`.json`). Raises
    :class:`CaseManifestError` on any schema or cross-field violation.
    """
    path = pathlib.Path(path)
    with open(path) as handle:
        if path.suffix in (".yaml", ".yml"):
            manifest = yaml.safe_load(handle)
        elif path.suffix == ".json":
            manifest = json.load(handle)
        else:
            raise CaseManifestError(f"{path}: unsupported manifest extension {path.suffix!r}")

    if not isinstance(manifest, dict):
        raise CaseManifestError(f"{path}: manifest must be a mapping")

    _validate_manifest(manifest, path)

    template_dir = (path.parent / manifest["template_directory"]).resolve()
    if not template_dir.is_dir():
        raise CaseManifestError(f"{path}: template_directory does not exist: {template_dir}")

    return Case(manifest=manifest, manifest_path=path, template_dir=template_dir)


def _validate_manifest(manifest, path):
    schema = _load_case_schema()
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(manifest), key=lambda e: list(map(str, e.absolute_path)))
    if errors:
        messages = "; ".join(
            f"{'/'.join(str(p) for p in e.absolute_path) or '<root>'}: {e.message}" for e in errors
        )
        raise CaseManifestError(f"{path}: {messages}")

    allowed = set(manifest["allowed_input_overrides"])
    illegal = allowed - GLOBALLY_SAFE_OVERRIDE_KEYS
    if illegal:
        raise CaseManifestError(
            f"{path}: allowed_input_overrides declares keys outside the "
            f"harness-wide safe set: {sorted(illegal)}"
        )
    if "ncell" not in allowed:
        raise CaseManifestError(
            f"{path}: allowed_input_overrides must include 'ncell' -- every "
            "case's declared size ladder is applied through it"
        )

    for variant in manifest["variants"]:
        bad = set(variant["overrides"]) - allowed
        if bad:
            raise CaseManifestError(
                f"{path}: variant {variant['id']!r} overrides keys not declared "
                f"in allowed_input_overrides: {sorted(bad)}"
            )

    dimensionality = manifest["dimensionality"]
    expected_len = {"2D": 2, "3D": 3}.get(dimensionality)
    for size in manifest["sizes"]:
        replication = size["replication"]
        if expected_len is not None and len(replication) != expected_len:
            raise CaseManifestError(
                f"{path}: size {size['id']!r} has a {len(replication)}-tuple "
                f"replication but dimensionality {dimensionality!r} requires "
                f"length {expected_len}"
            )
        if expected_len is None and len(replication) != len(manifest["scalable_dimensions"]):
            raise CaseManifestError(
                f"{path}: size {size['id']!r} replication length "
                f"{len(replication)} does not match scalable_dimensions "
                f"{manifest['scalable_dimensions']!r}"
            )


def discover_cases(cases_dir):
    """Load every ``case.yaml``/``case.json`` under ``cases_dir``.

    Returns ``{case_id: Case}``. Raises :class:`CaseManifestError` if two
    manifests declare the same ``id``.
    """
    cases_dir = pathlib.Path(cases_dir)
    found = {}
    for pattern in ("*/case.yaml", "*/case.yml", "*/case.json"):
        for manifest_path in sorted(cases_dir.glob(pattern)):
            case = load_case_manifest(manifest_path)
            if case.id in found:
                raise CaseManifestError(
                    f"duplicate case id {case.id!r}: {found[case.id].manifest_path} "
                    f"and {manifest_path}"
                )
            found[case.id] = case
    return found


def filter_cases_for_campaign(cases, *, authoritative):
    """Return the subset of ``cases`` legal for a campaign of this tier.

    Authoritative (FULL) campaigns must never admit an
    ``infrastructure_test_only`` case; it exists solely to exercise this
    framework's own machinery, not to report performance.
    """
    if not authoritative:
        return dict(cases)
    return {cid: case for cid, case in cases.items() if not case.infrastructure_test_only}


def require_not_infrastructure_only(case):
    """Raise if ``case`` may not contribute to an authoritative report."""
    if case.infrastructure_test_only:
        raise CaseManifestError(
            f"case {case.id!r} is infrastructure_test_only and cannot enter "
            "an authoritative campaign report"
        )


def compute_template_fingerprint(template_dir):
    """Deterministic SHA-256 fingerprint of every file under ``template_dir``.

    Depends only on the relative path and byte content of each tracked file,
    so it is stable across machines and reruns and changes if any template
    file is added, removed or edited.
    """
    template_dir = pathlib.Path(template_dir)
    files = sorted(p for p in template_dir.rglob("*") if p.is_file())
    hasher = hashlib.sha256()
    for file_path in files:
        rel = file_path.relative_to(template_dir).as_posix()
        hasher.update(rel.encode("utf-8"))
        hasher.update(b"\0")
        hasher.update(file_path.read_bytes())
        hasher.update(b"\0")
    return hasher.hexdigest()


def replication_to_ncell(case, size):
    manifest = case.manifest
    replication = size["replication"]
    if manifest["dimensionality"] == "2D":
        nx, ny = replication
        return [nx, ny, manifest["thickness_cells"]]
    if manifest["dimensionality"] == "3D":
        return list(replication)
    if len(replication) == 3:
        return list(replication)
    raise CaseManifestError(
        f"case {case.id!r} is case_specific with a {len(replication)}-tuple "
        "replication; ncell folding is only defined for length-3 tuples "
        "(override 'ncell' explicitly via extra_overrides instead)"
    )


def build_overrides(case, variant_id, size_id, *, extra_overrides=None):
    """Compute the full, allow-list-checked override set for one run.

    Combines the variant's declared overrides with the size's replication
    (folded into ``ncell``) and any additional overrides the caller
    requests. Every key -- from any source -- must be declared in the
    case's ``allowed_input_overrides``; anything else raises
    :class:`CaseManifestError` rather than being silently dropped or
    silently applied.
    """
    variant = case.resolve_variant(variant_id)
    size = case.resolve_size(size_id)
    allowed = set(case.manifest["allowed_input_overrides"])

    overrides = dict(variant["overrides"])
    overrides["ncell"] = replication_to_ncell(case, size)

    for key, value in (extra_overrides or {}).items():
        if key not in allowed:
            raise CaseManifestError(
                f"override {key!r} is not in allowed_input_overrides for case "
                f"{case.id!r}"
            )
        overrides[key] = value

    return overrides


def _read_line_keyword(line):
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return None
    return stripped.split(None, 1)[0]


def _format_override_line(key, value):
    if isinstance(value, (list, tuple)):
        return key + " " + " ".join(str(v) for v in value)
    return f"{key} {value}"


def _apply_overrides(work_dir, overrides):
    inpsd_path = work_dir / _INPUT_FILENAME
    if not inpsd_path.is_file():
        raise CaseManifestError(f"generated work directory has no {_INPUT_FILENAME}: {work_dir}")

    remaining = dict(overrides)
    out_lines = []
    for line in inpsd_path.read_text().splitlines():
        keyword = _read_line_keyword(line)
        if keyword is not None and keyword in remaining:
            out_lines.append(_format_override_line(keyword, remaining.pop(keyword)))
        else:
            out_lines.append(line)

    for key, value in remaining.items():
        out_lines.append(_format_override_line(key, value))

    inpsd_path.write_text("\n".join(out_lines) + "\n")


def generate_run_directory(case, variant_id, size_id, work_root, *, run_id=None, extra_overrides=None):
    """Copy ``case``'s template into a fresh work directory and apply overrides.

    The template directory is never executed in and never modified: this
    function only ever reads from it. Returns a :class:`RunDirectory`
    describing the generated directory, the template fingerprint used, and
    the overrides actually applied.
    """
    work_root = pathlib.Path(work_root)
    overrides = build_overrides(case, variant_id, size_id, extra_overrides=extra_overrides)
    fingerprint = compute_template_fingerprint(case.template_dir)

    run_dir_name = run_id or f"{case.id}__{variant_id}__{size_id}"
    run_dir = work_root / run_dir_name
    if run_dir.resolve() == case.template_dir.resolve():
        raise CaseManifestError("refusing to generate a run directory equal to the template directory")
    if run_dir.exists():
        raise CaseManifestError(f"run directory already exists: {run_dir}")

    work_root.mkdir(parents=True, exist_ok=True)
    shutil.copytree(case.template_dir, run_dir)
    _apply_overrides(run_dir, overrides)

    manifest_out = {
        "case_id": case.id,
        "variant_id": variant_id,
        "size_id": size_id,
        "template_directory": str(case.template_dir),
        "template_fingerprint": fingerprint,
        "applied_overrides": overrides,
    }
    (run_dir / "_benchmark_run.json").write_text(json.dumps(manifest_out, indent=2, sort_keys=True) + "\n")

    return RunDirectory(
        path=run_dir,
        case=case,
        variant_id=variant_id,
        size_id=size_id,
        template_fingerprint=fingerprint,
        overrides=overrides,
    )


def read_keyword(run_dir, keyword):
    """Return the first token following ``keyword`` in a run's inpsd.dat, or None."""
    inpsd_path = pathlib.Path(run_dir) / _INPUT_FILENAME
    for line in inpsd_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(None, 1)
        if len(parts) == 2 and parts[0].lower() == keyword.lower():
            return parts[1].split()[0]
    return None


def read_simid(run_dir):
    """Read the ``simid`` keyword out of a generated run directory's inpsd.dat."""
    simid = read_keyword(run_dir, "simid")
    if simid is None:
        raise CaseManifestError(f"{pathlib.Path(run_dir) / _INPUT_FILENAME}: no simid keyword found")
    return simid
