"""Source and build provenance (WP-04 section A).

Captures what a benchmark record's `git_commit`/`git_dirty`/`git_dirty_files`,
`binary_checksum`, `build_id`, `compiler`/`compiler_version`, `compile_flags`
and `cmake_options` fields need, from real, inspectable evidence:

* git state comes from `git rev-parse`/`git status --porcelain`;
* the requested precision and GPU backend come from the build's own
  `CMakeCache.txt` (`UPPASD_PRECISION`, `UPPASD_GPU_BACKEND`) -- never from
  the executable's filename. A build directory called ``build_gpu_fp32``
  proves nothing on its own; the cache entry is the actual configuration
  CMake recorded.
* compiler identity/version/flags come from the same cache, cross-checked by
  actually invoking the recorded compiler with `--version`.

All external calls (`git`, the compiler) go through an injectable ``run``
callable (defaulting to :func:`subprocess.run`) so tests can supply a mock
without touching the real toolchain -- see WP-04 section G.
"""

from __future__ import annotations

import hashlib
import pathlib
import subprocess

_ZERO_COMMIT = "0" * 40


class SourceProvenanceError(ValueError):
    """Source/build provenance could not be established from real evidence."""


def capture_git_provenance(repo_root=None, *, run=subprocess.run, timeout=10):
    """Return ``{git_commit, git_dirty, git_dirty_files}`` from real `git` state.

    Falls back to an honest all-zero commit and ``git_dirty=False`` (with an
    empty file list) only when `git` itself could not be consulted -- this is
    never presented as a clean tree; callers must still raise
    `metadata_incomplete` when this fallback is hit (see `provenance.py`).
    """
    cwd = str(repo_root) if repo_root is not None else None
    git_commit = _ZERO_COMMIT
    git_dirty = False
    git_dirty_files = []
    incomplete = False
    try:
        commit = run(
            ["git", "rev-parse", "HEAD"], cwd=cwd, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout,
        )
        if commit.returncode == 0 and commit.stdout.strip():
            git_commit = commit.stdout.strip()
        else:
            incomplete = True

        status = run(
            ["git", "status", "--porcelain"], cwd=cwd, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=timeout,
        )
        if status.returncode == 0:
            git_dirty_files = sorted(
                {line[3:] for line in status.stdout.splitlines() if line.strip()}
            )
            git_dirty = bool(git_dirty_files)
        else:
            incomplete = True
    except (OSError, subprocess.TimeoutExpired):
        incomplete = True

    return {
        "git_commit": git_commit,
        "git_dirty": git_dirty,
        "git_dirty_files": git_dirty_files,
        "metadata_incomplete": incomplete,
    }


def compute_binary_checksum(binary_path):
    """SHA-256 of the executable actually run, hex-encoded."""
    digest = hashlib.sha256()
    with open(binary_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_cmake_cache(cache_path):
    """Parse a `CMakeCache.txt` into ``{name: value}``, decoding booleans.

    Only ``NAME:TYPE=VALUE`` entry lines are considered; comments (``#``,
    ``//``) and blank lines are skipped. ``BOOL`` entries decode `ON`/`TRUE`/
    `YES`/`1` (case-insensitively) to `True`, `OFF`/`FALSE`/`NO`/`0`/empty to
    `False`; everything else is kept as its raw string.
    """
    entries = {}
    text = pathlib.Path(cache_path).read_text()
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        if ":" not in line or "=" not in line:
            continue
        name_type, _, value = line.partition("=")
        name, _, vtype = name_type.partition(":")
        name = name.strip()
        vtype = vtype.strip().upper()
        if not name:
            continue
        if vtype == "BOOL":
            entries[name] = value.strip().upper() in ("ON", "TRUE", "YES", "1")
        else:
            entries[name] = value
    return entries


# Cache entries worth carrying into `cmake_options`: the UppASD-specific
# switches plus the handful of generic CMake settings that materially affect
# what got built.
_INTERESTING_CACHE_PREFIXES = ("UPPASD_",)
_INTERESTING_CACHE_KEYS = (
    "CMAKE_BUILD_TYPE",
    "CMAKE_CXX_COMPILER",
    "CMAKE_C_COMPILER",
    "CMAKE_Fortran_COMPILER",
    "CMAKE_CUDA_COMPILER",
    "CMAKE_HIP_COMPILER",
    "CMAKE_CUDA_ARCHITECTURES",
    "CMAKE_HIP_ARCHITECTURES",
    "USE_CUDA",
    "USE_HIP",
)


def select_cmake_options(cache):
    """Filter a full cache dict down to the options worth recording."""
    return {
        name: value
        for name, value in cache.items()
        if name in _INTERESTING_CACHE_KEYS or name.startswith(_INTERESTING_CACHE_PREFIXES)
    }


def resolve_requested_precision(cache):
    """`UPPASD_PRECISION` as actually configured -- `SINGLE`/`DOUBLE`/`MIXED`.

    Raises rather than guessing when the cache does not declare it; a build
    directory's name is never consulted (section A: "Do not infer build
    configuration from executable names alone").
    """
    value = cache.get("UPPASD_PRECISION")
    if value not in ("SINGLE", "DOUBLE", "MIXED"):
        raise SourceProvenanceError(
            f"CMakeCache.txt does not declare a recognized UPPASD_PRECISION (got {value!r})"
        )
    return value


def resolve_backend(cache):
    """`(backend, gpu_backend)` from `UPPASD_GPU_BACKEND` -- `OFF`/`CUDA`/`HIP`."""
    value = cache.get("UPPASD_GPU_BACKEND", "OFF")
    if value in (None, "", "OFF"):
        return "CPU", None
    if value in ("CUDA", "HIP"):
        return "GPU", value
    raise SourceProvenanceError(f"CMakeCache.txt has an unrecognized UPPASD_GPU_BACKEND ({value!r})")


def resolve_compile_flags(cache, *, language):
    """Combine `CMAKE_<language>_FLAGS` with the active build type's flags.

    Returns a flat, whitespace-split list -- what actually reaches the
    compiler for a normal (non-per-target-override) build.
    """
    build_type = cache.get("CMAKE_BUILD_TYPE") or ""
    base = cache.get(f"CMAKE_{language}_FLAGS", "") or ""
    per_type = cache.get(f"CMAKE_{language}_FLAGS_{build_type.upper()}", "") if build_type else ""
    combined = f"{base} {per_type}".split()
    return combined


def resolve_compiler(cache, *, language, run=subprocess.run, timeout=10):
    """`(compiler, compiler_version)` for `CMAKE_<language>_COMPILER`.

    The compiler is actually invoked with `--version`; the first output line
    is kept verbatim rather than parsed into a bare version number, since
    compiler `--version` formats are not uniform across toolchains.
    """
    compiler_path = cache.get(f"CMAKE_{language}_COMPILER")
    if not compiler_path:
        raise SourceProvenanceError(f"CMakeCache.txt has no CMAKE_{language}_COMPILER entry")
    compiler = pathlib.Path(compiler_path).name
    try:
        result = run(
            [compiler_path, "--version"], text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout,
        )
        first_line = (result.stdout or "").splitlines()[0].strip() if result.stdout else ""
    except (OSError, subprocess.TimeoutExpired, IndexError):
        first_line = ""
    return compiler, (first_line or "unknown")


def build_source_context(binary_path, build_dir, *, repo_root=None, language="Fortran", run=subprocess.run):
    """Assemble the section-A provenance block for one (binary, build).

    ``build_dir`` must be the CMake build directory the binary was produced
    in (the one holding its `CMakeCache.txt`); ``language`` selects which
    per-language compiler/flags identify the build (UppASD's driver language
    is Fortran, so that is the default across CPU and GPU backends).

    Returns a dict with keys ``git_commit``, ``git_dirty``,
    ``git_dirty_files``, ``binary_checksum``, ``build_id``, ``compiler``,
    ``compiler_version``, ``compile_flags``, ``cmake_options``,
    ``requested_precision``, ``precision_support_state``, ``backend``,
    ``gpu_backend``, plus ``metadata_incomplete`` (True if any part of this
    had to fall back rather than reading real evidence).
    """
    build_dir = pathlib.Path(build_dir)
    cache_path = build_dir / "CMakeCache.txt"
    incomplete = not cache_path.is_file()

    binary_checksum = compute_binary_checksum(binary_path)
    git_info = capture_git_provenance(repo_root, run=run)
    incomplete = incomplete or git_info["metadata_incomplete"]

    if cache_path.is_file():
        cache = parse_cmake_cache(cache_path)
        requested_precision = resolve_requested_precision(cache)
        backend, gpu_backend = resolve_backend(cache)
        compiler, compiler_version = resolve_compiler(cache, language=language, run=run)
        compile_flags = resolve_compile_flags(cache, language=language)
        cmake_options = select_cmake_options(cache)
    else:
        requested_precision = "DOUBLE"
        backend, gpu_backend = "CPU", None
        compiler, compiler_version = "unaudited", "unaudited"
        compile_flags = []
        cmake_options = {}

    return {
        "git_commit": git_info["git_commit"],
        "git_dirty": git_info["git_dirty"],
        "git_dirty_files": git_info["git_dirty_files"],
        "binary_checksum": binary_checksum,
        "build_id": f"{build_dir.name}-{binary_checksum[:12]}",
        "compiler": compiler,
        "compiler_version": compiler_version,
        "compile_flags": compile_flags,
        "cmake_options": cmake_options,
        "requested_precision": requested_precision,
        "precision_support_state": "unsupported" if requested_precision == "MIXED" else "unaudited",
        "backend": backend,
        "gpu_backend": gpu_backend,
        "metadata_incomplete": incomplete,
    }
