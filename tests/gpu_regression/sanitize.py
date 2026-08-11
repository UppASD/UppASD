#!/usr/bin/env python3
"""Run CUDA compute-sanitizer over GPU regression cases.

Two independent targets share this wrapper:

1. The fast-copy measurement path (F8 / audit C8, the default target), whose
   event and host-callback choreography spans workStream and copyStream and is
   therefore not covered by output comparison alone.  Note that `USE_FAST_COPY`
   is not exposed as a CMake option, so the fast path only exists in a build
   configured with the define, e.g.:

    cmake -S . -B build_fastcopy -DCMAKE_BUILD_TYPE=DEBUG -DUSE_CUDA=ON \
        -DUSE_MKL=ON -DCMAKE_CUDA_FLAGS=-DUSE_FAST_COPY \
        -DCMAKE_CXX_FLAGS=-DUSE_FAST_COPY
    cmake --build build_fastcopy --parallel 8

   Fixture preparation is shared with run.py so the sanitized inputs match what
   the regression suite compares.  Only the GPU side is run: the CPU reference
   has no device work to sanitize.

2. RCG-05E's `--dilation-kernel` target: `dilateAdaptiveState`
   (`source/gpu_files/gpuAdaptiveRuntime.cpp`), sanitized directly through the
   `gpu_adaptive_runtime_tests` binary (which already launches it with nonzero
   per-axis buffer widths) via `--kernel-name kernel_substring=...`, no
   fixture/case machinery required.  This target is *not* covered by (1):
   `--case`/`--cases-file`/`DEFAULT_CASES` and this module's fast-copy
   docstring above describe only the USE_FAST_COPY path, and doing nothing
   further would have left the dilation kernel merely assumed clean rather
   than actually checked, which is exactly what RCG-05E's governing rule
   against clearing a race from source reading alone forbids.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import run as regression

ROOT = regression.ROOT

# racecheck is the tool named by the F8 acceptance criteria; the others catch the
# adjacent failure modes that the fast-copy path can produce (stale pinned reads,
# use-after-free of staging tensors, missing stream synchronization).
DEFAULT_TOOLS = ("racecheck", "memcheck", "synccheck", "initcheck")

# Cases that actually exercise measurement/correlation copy-out.  Keep this list
# short: sanitized runs are one to two orders of magnitude slower than native.
DEFAULT_CASES = ("bcc_scalar", "bcc_scalar_gpu_measurements", "bcc_convolution_gpu_measurements")

# RCG-05E: the dilation-kernel target's default binary and kernel filter.
# gpu_adaptive_runtime_tests is a plain CUDA/HIP test build (no USE_FAST_COPY
# define needed); its CG-09/CG-10 test drives GpuAdaptiveRuntime::
# proposeSelectorState with a nonzero per-axis bufferDilationBlocks, which is
# the only call site that ever launches dilateAdaptiveState
# (gpuAdaptiveRuntime.cpp:proposeSelectorState).
DEFAULT_DILATION_BINARY = ROOT / "build_gpu/bin/gpu_adaptive_runtime_tests"
DILATION_KERNEL_NAME = "dilateAdaptiveState"

ERROR_SUMMARY = re.compile(r"ERROR SUMMARY:\s+(\d+)\s+error", re.IGNORECASE)


def find_sanitizer(explicit: str | None) -> str:
    if explicit:
        return explicit
    found = shutil.which("compute-sanitizer")
    if found:
        return found
    # nvcc on PATH implies a toolkit layout with the sanitizer beside it.
    nvcc = shutil.which("nvcc")
    if nvcc:
        candidate = Path(nvcc).with_name("compute-sanitizer")
        if candidate.exists():
            return str(candidate)
    raise SystemExit("compute-sanitizer not found; pass --sanitizer with an explicit path")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--binary", type=Path,
                        help="GPU binary to sanitize (default: the USE_FAST_COPY build, or "
                             "gpu_adaptive_runtime_tests with --dilation-kernel)")
    parser.add_argument("--case", action="append", dest="cases", metavar="NAME",
                        help="case name from cases.json; repeatable (default: measurement cases); "
                             "not used with --dilation-kernel")
    parser.add_argument("--tool", action="append", dest="tools", metavar="TOOL",
                        help=f"compute-sanitizer tool; repeatable (default: {' '.join(DEFAULT_TOOLS)})")
    parser.add_argument("--steps", type=int,
                        help="override Nstep; lower values keep sanitized runs tractable; "
                             "not used with --dilation-kernel")
    parser.add_argument("--sanitizer", help="path to compute-sanitizer")
    parser.add_argument("--cases-file", type=Path, default=regression.DEFAULT_CASES)
    parser.add_argument("--workdir", type=Path, help="reuse this directory instead of a temp one")
    parser.add_argument("--keep", action="store_true", help="keep the work directory on success")
    parser.add_argument("--timeout", type=int, default=3600, help="per-run timeout in seconds")
    parser.add_argument("--dilation-kernel", action="store_true",
                        help="RCG-05E target: sanitize dilateAdaptiveState directly via "
                             "gpu_adaptive_runtime_tests, instead of the fast-copy fixture path")
    parser.add_argument("--kernel-name", default=DILATION_KERNEL_NAME,
                        help=f"compute-sanitizer --kernel-name kernel_substring filter used with "
                             f"--dilation-kernel (default: {DILATION_KERNEL_NAME})")
    return parser.parse_args()


def sanitize(sanitizer: str, tool: str, binary: Path, directory: Path, timeout: int,
             kernel_name: str | None = None) -> tuple[int, str]:
    """Run one case under one tool.  Returns (error_count, log_path)."""
    command = [sanitizer, "--tool", tool, "--error-exitcode", "1"]
    if tool == "racecheck":
        # Report the racing accesses, not just that a race occurred.
        command += ["--racecheck-report", "all"]
    if kernel_name:
        command += ["--kernel-name", f"kernel_substring={kernel_name}"]
        if tool == "initcheck":
            # initcheck's host-API cudaMemcpy-source check is not
            # kernel-scoped, so with a --kernel-name filter in effect it
            # still flags unrelated host buffers touched by other code paths
            # in the same process (observed in practice: 65 findings from
            # gpu_adaptive_runtime_tests' unrelated testPolarizationGate(),
            # with no relationship to the filtered kernel). Disabling it
            # keeps this tool's result attributable to the filtered kernel
            # specifically, matching what --kernel-name already promises for
            # every other tool and for initcheck's own device-side checks.
            command += ["--check-api-memory-access", "no"]
    command.append(str(binary))

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    try:
        completed = subprocess.run(command, cwd=directory, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   timeout=timeout)
        output, returncode = completed.stdout, completed.returncode
    except subprocess.TimeoutExpired as expired:
        output = (expired.stdout or "") + f"\n*** timed out after {timeout}s ***\n"
        returncode = -1

    log = directory / f"sanitize-{tool}.log"
    log.write_text(output)

    match = ERROR_SUMMARY.search(output)
    if match:
        return int(match.group(1)), str(log)
    # No summary line means the tool never got far enough to emit one.
    return (0 if returncode == 0 else 1), str(log)


def sanitize_dilation_kernel(args: argparse.Namespace, sanitizer: str) -> int:
    """RCG-05E target: run every requested tool against dilateAdaptiveState
    directly, with no fixture/case machinery -- gpu_adaptive_runtime_tests
    already exercises it with a nonzero per-axis width."""
    binary = (args.binary or DEFAULT_DILATION_BINARY).resolve()
    if not binary.exists():
        raise SystemExit(
            f"binary not found: {binary}\n"
            "Configure and build a CUDA or HIP test build first, e.g.:\n"
            "  cmake -S . -B build_gpu -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA\n"
            "  cmake --build build_gpu --target gpu_adaptive_runtime_tests")
    tools = args.tools or list(DEFAULT_TOOLS)

    if args.workdir:
        root = args.workdir
        root.mkdir(parents=True, exist_ok=True)
    else:
        root = Path(tempfile.mkdtemp(prefix="uppasd-gpu-sanitize-dilation-"))

    print(f"sanitizer:   {sanitizer}")
    print(f"binary:      {binary}")
    print(f"kernel-name: kernel_substring={args.kernel_name}")
    print(f"tools:       {' '.join(tools)}")

    failed = 0
    try:
        for tool in tools:
            directory = root / f"dilateAdaptiveState-{tool}"
            if directory.exists():
                shutil.rmtree(directory)
            directory.mkdir(parents=True)
            errors, log = sanitize(sanitizer, tool, binary, directory, args.timeout,
                                   kernel_name=args.kernel_name)
            if errors:
                failed += 1
                print(f"FAIL dilateAdaptiveState [{tool}]: {errors} error(s); see {log}")
            else:
                print(f"PASS dilateAdaptiveState [{tool}]: 0 errors")
        return 1 if failed else 0
    finally:
        if args.workdir or args.keep or failed:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


def main() -> int:
    args = parse_args()
    sanitizer = find_sanitizer(args.sanitizer)

    if args.dilation_kernel:
        return sanitize_dilation_kernel(args, sanitizer)

    if args.binary is None:
        args.binary = ROOT / "build_fastcopy/bin/sd.f95.cuda"
    # Runs execute in the case directory, so a relative --binary must be
    # resolved against the invocation directory before we change into it.
    args.binary = args.binary.resolve()
    if not args.binary.exists():
        raise SystemExit(f"binary not found: {args.binary}\n"
                         "Configure a USE_FAST_COPY build first (see this script's docstring).")

    defaults, all_cases = regression.load_cases(args.cases_file)
    wanted = args.cases or list(DEFAULT_CASES)
    by_name = {case["name"]: case for case in all_cases}
    unknown = [name for name in wanted if name not in by_name]
    if unknown:
        raise SystemExit(f"unknown case(s): {', '.join(unknown)}\n"
                         f"available: {', '.join(by_name)}")
    cases = [by_name[name] for name in wanted]
    tools = args.tools or list(DEFAULT_TOOLS)

    if args.workdir:
        root = args.workdir
        root.mkdir(parents=True, exist_ok=True)
    else:
        root = Path(tempfile.mkdtemp(prefix="uppasd-gpu-sanitize-"))

    print(f"sanitizer: {sanitizer}")
    print(f"binary:    {args.binary}")
    print(f"tools:     {' '.join(tools)}")

    failed = 0
    try:
        for case in cases:
            case = dict(case)
            if args.steps is not None:
                case["steps"] = args.steps
            for tool in tools:
                directory = root / f"{case['name']}-{tool}"
                if directory.exists():
                    shutil.rmtree(directory)
                try:
                    regression.prepare_run(case, defaults, gpu=True, destination=directory)
                    errors, log = sanitize(sanitizer, tool, args.binary, directory, args.timeout)
                except Exception as error:  # Keep going so one run gives a full table.
                    failed += 1
                    print(f"FAIL {case['name']} [{tool}]: {error}")
                    continue
                if errors:
                    failed += 1
                    print(f"FAIL {case['name']} [{tool}]: {errors} error(s); see {log}")
                else:
                    print(f"PASS {case['name']} [{tool}]: 0 errors")
        return 1 if failed else 0
    finally:
        if args.workdir or args.keep or failed:
            print(f"work directory: {root}")
        else:
            shutil.rmtree(root, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
