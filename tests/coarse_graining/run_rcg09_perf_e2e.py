#!/usr/bin/env python3
"""RCG-09 end-to-end performance harness: PERF-ATOMISTIC-PROD vs PERF-CG-SWEEP.

This is the headline measurement.  It runs the *ordinary* UppASD executable
twice on inputs that differ only in whether adaptive coarse graining is
enabled, so "identical supported physics, precision, geometry, ensemble count,
timestep, device and measurement scope" is a property of the inputs rather than
a claim about a test harness.  No internal staging, no direct kernel
invocation, no synthetic work.

Steady state is separated from setup by the two-point method: each case is run
at a short and a long ``Nstep`` and the per-step cost is the difference divided
by the step difference.  Everything that happens once -- process start, CUDA
context creation, input parsing, neighbour-list construction, adaptive topology
upload, output files -- cancels exactly, without needing to trust an internal
timer.

The device is checked before measuring.  RCG-08-FU3 recorded that every RCG-08
timing was taken while another user drove both GPUs to ~90% utilization with
thermal throttling, and required that RCG-09 not repeat that; the check is
refusable rather than advisory for exactly that reason.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PERF_ROOT = Path(__file__).resolve().parent / "perf"

TEXTURES = {
    "spiral": ("spiral_atomistic", "spiral_adaptive"),
    "uniform": ("uniform_atomistic", "uniform_adaptive"),
}


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def median_absolute_deviation(values: list[float]) -> float:
    if not values:
        return 0.0
    center = median(values)
    return median([abs(value - center) for value in values])


def nvidia_smi(query: str) -> list[str]:
    result = subprocess.run(
        ["nvidia-smi", f"--query-gpu={query}", "--format=csv,noheader,nounits"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(f"nvidia-smi failed: {result.stdout}")
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def compute_processes() -> list[str]:
    result = subprocess.run(
        ["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory",
         "--format=csv,noheader"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        return []
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def device_state(device: int) -> dict:
    fields = ("index,name,driver_version,utilization.gpu,memory.used,memory.total,"
              "temperature.gpu,clocks.current.sm,clocks.max.sm,"
              "clocks_event_reasons.active,persistence_mode,compute_mode")
    rows = nvidia_smi(fields)
    if device >= len(rows):
        raise RuntimeError(f"device {device} not present ({len(rows)} reported)")
    parts = [part.strip() for part in rows[device].split(",")]
    keys = ("index", "name", "driver_version", "utilization_percent", "memory_used_mib",
            "memory_total_mib", "temperature_c", "sm_clock_mhz", "sm_clock_max_mhz",
            "throttle_reasons", "persistence_mode", "compute_mode")
    return dict(zip(keys, parts))


def check_device_clean(device: int, max_utilization: float, max_temperature: float,
                       allow_dirty: bool) -> dict:
    state = device_state(device)
    processes = compute_processes()
    complaints = []
    try:
        utilization = float(state["utilization_percent"])
    except ValueError:
        utilization = float("nan")
        complaints.append("utilization unreadable")
    else:
        if utilization > max_utilization:
            complaints.append(
                f"utilization {utilization:.0f}% exceeds {max_utilization:.0f}%")
    try:
        temperature = float(state["temperature_c"])
    except ValueError:
        complaints.append("temperature unreadable")
    else:
        if temperature > max_temperature:
            complaints.append(
                f"temperature {temperature:.0f}C exceeds {max_temperature:.0f}C")
    throttle = state["throttle_reasons"]
    # 0x0 means "not throttled".  Anything else is a real clock event, and on
    # this host it has been thermal slowdown (0x20) caused by another tenant.
    if throttle not in ("0x0000000000000000", "0x0", "Not Active"):
        complaints.append(f"active clock event reasons {throttle}")
    if processes:
        complaints.append(f"{len(processes)} foreign compute process(es): " +
                          "; ".join(processes))
    state["clean"] = not complaints
    state["complaints"] = complaints
    state["foreign_processes"] = processes
    if complaints and not allow_dirty:
        print("rcg09-device-check result=DIRTY", flush=True)
        for complaint in complaints:
            print(f"rcg09-device-check complaint={complaint!r}", flush=True)
        print("rcg09-device-check note=refusing to measure; RCG-08-FU3 requires a "
              "clean device. Re-run with --allow-dirty-device only for smoke "
              "tests, whose numbers are not acceptance evidence.", flush=True)
        raise SystemExit(4)
    print(f"rcg09-device-check result={'CLEAN' if state['clean'] else 'DIRTY_ALLOWED'} "
          f"device={state['name']} driver={state['driver_version']} "
          f"utilization_percent={state['utilization_percent']} "
          f"temperature_c={state['temperature_c']} "
          f"sm_clock_mhz={state['sm_clock_mhz']}/{state['sm_clock_max_mhz']} "
          f"throttle_reasons={throttle} persistence={state['persistence_mode']} "
          f"compute_mode={state['compute_mode']}", flush=True)
    for complaint in state["complaints"]:
        print(f"rcg09-device-check allowed_complaint={complaint!r}", flush=True)
    return state


def cache_value(cache: Path, key: str) -> str:
    if not cache.is_file():
        return ""
    pattern = re.compile(rf"^{re.escape(key)}:[^=]*=(.*)$")
    for line in cache.read_text(errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            return match.group(1).strip()
    return ""


def build_metadata(binary: Path) -> dict:
    cache = binary.parent.parent / "CMakeCache.txt"
    return {
        "binary": str(binary),
        "build_dir": str(cache.parent),
        "backend": cache_value(cache, "UPPASD_GPU_BACKEND"),
        "precision": cache_value(cache, "UPPASD_PRECISION") or "DOUBLE",
        "build_type": cache_value(cache, "CMAKE_BUILD_TYPE"),
        "cuda_architectures": cache_value(cache, "CMAKE_CUDA_ARCHITECTURES"),
        "cuda_compiler": cache_value(cache, "CMAKE_CUDA_COMPILER"),
        "cxx_compiler": cache_value(cache, "CMAKE_CXX_COMPILER"),
        "fortran_compiler": cache_value(cache, "CMAKE_Fortran_COMPILER"),
        "cuda_flags_release": cache_value(cache, "CMAKE_CUDA_FLAGS_RELEASE"),
        "cxx_flags_release": cache_value(cache, "CMAKE_CXX_FLAGS_RELEASE"),
        "fortran_flags_release": cache_value(cache, "CMAKE_Fortran_FLAGS_RELEASE"),
    }


def host_metadata() -> dict:
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except AttributeError:  # pragma: no cover - non-Linux
        affinity = []
    model = ""
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(errors="replace").splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {
        "platform": platform.platform(),
        "cpu_model": model,
        "cpu_affinity": affinity,
        "cpu_affinity_count": len(affinity),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "unset"),
        "omp_proc_bind": os.environ.get("OMP_PROC_BIND", "unset"),
        "omp_places": os.environ.get("OMP_PLACES", "unset"),
    }


ACTIVE_RE = re.compile(
    r"AdaptiveCG initial active_atoms=(\d+) active_blocks=(\d+) "
    r"interface_atoms=(\d+) device_bytes=(\d+)")
PHASE_TIMING_RE = re.compile(r"AdaptiveCG per-phase device timing (\S+)")


def run_case(binary: Path, case: Path, nstep: int, workspace: Path,
             phase_timing: bool, device: int, timeout: float) -> dict:
    """Run one case at one Nstep in a throwaway copy of its directory."""
    # The whole perf tree is copied, not just the case: the cases share
    # ../posfile, ../momfile and ../jfile so that both arms provably read the
    # same geometry and exchange from one tracked source.
    root = workspace / f"{case.name}_n{nstep}_{'timed' if phase_timing else 'lean'}"
    if root.exists():
        shutil.rmtree(root)
    shutil.copytree(PERF_ROOT, root)
    target = root / case.name
    inp = target / "inpsd.dat"
    text = inp.read_text()
    if not re.search(r"^Nstep\s+\d+\s*$", text, flags=re.MULTILINE):
        raise RuntimeError(f"{inp} has no Nstep line to rewrite")
    inp.write_text(re.sub(r"^Nstep\s+\d+\s*$", f"Nstep {nstep}",
                          text, flags=re.MULTILINE))

    environment = dict(os.environ)
    environment["CUDA_VISIBLE_DEVICES"] = str(device)
    if not phase_timing:
        environment["UPPASD_ADAPTIVE_PHASE_TIMING"] = "0"
    else:
        environment.pop("UPPASD_ADAPTIVE_PHASE_TIMING", None)

    started = time.perf_counter()
    result = subprocess.run(
        [str(binary)],
        cwd=target,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        check=False,
    )
    wall = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(
            f"{case.name} at Nstep={nstep} exited {result.returncode}\n"
            f"{result.stdout[-4000:]}")
    active = ACTIVE_RE.search(result.stdout)
    timing_state = PHASE_TIMING_RE.search(result.stdout)
    return {
        "case": case.name,
        "nstep": nstep,
        "wall_s": wall,
        "phase_timing_requested": phase_timing,
        "phase_timing_reported": timing_state.group(1) if timing_state else "n/a",
        "active_atoms": int(active.group(1)) if active else None,
        "active_blocks": int(active.group(2)) if active else None,
        "interface_atoms": int(active.group(3)) if active else None,
        "device_bytes": int(active.group(4)) if active else None,
        "stdout_tail": result.stdout[-2000:],
    }


def per_step_us(short: dict, long: dict) -> float:
    steps = long["nstep"] - short["nstep"]
    if steps <= 0:
        raise ValueError("long run must have more steps than short run")
    return 1.0e6 * (long["wall_s"] - short["wall_s"]) / steps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path, help="UppASD GPU executable")
    parser.add_argument("--texture", action="append", choices=sorted(TEXTURES),
                        help="texture to measure (default: all)")
    parser.add_argument("--rounds", type=int, default=5,
                        help="interleaved A/B rounds (default 5)")
    parser.add_argument("--short-steps", type=int, default=20)
    parser.add_argument("--long-steps", type=int, default=220)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--timeout", type=float, default=1800.0)
    parser.add_argument("--crossover-margin-percent", type=float, default=2.0)
    parser.add_argument("--max-utilization-percent", type=float, default=5.0)
    parser.add_argument("--max-temperature-c", type=float, default=60.0)
    parser.add_argument("--allow-dirty-device", action="store_true",
                        help="smoke-test escape hatch; results are NOT acceptance evidence")
    parser.add_argument("--json", type=Path, help="write raw samples here")
    arguments = parser.parse_args()

    binary = arguments.binary.resolve()
    if not binary.is_file():
        print(f"missing binary {binary}", file=sys.stderr)
        return 1
    textures = arguments.texture or sorted(TEXTURES)
    margin = arguments.crossover_margin_percent / 100.0

    build = build_metadata(binary)
    host = host_metadata()
    print("rcg09-e2e-benchmark " + " ".join(
        f"{key}={value!r}" for key, value in build.items()), flush=True)
    print("rcg09-e2e-host " + " ".join(
        f"{key}={value!r}" for key, value in host.items()), flush=True)
    device = check_device_clean(arguments.device,
                               arguments.max_utilization_percent,
                               arguments.max_temperature_c,
                               arguments.allow_dirty_device)

    record: dict = {
        "build": build,
        "host": host,
        "device_before": device,
        "short_steps": arguments.short_steps,
        "long_steps": arguments.long_steps,
        "rounds": arguments.rounds,
        "crossover_margin": margin,
        "textures": {},
    }
    exit_code = 0

    with tempfile.TemporaryDirectory(prefix="rcg09-perf-") as scratch:
        workspace = Path(scratch)
        for texture in textures:
            atomistic_name, adaptive_name = TEXTURES[texture]
            atomistic_case = PERF_ROOT / atomistic_name
            adaptive_case = PERF_ROOT / adaptive_name
            for case in (atomistic_case, adaptive_case):
                if not (case / "inpsd.dat").is_file():
                    raise SystemExit(f"missing fixture {case}")

            atomistic_samples: list[float] = []
            adaptive_samples: list[float] = []
            runs: list[dict] = []
            # Interleaved and order-alternating: adjacent samples of the two
            # arms see the same clock and thermal state, so the ratio survives
            # drift even though the absolute numbers would not.
            for round_index in range(arguments.rounds):
                order = [("atomistic", atomistic_case), ("adaptive", adaptive_case)]
                if round_index % 2:
                    order.reverse()
                for arm, case in order:
                    short = run_case(binary, case, arguments.short_steps, workspace,
                                     False, arguments.device, arguments.timeout)
                    long = run_case(binary, case, arguments.long_steps, workspace,
                                    False, arguments.device, arguments.timeout)
                    step_us = per_step_us(short, long)
                    runs.extend([short, long])
                    (atomistic_samples if arm == "atomistic"
                     else adaptive_samples).append(step_us)
                    print(f"rcg09-e2e-sample texture={texture} arm={arm} "
                          f"round={round_index} short_wall_s={short['wall_s']:.6f} "
                          f"long_wall_s={long['wall_s']:.6f} "
                          f"step_us={step_us:.6f} "
                          f"phase_timing={short['phase_timing_reported']}",
                          flush=True)

            # One instrumented adaptive pass: the per-phase breakdown, and the
            # cost of obtaining it (RCG-08-FU2).
            instrumented_short = run_case(binary, adaptive_case, arguments.short_steps,
                                          workspace, True, arguments.device,
                                          arguments.timeout)
            instrumented_long = run_case(binary, adaptive_case, arguments.long_steps,
                                         workspace, True, arguments.device,
                                         arguments.timeout)
            instrumented_us = per_step_us(instrumented_short, instrumented_long)
            runs.extend([instrumented_short, instrumented_long])

            atomistic_median = median(atomistic_samples)
            adaptive_median = median(adaptive_samples)
            atomistic_mad = median_absolute_deviation(atomistic_samples)
            adaptive_mad = median_absolute_deviation(adaptive_samples)
            uncertainty = 3.0 * (atomistic_mad + adaptive_mad)
            target = atomistic_median * (1.0 - margin)
            crossed = adaptive_median + uncertainty < target
            speedup = (atomistic_median / adaptive_median) if adaptive_median else 0.0

            print(f"rcg09-e2e-result texture={texture} "
                  f"atomistic_step_us={atomistic_median:.6f} "
                  f"atomistic_mad_us={atomistic_mad:.6f} "
                  f"adaptive_step_us={adaptive_median:.6f} "
                  f"adaptive_mad_us={adaptive_mad:.6f} "
                  f"speedup={speedup:.6f} uncertainty_us={uncertainty:.6f} "
                  f"margin_percent={arguments.crossover_margin_percent:.3f} "
                  f"instrumented_adaptive_step_us={instrumented_us:.6f} "
                  f"instrumentation_overhead_percent="
                  f"{100.0 * (instrumented_us - adaptive_median) / adaptive_median:.6f} "
                  f"active_atoms={instrumented_short['active_atoms']} "
                  f"active_blocks={instrumented_short['active_blocks']} "
                  f"interface_atoms={instrumented_short['interface_atoms']} "
                  f"adaptive_device_bytes={instrumented_short['device_bytes']} "
                  f"crossover={'PASS' if crossed else 'NOT_OBSERVED'}", flush=True)
            if not crossed:
                exit_code = max(exit_code, 3)

            record["textures"][texture] = {
                "atomistic_samples_us": atomistic_samples,
                "adaptive_samples_us": adaptive_samples,
                "atomistic_median_us": atomistic_median,
                "adaptive_median_us": adaptive_median,
                "atomistic_mad_us": atomistic_mad,
                "adaptive_mad_us": adaptive_mad,
                "uncertainty_us": uncertainty,
                "speedup": speedup,
                "crossover": bool(crossed),
                "instrumented_adaptive_step_us": instrumented_us,
                "runs": runs,
            }

    record["device_after"] = device_state(arguments.device)
    if arguments.json:
        arguments.json.write_text(json.dumps(record, indent=2, sort_keys=True))
        print(f"rcg09-e2e-artifact path={arguments.json}", flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
