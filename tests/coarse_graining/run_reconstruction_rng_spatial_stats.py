#!/usr/bin/env python3
"""RCG-06D: reconstruction RNG spatial-statistics evidence (F-20).

F-20 is "the tuple-seeded MINSTD reconstruction can correlate nearby
refinements." The parent blueprint's instruction is to measure this, not
assume it, and replace the generator only if the *accepted* nonzero-cone
model needs stronger independence -- otherwise document why zero-cone
production is unaffected and leave nonzero-cone acceptance open. No tracked
fixture uses ``cg_reconstruction CONSTRAINED_CONE`` with a nonzero cone
angle (every one uses ``ALIGNED``, the production default), so there is no
"accepted nonzero-cone model" today; this fixture supplies the measurement
that a future acceptance decision would need.

This driver runs ``test_reconstruction_rng_spatial_stats.f90`` (the CPU
binary, always) and, if supplied, ``test_reconstruction_rng_gpu_parity.cpp``
(the GPU binary), both of which emit raw draws from the literal production
generator (``deterministic_reconstruction_seed`` + ``next_uniform`` on CPU;
the identical formula via the shared ``gpuAdaptiveReconstructionRng.hpp`` on
GPU) over a sweep of block indices at several fixed
(global_seed,channel,ensemble,epoch) tuples.

Two independent things are checked:

1. **Spatial correlation is real and stable, not measurement noise.**  For
   each fixed (global_seed,channel,ensemble,epoch) tuple, the first draw's
   value as a function of block index is an affine (mod-prime) sequence by
   construction -- adjacent block seeds differ by a fixed additive offset
   before being run through one multiplicative LCG step, so their first
   outputs differ by a fixed additive offset too (mod the prime modulus).
   This produces a genuine, deterministic, non-decaying lag correlation
   between nearby blocks' first draws -- measured here, not asserted a
   priori. A block-shuffled negative control (same draws, block identity
   discarded) is required to show near-zero correlation, proving the
   fixture is discriminating structure from an artifact of the estimator
   itself.
2. **CPU/GPU parity**, if a GPU binary is supplied: before this slice, the
   GPU device generator used multiplier 16807 (the original Lehmer/"MINSTD"
   constant) while the CPU generator used 48271 (the revised Park-Miller
   "minimal standard" constant) against the same modulus and seed formula,
   so identical tuple seeds produced two unrelated draw sequences. This
   slice fixed GPU to 48271 unconditionally (see
   source/gpu_files/gpuAdaptiveReconstructionRng.hpp); this check is the
   regression test for that fix.

This fixture deliberately does NOT assert an upper bound on the sequential
spatial correlation's magnitude, and does NOT replace the generator: F-20's
own instruction is to measure and record a scope decision, not to fix
nonzero-cone independence pre-emptively while it remains unaccepted for
production. The zero-cone-is-unaffected claim is proven separately, by
direct construction, in
tests/coarse_graining/test_adaptive_hybrid_solver.f90 (two different
global_seed values produce byte-identical cone_angle_rad=0 output) -- not
repeated here.
"""

from __future__ import annotations

import argparse
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

Key = Tuple[int, int, int, int]  # global_seed, channel, ensemble, epoch
GroupMap = Dict[Key, List[Tuple[int, float, float]]]


class SpatialStatsError(AssertionError):
    pass


def _run(binary: Path) -> str:
    result = subprocess.run(
        [str(binary)], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        timeout=120, check=False,
    )
    if result.returncode == 77:
        return ""
    if result.returncode != 0:
        raise SpatialStatsError(
            f"{binary} exited {result.returncode}:\n{result.stdout}"
        )
    return result.stdout


def _parse(stdout: str) -> GroupMap:
    """Group rows by (global_seed,channel,ensemble,epoch) -> [(block,u1,u2), ...]."""
    groups: GroupMap = defaultdict(list)
    for line in stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) != 7:
            continue
        global_seed, channel, ensemble, epoch, block = (int(x) for x in parts[:5])
        u1, u2 = float(parts[5]), float(parts[6])
        groups[(global_seed, channel, ensemble, epoch)].append((block, u1, u2))
    if not groups:
        raise SpatialStatsError("no data rows parsed from binary output")
    for key in groups:
        groups[key].sort(key=lambda row: row[0])
    return groups


def _pearson(x: List[float], y: List[float]) -> float:
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    vx = sum((xi - mx) ** 2 for xi in x)
    vy = sum((yi - my) ** 2 for yi in y)
    denom = (vx * vy) ** 0.5
    if denom <= 0.0:
        return 0.0
    return cov / denom


def _lag_correlation(values: List[float], lag: int) -> float:
    return _pearson(values[:-lag], values[lag:])


def analyze_spatial_correlation(groups: GroupMap) -> None:
    lags = (1, 2, 3, 4)
    sequential_r: Dict[int, List[float]] = {lag: [] for lag in lags}
    shuffled_r: List[float] = []
    rng = random.Random(20260812)

    for key, rows in groups.items():
        u1 = [row[1] for row in rows]
        for lag in lags:
            sequential_r[lag].append(_lag_correlation(u1, lag))
        shuffled = list(u1)
        rng.shuffle(shuffled)
        shuffled_r.append(_lag_correlation(shuffled, 1))

    print(f"RECONSTRUCTION-RNG-SPATIAL-STATS: {len(groups)} (global_seed,channel,ensemble,epoch) "
          f"groups, {len(next(iter(groups.values())))} blocks each")
    print(f"{'lag':>4} {'mean_r':>10} {'min_r':>10} {'max_r':>10}")
    for lag in lags:
        values = sequential_r[lag]
        mean_r = sum(values) / len(values)
        print(f"{lag:>4} {mean_r:>10.4f} {min(values):>10.4f} {max(values):>10.4f}")

    mean_lag1 = sum(sequential_r[1]) / len(sequential_r[1])
    mean_shuffled = sum(shuffled_r) / len(shuffled_r)
    print(f"block-shuffled negative control, lag-1: mean_r={mean_shuffled:.4f} "
          f"(expected near zero if the fixture is discriminating structure, "
          f"not manufacturing it)")

    # Discriminating-test check: real structure must be measured, not a flat
    # near-zero curve indistinguishable from noise (which would mean this
    # fixture is not actually exercising F-20's claim). The magnitude and
    # sign are not asserted beyond "clearly nonzero" -- the fixed structural
    # threshold below is derived from this measurement, not assumed a
    # priori; see the evidence doc for the observed numbers and analytic
    # explanation (affine-in-block-index seed construction).
    if not (abs(mean_lag1) > 0.15):
        raise SpatialStatsError(
            f"sequential lag-1 correlation mean |{mean_lag1:.4f}| is not "
            "clearly nonzero -- fixture is not discriminating, cannot "
            "support the F-20 spatial-correlation finding"
        )

    # Negative control: the shuffled sequence must show only sampling noise,
    # not the same structure -- proving the sequential result above is a
    # property of block adjacency, not an artifact of the Pearson estimator
    # applied to this generator's numbers in general.
    # Independence's 95% CI half-width is ~1.96/sqrt(n); this budget gives
    # ~2.5x that margin so ordinary sampling noise never trips the check.
    noise_budget = 5.0 / (len(next(iter(groups.values()))) ** 0.5)
    if abs(mean_shuffled) > noise_budget:
        raise SpatialStatsError(
            f"block-shuffled negative control mean_r={mean_shuffled:.4f} "
            f"exceeds the sampling-noise budget {noise_budget:.4f} -- "
            "the estimator itself may be biased, not just measuring "
            "block-adjacency structure"
        )

    print(f"RECONSTRUCTION-RNG-SPATIAL-STATS: PASS -- sequential lag-1 correlation "
          f"(mean={mean_lag1:.4f}) is real, stable across "
          f"(global_seed,channel,ensemble,epoch) tuples, and clearly distinct from the "
          f"block-shuffled negative control (mean={mean_shuffled:.4f}); scope decision "
          "recorded in docs/RCG-06_MEMORY_TIMING_PRECISION_EVIDENCE.md (RCG-06D)")


def _flatten(groups: GroupMap) -> Dict[Tuple[Key, int], Tuple[float, float]]:
    flat: Dict[Tuple[Key, int], Tuple[float, float]] = {}
    for key, rows in groups.items():
        for block, u1, u2 in rows:
            flat[(key, block)] = (u1, u2)
    return flat


def check_gpu_parity(cpu_groups: GroupMap, gpu_groups: GroupMap) -> None:
    cpu_rows = _flatten(cpu_groups)
    gpu_rows = _flatten(gpu_groups)
    missing = set(cpu_rows) - set(gpu_rows)
    if missing:
        raise SpatialStatsError(
            f"GPU output is missing {len(missing)} rows the CPU driver produced, "
            f"e.g. {sorted(missing)[:3]}"
        )
    # Loose enough to tolerate FP32 rounding in a SINGLE_PREC build, tight
    # enough that a genuine multiplier mismatch (unrelated draws, O(1)
    # differences) cannot pass silently.
    tolerance = 1.0e-5
    worst = 0.0
    worst_key = None
    for key, (cpu_u1, cpu_u2) in cpu_rows.items():
        gpu_u1, gpu_u2 = gpu_rows[key]
        diff = max(abs(cpu_u1 - gpu_u1), abs(cpu_u2 - gpu_u2))
        if diff > worst:
            worst = diff
            worst_key = key
    if worst > tolerance:
        raise SpatialStatsError(
            f"CPU/GPU reconstruction RNG parity failed: worst |diff|={worst:.6e} at "
            f"(global_seed,channel,ensemble,epoch,block)={worst_key}, exceeds tolerance "
            f"{tolerance:.1e} -- generators have diverged (multiplier mismatch?)"
        )
    print(f"CPU/GPU reconstruction RNG parity: {len(cpu_rows)} rows compared, "
          f"worst |diff|={worst:.3e} (tolerance {tolerance:.1e}) -- PASS")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cpu-binary", required=True, type=Path)
    parser.add_argument("--gpu-binary", type=Path)
    args = parser.parse_args()

    cpu_stdout = _run(args.cpu_binary)
    cpu_groups = _parse(cpu_stdout)
    analyze_spatial_correlation(cpu_groups)

    if args.gpu_binary:
        gpu_stdout = _run(args.gpu_binary)
        if not gpu_stdout:
            print("GPU binary reported no device (returncode 77); GPU parity check skipped")
        else:
            gpu_groups = _parse(gpu_stdout)
            check_gpu_parity(cpu_groups, gpu_groups)
    else:
        print("GPU binary not supplied; CPU/GPU parity check not run (deferred, not blocking)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
