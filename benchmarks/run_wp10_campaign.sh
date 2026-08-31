#!/usr/bin/env bash
# WP-10 authoritative campaign launcher.
#
# Meant to be started manually during off-hours on this shared host (both
# GPUs and all 8 physical cores are held for the campaign's duration) and
# left to run unattended. Estimated total wall time on this host: roughly
# 10-11 hours (dominated by the largest size tier's fixed/setup cost; see
# benchmarks/campaigns/wp10_2026-08-28.yaml's description for how that
# estimate was derived).
#
# Usage:
#   benchmarks/run_wp10_campaign.sh
#
# Detaches from the terminal (setsid + nohup), so it is safe to close the
# terminal/SSH session after launch. Progress can be followed with:
#   tail -f benchmarks/results/wp10_2026-08-28/driver_log.jsonl
# Completion is marked by a {"kind": "campaign_done", ...} line at the end
# of that file. A plain-text launcher log also goes to
# benchmarks/results/wp10_2026-08-28_launcher.log.
#
# Safe to re-run: it wipes only this campaign's own work/results directories
# before starting (nothing outside benchmarks/work/wp10_2026-08-28 and
# benchmarks/results/wp10_2026-08-28 is touched), so an interrupted or
# aborted attempt never leaves stale/partial data mixed into a later run.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BENCH_DIR="$REPO_ROOT/benchmarks"
CAMPAIGN_ID="wp10_2026-08-28"
CAMPAIGN_MANIFEST="$BENCH_DIR/campaigns/${CAMPAIGN_ID}.yaml"
RESULTS_DIR="$BENCH_DIR/results/${CAMPAIGN_ID}"
WORK_ROOT="$BENCH_DIR/work/${CAMPAIGN_ID}"
LAUNCHER_LOG="$BENCH_DIR/results/${CAMPAIGN_ID}_launcher.log"

CPU_BUILD="$REPO_ROOT/build_wp10_cpu"
CPU_BINARY="$CPU_BUILD/bin/sd.f95"
GPU_FP32_BUILD="$REPO_ROOT/build_wp10_cuda_fp32"
GPU_FP32_BINARY="$GPU_FP32_BUILD/bin/sd.f95.cuda"
GPU_FP64_BUILD="$REPO_ROOT/build_wp10_cuda_fp64"
GPU_FP64_BINARY="$GPU_FP64_BUILD/bin/sd.f95.cuda"

MACHINE_ID="dev-host-i9-11900-a4000x2"
GPU_INDEX=0
# The repo lives on /home, a shared partition observed at 100% nominal
# capacity with only ~20GiB free (df's Avail column) -- separate from and
# much tighter than /tmp's headroom. B04_dhcpNd's workload-metadata probe
# alone can transiently write several GiB (see the campaign manifest's
# description), so this threshold is deliberately conservative and the
# preflight also prints the live disk state -- read it before trusting the
# rest of an off-hours run to unattended cleanup.
MIN_FREE_DISK_GIB=12

fail() { echo "ERROR: $*" >&2; exit 1; }

echo "== WP-10 campaign preflight =="

[[ -f "$CAMPAIGN_MANIFEST" ]] || fail "campaign manifest not found: $CAMPAIGN_MANIFEST"
[[ -x "$CPU_BINARY" ]] || fail "CPU binary not found/executable: $CPU_BINARY (build build_wp10_cpu first)"
[[ -x "$GPU_FP32_BINARY" ]] || fail "CUDA SINGLE binary not found/executable: $GPU_FP32_BINARY"
[[ -x "$GPU_FP64_BINARY" ]] || fail "CUDA DOUBLE binary not found/executable: $GPU_FP64_BINARY"

# Refuse to run these binaries against a tree whose *physics source* has
# diverged from the commit they were built from (blueprint section A:
# never mix builds from different revisions in one performance curve).
# Checks source/ and CMakeLists.txt specifically, not an exact HEAD match:
# an exact-commit check breaks every time an unrelated benchmarks/tooling
# commit lands after this one -- what actually invalidates a build is a
# change to the compiled source, not the repository's commit count moving
# on (found the hard way: this check's original exact-HEAD-match form was
# already stale by the time run_dhcpnd_convolution_campaign.sh was added).
BUILD_COMMIT="ad44f260aa918aca2304d2bf22bd544e11e199bc"
CURRENT_COMMIT="$(git -C "$REPO_ROOT" rev-parse HEAD)"
if ! git -C "$REPO_ROOT" diff --quiet "$BUILD_COMMIT" HEAD -- source/ CMakeLists.txt; then
  fail "source/ or CMakeLists.txt differs between the build commit ($BUILD_COMMIT) and HEAD ($CURRENT_COMMIT) -- the existing build_wp10_* binaries no longer match the checked-out physics source. Rebuild from HEAD, or check out $BUILD_COMMIT."
fi

# GPU courtesy check: this host's two A4000s are shared with other users.
# Refuse to start an hours-long exclusive campaign if someone else is
# already actively using either device.
GPU_BUSY="$(nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader 2>/dev/null || true)"
if [[ -n "$GPU_BUSY" ]]; then
  fail "GPU(s) already have active compute processes:
$GPU_BUSY
Refusing to start an hours-long exclusive campaign while the GPU is in use. Re-run once it is idle."
fi

df -h "$REPO_ROOT"
FREE_GIB="$(df --output=avail -BG "$REPO_ROOT" | tail -1 | tr -dc '0-9')"
if (( FREE_GIB < MIN_FREE_DISK_GIB )); then
  fail "only ${FREE_GIB}GiB free on the filesystem holding $REPO_ROOT (shown above); want at least ${MIN_FREE_DISK_GIB}GiB headroom (B04's workload-metadata probe alone can transiently write several GiB). This is a shared partition -- check what else is using it before retrying, or free some space."
fi

echo "Preflight OK: commit $CURRENT_COMMIT, GPU idle, ${FREE_GIB}GiB free disk on the partition holding $REPO_ROOT."

echo "== Wiping this campaign's own work/results directories (nothing else is touched) =="
rm -rf "$RESULTS_DIR" "$WORK_ROOT"
mkdir -p "$RESULTS_DIR" "$WORK_ROOT"

echo "== Launching (detached; safe to close this terminal) =="
cd "$BENCH_DIR"
setsid nohup python3 harness/wp10_driver.py \
  --campaign "$CAMPAIGN_MANIFEST" \
  --results-dir "$RESULTS_DIR" \
  --work-root "$WORK_ROOT" \
  --machine-id "$MACHINE_ID" \
  --cpu-binary "$CPU_BINARY" \
  --cpu-build-dir "$CPU_BUILD" \
  --gpu-fp32-binary "$GPU_FP32_BINARY" \
  --gpu-fp32-build-dir "$GPU_FP32_BUILD" \
  --gpu-fp64-binary "$GPU_FP64_BINARY" \
  --gpu-fp64-build-dir "$GPU_FP64_BUILD" \
  </dev/null >"$LAUNCHER_LOG" 2>&1 &

DRIVER_PID=$!
disown
echo "Started, PID $DRIVER_PID."
echo "Follow progress:   tail -f $RESULTS_DIR/driver_log.jsonl"
echo "Launcher stderr/out: $LAUNCHER_LOG"
echo "Stop it early:      kill $DRIVER_PID   (and any lingering sd.f95* child: pkill -f wp10_driver.py)"
echo "Done when the last line of driver_log.jsonl has \"kind\": \"campaign_done\"."
