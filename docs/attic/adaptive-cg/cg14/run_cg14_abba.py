#!/usr/bin/env python3
"""CG-14 ABBA-interleaved measurement.

Protocol follows RCG-09C section 2: at each size the first invocation is a
discarded warm-up, then the baseline and live arms alternate in ABBA order,
`--reps-per-arm` runs per arm, medians reported.  The production-oracle row is
identical code in both binaries and is the check that both arms saw the same
device.  nvidia-smi state is sampled before and after every size.
"""
import argparse, json, re, statistics, subprocess, sys, time

BASE = "/tmp/cg14_baseline/build_cg14_base_cuda/bin/gpu_adaptive_runtime_benchmark"
LIVE = "/home/andersb/SD/UppASD_gpu_hip_cu/build_cg14_cuda/bin/gpu_adaptive_runtime_benchmark"

SWEEP = re.compile(r"adaptive-sweep requested_fraction=([\d.]+) .*?active_blocks=(\d+) "
                   r".*?step_wall_us=([\d.]+)")
PHASE = re.compile(r"adaptive-sweep-phases requested_fraction=([\d.]+) atomistic_us=([\d.]+) "
                   r"coarse_us=([\d.]+) interface_us=([\d.]+) .*?integration_us=([\d.]+)")
PROD = re.compile(r"production-atomistic .*?step_wall_us=([\d.]+)")
CROSS_NO = re.compile(r"production-crossover result=(NOT_OBSERVED) "
                      r"production_step_wall_us=([\d.]+) "
                      r"best_adaptive_step_wall_us=([\d.]+)")
# The benchmark prints a different line when its conservative crossover test
# passes: fraction, the adaptive step at that point, production, and speedup.
CROSS_YES = re.compile(r"production-crossover result=(PASS) requested_fraction=([\d.]+) "
                       r".*?adaptive_step_wall_us=([\d.]+) "
                       r"production_step_wall_us=([\d.]+) speedup=([\d.]+)")
SELFREF = re.compile(r"adaptive-self-reference .*?zero_fine_step_wall_us=([\d.]+)")


def smi():
    out = subprocess.run(
        ["nvidia-smi", "--id=0", "--query-gpu=utilization.gpu,memory.used,temperature.gpu,"
         "clocks.sm,clocks_throttle_reasons.active", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    procs = subprocess.run(
        ["nvidia-smi", "--id=0", "--query-compute-apps=pid,used_memory", "--format=csv,noheader"],
        capture_output=True, text=True).stdout.strip()
    return {"gpu": out, "compute_apps": procs or "none"}


def run(binary, blocks, apb, texture):
    cmd = [binary, "--blocks", str(blocks), "--atoms-per-block", str(apb),
           "--warmup", "1", "--iterations", "3", "--repetitions", "3",
           "--texture", texture]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True,
                       env={"CUDA_VISIBLE_DEVICES": "0", "PATH": "/usr/bin:/bin"})
    if p.returncode != 0:
        raise SystemExit(f"benchmark failed: {' '.join(cmd)}\n{p.stdout[-3000:]}{p.stderr[-2000:]}")
    o = p.stdout
    points = {}
    for m in SWEEP.finditer(o):
        points.setdefault(m.group(1), {})["step_wall_us"] = float(m.group(3))
        points[m.group(1)]["active_blocks"] = int(m.group(2))
    for m in PHASE.finditer(o):
        d = points.setdefault(m.group(1), {})
        d["atomistic_us"] = float(m.group(2)); d["coarse_us"] = float(m.group(3))
        d["interface_us"] = float(m.group(4)); d["integration_us"] = float(m.group(5))
    prod = [float(x) for x in PROD.findall(o)]
    no, yes = CROSS_NO.search(o), CROSS_YES.search(o)
    sr = SELFREF.search(o)
    out = {"points": points, "production_us": prod,
           "zero_fine_step_wall_us": float(sr.group(1)) if sr else None,
           "seconds": round(time.time() - t0, 1)}
    if yes:
        out.update(crossover="PASS",
                   crossover_fraction=float(yes.group(2)),
                   crossover_adaptive_step_wall_us=float(yes.group(3)),
                   production_step_wall_us=float(yes.group(4)),
                   crossover_speedup=float(yes.group(5)),
                   best_adaptive_step_wall_us=float(yes.group(3)))
    elif no:
        out.update(crossover="NOT_OBSERVED",
                   production_step_wall_us=float(no.group(2)),
                   best_adaptive_step_wall_us=float(no.group(3)))
    else:
        raise SystemExit("no production-crossover line in benchmark output")
    return out


def median_of(runs, fraction, key):
    vals = [r["points"][fraction][key] for r in runs
            if fraction in r["points"] and key in r["points"][fraction]]
    return statistics.median(vals) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", required=True, help="blocks:atomsperblock,...")
    ap.add_argument("--reps-per-arm", type=int, default=4)
    ap.add_argument("--texture", default="spiral")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    results = {"texture": a.texture, "reps_per_arm": a.reps_per_arm, "sizes": []}
    for spec in a.sizes.split(","):
        blocks, apb = (int(x) for x in spec.split(":"))
        entry = {"blocks": blocks, "atoms_per_block": apb,
                 "device_before": smi(), "runs": {"base": [], "live": []}}
        print(f"[{spec}] warm-up", flush=True)
        run(LIVE, blocks, apb, a.texture)                       # discarded
        order = []
        for _ in range(a.reps_per_arm // 2):
            order += ["base", "live", "live", "base"]           # ABBA
        if a.reps_per_arm % 2:
            order += ["base", "live"]
        for i, arm in enumerate(order):
            r = run(BASE if arm == "base" else LIVE, blocks, apb, a.texture)
            entry["runs"][arm].append(r)
            print(f"[{spec}] {i+1}/{len(order)} {arm:4s} {r['seconds']}s "
                  f"prod={r['production_step_wall_us']:.1f} "
                  f"best={r['best_adaptive_step_wall_us']:.1f}", flush=True)
        entry["device_after"] = smi()

        fractions = sorted({f for arm in entry["runs"] for r in entry["runs"][arm]
                            for f in r["points"]}, key=float, reverse=True)
        entry["table"] = []
        for f in fractions:
            row = {"fraction": float(f),
                   "active_blocks": entry["runs"]["live"][0]["points"][f]["active_blocks"]}
            for key in ("atomistic_us", "coarse_us", "interface_us",
                        "integration_us", "step_wall_us"):
                row[f"base_{key}"] = median_of(entry["runs"]["base"], f, key)
                row[f"live_{key}"] = median_of(entry["runs"]["live"], f, key)
            entry["table"].append(row)
        for arm in ("base", "live"):
            entry[f"{arm}_production_step_wall_us"] = statistics.median(
                [r["production_step_wall_us"] for r in entry["runs"][arm]])
            entry[f"{arm}_best_adaptive_step_wall_us"] = statistics.median(
                [r["best_adaptive_step_wall_us"] for r in entry["runs"][arm]])
            entry[f"{arm}_crossover"] = sorted({r["crossover"] for r in entry["runs"][arm]})
            entry[f"{arm}_zero_fine_step_wall_us"] = statistics.median(
                [r["zero_fine_step_wall_us"] for r in entry["runs"][arm]])
            sp = [r.get("crossover_speedup") for r in entry["runs"][arm]
                  if r.get("crossover_speedup")]
            entry[f"{arm}_crossover_speedup"] = statistics.median(sp) if sp else None
            fr = [r.get("crossover_fraction") for r in entry["runs"][arm]
                  if r.get("crossover_fraction") is not None]
            entry[f"{arm}_crossover_fraction"] = sorted(set(fr)) if fr else None
        results["sizes"].append(entry)
        json.dump(results, open(a.out, "w"), indent=1)
    json.dump(results, open(a.out, "w"), indent=1)
    print("written", a.out)


if __name__ == "__main__":
    main()
