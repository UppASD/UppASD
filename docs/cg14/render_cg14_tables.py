#!/usr/bin/env python3
"""Render the CG-14 ABBA results as the tables the evidence document needs."""
import json, sys

for path in sys.argv[1:]:
    d = json.load(open(path))
    for e in d["sizes"]:
        b, apb = e["blocks"], e["atoms_per_block"]
        print(f"\nblocks={b}  atoms_per_block={apb}  atoms={b*apb}  "
              f"reps_per_arm={d['reps_per_arm']}  texture={d['texture']}")
        print(f"  production oracle step_wall_us: baseline "
              f"{e['base_production_step_wall_us']:.1f}, live "
              f"{e['live_production_step_wall_us']:.1f}")
        print(f"  device before: {e['device_before']['gpu']}  apps={e['device_before']['compute_apps']}")
        print(f"  device after : {e['device_after']['gpu']}  apps={e['device_after']['compute_apps']}")
        print()
        print("  fine   coarse   |     coarse us      |    step wall us    | coarse   step")
        print("  frac   blocks   |  base      live    |  base      live    | speedup  speedup")
        for r in e["table"]:
            cs = r["base_coarse_us"] / r["live_coarse_us"] if r["live_coarse_us"] else float("nan")
            ss = r["base_step_wall_us"] / r["live_step_wall_us"] if r["live_step_wall_us"] else float("nan")
            print(f"  {r['fraction']:.4f} {r['active_blocks']:7d}   "
                  f"| {r['base_coarse_us']:8.1f} {r['live_coarse_us']:8.1f}  "
                  f"| {r['base_step_wall_us']:8.1f} {r['live_step_wall_us']:8.1f}  "
                  f"| {cs:6.2f}x  {ss:5.2f}x")
        allc = [r for r in e["table"] if r["fraction"] == 0.0]
        if allc:
            r = allc[0]
            nb = r["active_blocks"]
            print(f"\n  all-coarse coarse-phase constant: baseline "
                  f"{1000*r['base_coarse_us']/nb:.2f} ns/block -> live "
                  f"{1000*r['live_coarse_us']/nb:.2f} ns/block")
        print(f"  crossover: baseline {e['base_crossover']}, live {e['live_crossover']}")
        if e.get('live_crossover_speedup'):
            print(f"    live crossover at fine fractions {e['live_crossover_fraction']}, "
                  f"speedup {e['live_crossover_speedup']:.3f}x")
        if e.get('base_crossover_speedup'):
            print(f"    baseline crossover at fine fractions {e['base_crossover_fraction']}, "
                  f"speedup {e['base_crossover_speedup']:.3f}x")
        print(f"  all-coarse (zero-fine) step: baseline {e['base_zero_fine_step_wall_us']:.1f} us, "
              f"live {e['live_zero_fine_step_wall_us']:.1f} us")
        print(f"  best adaptive step: baseline {e['base_best_adaptive_step_wall_us']:.1f} us, "
              f"live {e['live_best_adaptive_step_wall_us']:.1f} us "
              f"(production {e['live_production_step_wall_us']:.1f} us)")
