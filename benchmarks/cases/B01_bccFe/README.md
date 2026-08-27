# B01_bccFe — bcc Fe, medium-range interactions

Admission record for WP-08A (blueprint §3 "B01 — bcc Fe: medium-range
interactions"). This case is the framework's central 3D production
reference: primary OpenMP scaling reference, primary CPU/GPU crossover
reference, finite-temperature Depondt/Langevin benchmark.

## A. Template audit

Source: `examples/SimpleSystems/bccFeCPU/{inpsd.dat,jfile,posfile,momfile,kfile}`,
copied byte-for-byte into `template/` (checksums verified equal; the
examples-tree copy is untouched). `examples/SimpleSystems/bccFeGPU` is the
same deck with two differences, both handled below (`do_proj_avrg`, an
output option not needed here, and `gpu_mode`, now a per-run override —
see "Backend dispatch").

| Property | Value | Source |
| --- | --- | --- |
| Lattice | bcc, conventional two-atom cubic-cell basis: `(0,0,0)` and `(0.5,0.5,0.5)` | `posfile` |
| Basis/moment magnitude | `2.23 μ_B` on both basis atoms, initial direction `(1,0,0)` (`Initmag 3`) | `momfile`, `inpsd.dat` |
| Boundary conditions | Fully periodic, `BC P P P` | `inpsd.dat` |
| Symmetry | Cubic, `Sym 1` (jfile's single declared pair `1 1` is expanded to both basis atoms and the full point group) | `inpsd.dat`, `jfile` |
| Exchange file | `jfile`, 8 shells, cutoff radius 2.0 (reduced units) | `jfile` |
| Directed interactions / atom | **96**, measured directly (see §D) — not estimated | `struct.<simid>.out` |
| Anisotropy | **Active.** `anisotropy ./kfile` sets `do_anisotropy=1` (`source/Input/inputhandler.f90` case `'anisotropy'`). `anisotropytype=7` on both basis atoms = combined uniaxial+cubic (`source/Hamiltonian/energy.f90:299-306`): K1=-0.02, K2=0.0, easy axis `(0,1,0)`, cubic-term weight (`sb`/ratio_k) 0.1 | `kfile` |
| External field | Off — `hfield 0.00000 0.00000 0.00000` | `inpsd.dat` |
| Dipole | Off — no `dipole`/FFT keyword present | `inpsd.dat` |
| Temperature (maintainer default) | **100 K** (`Temp 100.00`) — used verbatim by `bcc_fe_finite_t`, not invented | `inpsd.dat` |
| Timestep | `1.000e-17 s` | `inpsd.dat` |
| Integrator | Depondt, `SDEalgh 5` | `inpsd.dat` |
| Damping | `0.50` | `inpsd.dat` |
| Ensembles | `Mensemble 10` (ten independent replicas per run; fixed by the template, not in the override allow-list) | `inpsd.dat` |
| Initial phase | `ip_mode N` — `'N'` matches none of the branches in `source/uppasd.f90::run_initial_phase`, so **no initial phase executes**; the `ip_mcanneal` block present in the file is dead input for this template | `inpsd.dat`, `source/uppasd.f90:176-235` |
| `do_prnstruct` | `1` (already set in the maintainer template — required for the workload-metadata parser below) | `inpsd.dat` |

### Backend dispatch (`gpu_mode`)

The maintainer's CPU and GPU decks (`bccFeCPU`/`bccFeGPU`) are otherwise
identical and differ only in `gpu_mode` (0 vs 1). Tracing the dispatch
(`source/Input/inputhandler.f90`'s `gpu_mode` case sets `do_gpu`;
`source/uppasd.f90::run_initial_phase`/`run_measurement_phase` call
`sd_iphaseGPU`/`sd_mphaseGPU` only when `do_gpu=='Y'`, else the ordinary
Fortran CPU path) shows this can't be baked into one fixed template value:
on a build without CUDA/HIP compiled in, the GPU driver resolves to
`source/Tools/nocuda.f90`'s stub module, whose subroutines silently
`return` with no work done — `gpu_mode 1` baked into a CPU-only run would
silently skip the simulation instead of failing loudly. Since this case's
`case_id`/`variant_id`/`size_id` identity is shared across backend records
(TERMINOLOGY.md: RUN CONFIGURATION, not CASE, carries `backend`), `gpu_mode`
is now a per-run override (`harness/cases.py::GLOBALLY_SAFE_OVERRIDE_KEYS`,
extended by this admission) instead of template content: the template's own
baseline stays `0` (safe on every build), and a GPU sample passes
`extra_overrides={"gpu_mode": 1}` to `cases.generate_run_directory`/
`runner.run_sample`. This is a backend-dispatch switch, not a
Hamiltonian/lattice/moment parameter, so it does not touch model physics.

## B. Size ladder

`Nx=Ny=Nz=n` isotropic bulk cubes (two-atom basis, `natom = 2n^3`), chosen
as the closest legal integer replication to each target in the blueprint's
~4x ladder (4k/16k/64k/256k/1M/4M). Geometry is never distorted to force an
exact atom count.

| `size_id` | `n` | `natom` | vs. target | Note |
| --- | --- | --- | --- | --- |
| `13x13x13` | 13 | 4,394 | target 4,096, +7.3% | closest integer cube; 12 undershoots by 15.6% |
| `20x20x20` | 20 | 16,000 | target 16,384, -2.3% | |
| `32x32x32` | 32 | 65,536 | target 65,536, **exact** | this project's established historical reference size (`docs/CGP-07_ARCHITECTURE_REBASELINE_EVIDENCE.md`, `tests/gpu_regression/bench.py`'s `bcc_medium`) |
| `51x51x51` | 51 | 265,302 | target 262,144, +1.2% | |
| `81x81x81` | 81 | 1,062,882 | target 1,048,576, +1.4% | |
| `128x128x128` | 128 | 4,194,304 | target 4,194,304, **exact** | |

## C. Variants

| `variant_id` | `temp` override | `tseed` override | Purpose |
| --- | --- | --- | --- |
| `bcc_fe_t0` | 0.0 | 1 | Deterministic zero-temperature dynamics; required by `harness.gpu_sanity` (`temperature==0`) for CPU/GPU consistency checks |
| `bcc_fe_finite_t` | 100.0 (maintainer's own baseline value, not invented) | 2 | Finite-temperature Langevin production dynamics |

## D. Scaling validation

Ran the real executable (`do_prnstruct 1`, `neighbor_list_from_struct_output`)
at four sizes spanning a 60x atom-count range and read `struct.<simid>.out`
directly — never estimated:

| `size_id` | `natom` | `directed_interactions` | `directed_interactions / natom` |
| --- | --- | --- | --- |
| `13x13x13` | 4,394 | 421,824 | 96.0 |
| `20x20x20` | 16,000 | 1,536,000 | 96.0 |
| `32x32x32` | 65,536 | 6,291,456 | 96.0 |
| `51x51x51` | 265,302 | 25,468,992 | 96.0 |

`mean_neighbors == max_neighbors == 96` at every tested size: full periodic
boundaries mean every atom is bulk (no surface truncation), so the
per-atom interaction count is exactly constant across the ladder — the
same Hamiltonian, same 8-shell exchange list, same cutoff, only the atom
count changes. 96 is exactly the sum of the cubic-symmetry-expanded
coordination numbers of the 8 shells in `jfile` (8+6+12+24+6+24+8+8=96,
including the four numerically-negligible-but-present outer shells) —
confirms the parser is counting the true production topology, not an
approximation. (An earlier synthetic schema-fixture,
`benchmarks/tests/fixtures/valid/cpu_bccfe_dynamics_only.json`, assumed 58
neighbours for illustration only; it was never a measurement of this
template and this admission supersedes it as the real value.)

`81x81x81` and `128x128x128` were not executed (impractical for an
admission-time sanity check); at the confirmed constant-96 topology,
`directed_interactions` at those sizes is `natom * 96` =
102,036,672 and 402,653,184 respectively — a prediction from the
confirmed scaling law, not an independent measurement.

## E. Sanity runs

All four runs used `13x13x13` (4,394 atoms), `do_prnstruct 1`, `Nstep 200`;
CPU on `build_cpu/bin/sd.f95` (`UPPASD_GPU_BACKEND=OFF`, `UPPASD_PRECISION=DOUBLE`),
GPU on `build_gpu/bin/sd.f95.cuda` (`UPPASD_GPU_BACKEND=CUDA`,
`UPPASD_PRECISION=DOUBLE`, `extra_overrides={"gpu_mode": 1}`).

| Run | Result |
| --- | --- |
| CPU, `bcc_fe_t0` | **PASS** — completes, `natom`/`directed_interactions` match §D, energy finite and constant |
| CUDA, `bcc_fe_t0` | **PASS** — completes, identical `natom`/`directed_interactions`, identical final moments and energy to the CPU run |
| CPU, `bcc_fe_finite_t` | **PASS** — completes, `|m|` conserved at 2.23 on every atom, energy relaxes from -12.410 toward -12.219 over 100 steps |
| CUDA, `bcc_fe_finite_t` | **PASS** — completes, `|m|` conserved at 2.23 on every atom, energy relaxes from -12.410 toward -12.154 over 200 steps |

Notes:

* **The `bcc_fe_t0` comparison is a weak discriminator.** With `Initmag 3`
  giving a uniform initial state along `(1,0,0)`, both the exchange field
  (uniform-state Heisenberg exchange is always parallel to the moment,
  giving zero torque for *any* uniform direction) and the combined
  anisotropy field (zero for a state exactly aligned with a cubic axis,
  since both the uniaxial term, easy axis `(0,1,0)` orthogonal to `(1,0,0)`,
  and the cubic term vanish there) are exactly zero at `t=0` — this initial
  condition is an exact fixed point of the T=0 dynamics, confirmed by
  `Ani=0.0` and unchanged moments/energy through all 200 steps on *both*
  backends. The CPU/GPU agreement is real (same physically-null result,
  reached by two independent code paths) but does not exercise the
  per-step LLG update the way a non-degenerate initial condition would. A
  later work package wanting a stronger bit-level CPU/GPU LLG check should
  use a size/variant with a perturbed or random initial condition — not
  possible here without touching `momfile`/`Initmag`, both outside this
  case's `allowed_input_overrides` by design (§A).
* **`restart.<simid>.out` reports a different final iteration index per
  backend** at fixed T=0/finite-T sizing: CPU reports `mstep=200`
  (the real final step), GPU reports `-1` on every row. This is a
  pre-existing production behaviour (`sd_mphaseGPU`,
  `source/sd_driver.f90:1027-1112`, never updates the Fortran `mstep`
  module variable the shared restart-writer reads — the GPU loop runs
  entirely inside `gpuSim_gpuRunSimulation`), not something introduced or
  fixed by this admission. It is cosmetic: the moment/energy *values*
  written are correct and were what was compared above.
* Finite-T CPU and GPU necessarily diverge in their specific per-atom
  moment realizations (independent RNG streams per backend) — only
  magnitude conservation (`|m|=2.23` everywhere) and plausible energy
  relaxation were checked, consistent with `harness.gpu_sanity`'s own
  restriction of strict moment comparison to `temperature==0`.
* `totenergy.<simid>.out` print cadence differs slightly at the very last
  step between backends (CPU's last row is `100`, not `200`, for a
  `Nstep 200` run; GPU prints both `100` and `200`) — a pre-existing
  end-of-loop energy-print convention difference, not a validity concern
  for either run.

## Checklist

- [x] Template provenance recorded (§A; verbatim copy of `examples/SimpleSystems/bccFeCPU`, checksums verified).
- [x] Original template remains untouched (`examples/SimpleSystems/bccFeCPU` unmodified; diffed against `template/` above).
- [x] Medium interaction cloud characterized (§A/§D: 8 shells, cutoff 2.0, 96 directed interactions/atom).
- [x] 3D size ladder created (§B, six sizes 4,394 → 4,194,304 atoms).
- [x] T=0 variant created (`bcc_fe_t0`).
- [x] finite-T variant created (`bcc_fe_finite_t`, maintainer's own 100 K).
- [x] Atom/bond counts validated (§D, measured at four sizes spanning 60x, constant 96 interactions/atom).
- [x] CPU sanity passes (§E).
- [x] CUDA sanity passes (§E).
- [x] dynamics-only profile verified (`measurement_profiles.DYNAMICS_ONLY` overrides `avrg_step`/`cumu_step`/`tottraj_step`/`ene_step` past `Nstep`, all declared in `allowed_input_overrides`; exercised implicitly since every sanity run's cadence keys are within this case's allow-list).
