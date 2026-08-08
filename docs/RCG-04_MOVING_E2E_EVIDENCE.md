# RCG-04 moving production-e2e evidence

**Status: RCG-04A only. Inventory and evidence-contract slice — no fixture,
generator, parser, or production diagnostic was added or changed.** This
document is created fresh; no prior RCG-04 evidence document existed.

**Base commit:** `e382423cec73537aad4023bcf6f7a9d78d5bc444` ("RCG-03: close,
with HIP and independent review deferred by Human decision"), which is the
exact commit named by `docs/RCG-04_MOVING_E2E_PROMPT_PACK.md` section 4 as
the required starting point. `git status --short` at session start showed no
staged or modified tracked files other than one pre-existing anomaly,
recorded honestly below; RCG-04A adds only this document.

**Pre-existing worktree anomaly (not caused by, and not fixed by, this
slice):** `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` was
already modified in the working tree before this session made any edit — a
verbatim copy of this session's own RCG-04A prompt text is spliced into the
middle of a sentence in the RCG-03 evidence section (inside the
"already-coarse unsafe blocks refine at an accepted synchronization point..."
paragraph). This was visible on the very first read of that file, before any
tool other than `Read` had been used. Per the user's explicit instruction
(2026-08-08) this file is left untouched as an unrelated worktree change; it
is the reason a fresh build's `git describe` reports a `-dirty` suffix in
this document's evidence below, unrelated to RCG-04A.

Dependencies **RCG-02** (`docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md`) and
**RCG-03** (`docs/RCG-03_POLARIZATION_ANISOTROPY_EVIDENCE.md`) are both
closed by Human decision (2026-08-08), each with HIP execution evidence and
a separate independent adversarial review explicitly deferred, not open
correctness questions. Both governing blueprints
(`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`,
`docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`) and the RCG-04 prompt pack
(`docs/RCG-04_MOVING_E2E_PROMPT_PACK.md`) were read in full before this
inventory was produced.

---

## 1. Scope of this slice

RCG-04A is an audit slice. It:

- inventories `tests/coarse_graining/run_production_e2e.py`, every directory
  it launches, the corresponding `inpsd.dat` files, `CMakeLists.txt`
  registration, and the production routines that print `AdaptiveCG:`/`Gpu:
  AdaptiveCG` diagnostics;
- classifies each current case's actual initial state and motion evidence,
  honestly, including cases whose *label* implies dynamics but whose *state*
  does not;
- maps every existing assertion to the claim it can genuinely support;
- defines the canonical moving-evidence record schema for RCG-04B onward;
- records a fresh CPU configure/build/test run and the fixture-dependency
  audit.

It does **not** implement moving-state generators, trajectory parsers,
physical oracles, new production diagnostics, or final tolerances. Those are
RCG-04B (generators), RCG-04C (parsing/observables), and RCG-04D–J
(fixtures, precision budgets, closure).

---

## 2. Harness and registration inventory

### 2.1 `run_production_e2e.py`

A single Python script (`tests/coarse_graining/run_production_e2e.py`, 434
lines) drives the real `sd.f95`/`sd.f95.cuda`/`sd.f95.hip` executable
directly against ordinary `inpsd.dat` inputs (no internal staging — matches
the CG-10.5 "must be exercised through the ordinary executable" contract).
It is registered as exactly one CTest target:

```text
add_test(NAME adaptive-cg-production-e2e
   COMMAND python3 tests/coarse_graining/run_production_e2e.py
      --binary $<TARGET_FILE:${UppASD_EXE}>
      [--gpu-binary $<TARGET_FILE:${UppASD_EXE}> if USE_CUDA/USE_HIP])
LABELS "coarse-graining;cg13;cg13-{cpu,cuda,hip};production"
ENVIRONMENT "OMP_NUM_THREADS=1"
```

`CMakeLists.txt` (~line 1013) registers this alongside the separate
`adaptive-cg-setup-rejection-matrix` (a different Python harness,
`run_setup_rejection_matrix.py`, that mutates `static_mixed`'s `inpsd.dat` to
prove setup-time capability rejections — out of RCG-04's moving-dynamics
scope, mentioned here only for completeness) and
`adaptive-cg-fixture-dependencies` (the tracked-fixture audit). There is no
per-fixture CTest registration inside `run_production_e2e.py`: all ~37
launches happen inside one Python process and one pass/fail CTest result, so
a single fixture regression fails the entire `adaptive-cg-production-e2e`
test, not a named sub-case.

The harness helper functions define what "passing" actually means today:

- `metric(output, name)` — regex-extracts one integer `name=<int>` value
  (used for `active_atom_updates`, `accepted_transitions`, etc.);
- `float_metric(output, name)` — regex-extracts the *last* matching
  `name=<float>` occurrence (used for energies, field checksums, direction
  sums/norms);
- `final_state(output)` — parses only the line printed with
  `label=final`/`final_state`, i.e. the block-resolution vector at the *last*
  step, never an intermediate step;
- `restart_state(case)` — parses `restart.<simid>.out`, which
  `source/System/restart.f90:prn_mag_conf_iter` (format
  `i8,i8,i8,2x,es16.8×4` = iter, ensemble, atom, `|M|`, `Mx`, `My`, `Mz`)
  writes **once, for the final iteration only** — it is a final-state
  snapshot file, not a trajectory, regardless of `do_tottraj`;
- `close(reference, candidate, rel, abs)` — a thin `math.isclose` wrapper.

No function in this file reads `moment.<simid>.out` (the ordinary
per-step/per-atom trajectory file written by
`source/Measurement/prn_trajectories.f90:print_trajectories` when
`do_tottraj='Y'`). Every fixture in the suite sets `do_tottraj N`
(confirmed by grep across all 38 tracked `inpsd.dat` files: zero matches for
`do_tottraj Y`). **No case in the current suite records or parses any
per-step, per-atom state.** The only time-resolved production output
consumed anywhere is the transition log printed by
`print_transition_events` (`adaptivecgproduction.f90:1448`), which the
harness never parses today (it only regex-scans for the aggregate
`accepted_transitions=`/`rejected_transitions=` counts from
`print_adaptive_cg_summary`).

### 2.2 Production diagnostic routines (`source/CoarseGraining/adaptivecgproduction.f90`)

- `print_resolved_configuration` (line 1359) — one-time setup banner:
  operator/selector/channel-mode names, block size, `atoms=`/`blocks=`,
  `initial_coarse=`, `interface_atoms=`, and `initial_state_source=` (which
  records whether the run started from `initmag` or a completed
  atomistic-IP handoff). This is what
  `"AdaptiveCG: capability accepted"` and `"initial_state_source=..."` come
  from.
- `print_adaptive_cg_summary` (line 1402) — printed once, after the whole
  run: `completed_steps`, `active_atom_updates`, `active_block_updates`,
  `accepted_transitions`, `rejected_transitions`, `baseline_atom_updates`,
  `reduced_atom_updates`, eight named `last_energy_j` terms
  (`atomistic_bilinear`, `atomistic_onsite`, `coarse_exchange`,
  `coarse_spiralization`, `coarse_anisotropy`, `coarse_external`,
  `coarse_dipole`, and `last_total_energy_j`), four named
  `last_field_checksums_t` sums (`atom_sum`, `atom_norm2`, `coarse_sum`,
  `coarse_norm2`), `owned_cpu_bytes`, and finally
  `trajectory_checksums direction_sum=sum(emom)
  direction_norm2=sum(emom*emom)` — **global aggregate sums over every atom
  and ensemble at the final step only**. This is the line that gives the
  misleading key name `trajectory_checksums`; it is not a trajectory in any
  sense that survives collision (two very different final configurations can
  share the same `direction_sum`).
- `print_transition_events` (line 1448) — per-event log line with step,
  block, old/new state, accepted flag, reason, before/after/jump energy, and
  (since RCG-03) `polarization_ratio`. Exists and is exercised at runtime,
  but `run_production_e2e.py` never parses an individual event, only the
  aggregate counts from the summary line above.
- `print_resolution_state(label, step)` (line 1469) — prints one
  comma-joined vector of every block's resolution state, called with
  `label='final'` only from the production driver path exercised by this
  harness; `final_state()` is the only regex that reads it.

The GPU path (`source/gpu_files/gpuSimulation.cpp:807-980`) mirrors this
exactly under a `Gpu: AdaptiveCG` prefix, including its own aggregate-only
`trajectory_checksums direction_sum=... direction_norm2=...` line (line
968) and its own `final_state values=` / `resolution_counts` /
`polarization_ratio values=` lines — no GPU path prints or could be parsed
for a per-step trajectory today either.

**Conclusion:** every "trajectory" claim in the current suite name or
comment (`trajectory_checksums`, `parity_*` fixture names) is actually a
same-step aggregate-sum or final-restart-state comparison. This is the
concrete finding the RCG-04 prompt pack anticipated in section 3.4/6.2:
"`direction_sum` and `direction_norm2` may remain diagnostics, but they are
not sufficient parity or dynamics oracles" and "do not accept only global
direction sums, squared norms, or final-state checksums as complete
trajectory evidence." The current suite does exactly what that clause
forbids being treated as sufficient.

---

## 3. Fixture-level inventory

All fixtures share `posfile` (2 atoms/cell, positions `(0,0,0)` and
`(0.5,0.5,0.5)`), `momfile` (both atoms `|M|=2.23` pointing `(1,0,0)`,
`Landeg=2.0`), and `jfile` (4 isotropic-exchange shells, no stated
antisymmetric term) unless a row states otherwise. `ncell 6 2 2` with
`block_size 1×2×2` gives 48 atoms in 12 blocks (4 atoms/block) everywhere
except `dmi_anisotropy_mixed[_gpu]` (`ncell 10 2 2`, 80 atoms, 20 blocks) and
the rejection cases. **No fixture in the suite sets `do_tottraj Y`.**
**No fixture computes or asserts an initial torque or field norm.**

Legend for **State class**: `U` = uniform, exactly aligned with `momfile`
(a stationary fixed point of exchange-only dynamics — see §3.5); `U+K` =
uniform, aligned with both exchange and anisotropy easy axis (still a fixed
point); `R` = seeded-random (`Initmag 1`, deterministic given `mseed`, but
not an independently characterized/nonzero-torque construction); `IP` =
relaxed via an initial-phase pass before CG production starts (state after
relaxation is not independently recorded); `SP` = `Initmag 8` analytic spin
spiral (flagged — see §3.5).

### 3.1 Feature-off / static-ownership family

| Fixture (path) | Init & steps | State | Terms & ownership | Observables asserted | Neg. control | Claim category |
| --- | --- | --- | --- | --- | --- | --- |
| `feature_off` (`e2e/feature_off`) | `Initmag 3`, `Nstep 2`, no CG | U | exchange only; CG disabled | `returncode==0`; `"AdaptiveCG:"` absent from stdout; restart state vs `static_all_fine` (48 rows, max abs diff ≤1e-10) | none | **zero-torque smoke**, doing double duty as a weak off/all-fine identity claim |
| `static_all_fine` (`e2e/static_all_fine`) | `Initmag 3`, `Nstep 2`, `cg_operator TENSOR`, STATIC, no mask file (all fine) | U | exchange only; all-fine | `returncode==0`; `"capability accepted"`; `active_atom_updates == baseline_atom_updates`; restart vs `feature_off` | none | **setup smoke** (ownership counting) + **zero-torque smoke** |
| `static_all_coarse` (`e2e/static_all_coarse`) | `Initmag 3`, `Nstep 2`, TENSOR, STATIC, `mask.dat` all-coarse | U | exchange only; all-coarse | `coarse_updates < mixed_updates < fine_updates`; `reduced_atom_updates > 0` | none | **setup smoke** (ownership counting) + **zero-torque smoke** |
| `static_mixed` (`e2e/static_mixed`) | `Initmag 3`, `Nstep 2`, PROJECTED, STATIC, `mask.dat`: block 1 `FINE`, rest default `COARSE` (11/12 blocks) | U | exchange only; genuine fine+coarse+interface split | same ordering assertion above; `"capability accepted"` | none | **setup smoke** (ownership counting, only fixture with a real 3-way ownership split) + **zero-torque smoke** |
| `examples/AdaptiveCoarseGraining/static_mixed` | identical construction to `e2e/static_mixed` | U | exchange only; mixed | `returncode==0`; `"capability accepted"`; `"last_energy_j"` present | none | **setup smoke** (user-facing input regression only) |

### 3.2 Adaptive-selector-mechanics family

All six `cg_refine_threshold`/`cg_coarsen_threshold` pairs in this family
are either **saturated** (`2.00`/`1.99`, i.e. above the `MAX_ANGLE`
selector's physical ceiling, so the criterion always requests coarsening
regardless of the actual state) or set to a fixed small value
(`0.20`/`0.05` for `adaptive_mixed`/`examples/adaptive`) evaluated against a
**uniform, zero-misalignment** state. In both cases the observed
`accepted_transitions > 0` is a direct, provable consequence of the
selector's initialization/threshold logic acting on a trivial input, not of
any texture moving during the run — see §3.5.

| Fixture (path) | Init & steps | State | Terms & thresholds | Observables asserted | Neg. control | Claim category |
| --- | --- | --- | --- | --- | --- | --- |
| `adaptive_mixed` (`e2e/adaptive_mixed`) | `Initmag 3`, `Nstep 2`, `refine 0.20`/`coarsen 0.05` | U | exchange only; ADAPTIVE, no static mask | `accepted_transitions > 0` | none | **zero-torque smoke** (mislabelled by directory name as adaptive-dynamics evidence) |
| `examples/AdaptiveCoarseGraining/adaptive` | same construction | U | same | `"capability accepted"`, `"last_energy_j"` present | none | **setup smoke** |
| `parity_adaptive_cpu`/`parity_adaptive_gpu` (`e2e/parity_adaptive_{cpu,gpu}`) | `Initmag 3` (changed from random by RCG-03 specifically to avoid the polarization gate), `Nstep 3`, saturated `2.00`/`1.99` | U | exchange only; ADAPTIVE | `accepted_transitions > 0`; named energies/field checksums/`final_state`/restart-state CPU vs GPU (rel 5e-4–8e-4, restart abs ≤8e-5) | none | **parity** (final-state + aggregate only; zero-torque — see §3.5) |
| `initial_phase_mc`/`initial_phase_heat_bath` (`ip_mode M`/`H`) | `Initmag 3` start, `ip_mcanneal 1` at 40 steps/`T=100`, then CG `Nstep 2`, saturated thresholds | U → **IP anneal at finite T (genuine stochastic motion, unrecorded)** → CG stage | exchange only | ordering (`"Performing MC initial phase:"` before `"capability accepted"`); `handoff_state_validated`; `direction_norm2≈48` (unit-norm check); `|direction_sum|<40` (loose spread bound) | none | **setup smoke** (IP→CG handoff plumbing only; the anneal's own motion is real but is not independently measured, and the following CG-stage dynamics is not exercised meaningfully under saturated thresholds) |
| `initial_phase_q`/`initial_phase_y`/`initial_phase_z` (`ip_mode Q`/`Y`/`Z`) | `Initmag 3` start, quasi-Newton line-search minimization vs a spiral `qfile`, then CG `Nstep 2`, `refine 0.80`/`coarsen 0.20` | U → **IP minimization (unrecorded)** → CG stage | exchange only | ordering; `handoff_state_validated`; `direction_norm2≈48`; `|direction_sum|<40` | none | **setup smoke** |
| `initial_phase_g` (`ip_mode G`) | `Initmag 2`, `theta0=1.2 rad` (deliberately tilted, not uniform), energy minimization (`min_itrmax 4`), CG `Nstep 2`, saturated thresholds | tilted → **IP minimization (unrecorded)** → CG stage | exchange only | ordering; `handoff_state_validated`; `direction_norm2≈48` | none | **setup smoke** |
| `initial_phase_sd_cpu`/`initial_phase_sd_gpu` (`ip_mode S`) | `Initmag 2`, `theta0=1.2`, `phi0=2π`, one atomistic SD relaxation phase, CG `Nstep 2`, saturated thresholds | tilted → **IP SD relaxation (unrecorded)** → CG stage | exchange only | ordering; `initial_state_source=completed_atomistic_sd`; `accepted_transitions > 0`; `direction_norm2≈48` (CPU tol 2e-12; not compared CPU/GPU) | none | **setup smoke** |
| `initial_phase_mc_gpu`/`initial_phase_q_gpu` | GPU mirrors of the CPU M/S cases above | same as CPU counterpart | same | GPU-specific banner strings (`"GpuMCSimulation: MC initial phase starting"`, `"GpuSDSimulation: SD initial phase starting"`); ordering; `handoff_state_validated` | none | **setup smoke** |

### 3.3 DMI/anisotropy family

| Fixture (path) | Init & steps | State | Terms & ownership | Observables asserted | Neg. control | Claim category |
| --- | --- | --- | --- | --- | --- | --- |
| `dmi_anisotropy_mixed`/`dmi_anisotropy_mixed_gpu` (`ncell 10 2 2`, 80 atoms) | `Initmag 1` (`mseed 73129`), `Nstep 2`, PROJECTED, STATIC, `mask.dat` (mixed) | R | exchange + DMI (`dmfile_cg`, directed ± pairs) + uniaxial anisotropy (`kfile_cg`, easy axis `(0,1,0)`, not aligned with the random state); fine+coarse+interface | CPU: `abs(term) > 1e-30` for `atomistic_onsite`, `coarse_spiralization`, `coarse_anisotropy` only (**nonzero-magnitude check, no sign, no oracle**); GPU: aggregate CPU/GPU parity of 8 energies + 4 field checksums + `final_state` (rel 5e-4/8e-4) | none | **currently unsupported claim** (DMI/anisotropy *dynamics or chirality*; only nonzero-term and weak aggregate CPU/GPU parity are actually asserted) |
| `anisotropy_uniform_fine`/`anisotropy_uniform_coarse` (`kfile_cg_x`, easy axis `(1,0,0)` — aligned with `momfile`) | `Initmag 3`, `Nstep 1`, TENSOR, STATIC (all-fine / all-coarse via `mask.dat`) | U+K | exchange + uniaxial anisotropy, aligned (stationary) | `atomistic_onsite`/`coarse_anisotropy` == closed-form `24·(-0.002-0.003)·2.179872325e-21` J to 2e-12 relative; opposite-ownership term `< 1e-32` | none | **setup smoke** with a genuine independent analytic oracle for a *static* energy value (the strongest oracle-backed check in the whole suite) — explicitly **not** a dynamics or chirality claim; `Nstep 1` performs no integration |

### 3.4 Polarization-gate negative-control family (RCG-03 regression, retained)

| Fixture (path) | Init & steps | State | Terms & thresholds | Observables asserted | Neg. control | Claim category |
| --- | --- | --- | --- | --- | --- | --- |
| `polarization_gate_cpu`/`polarization_gate_gpu` | `Initmag 1` (`mseed 92731`), `Nstep 3`, saturated `2.00`/`1.99` misalignment thresholds **plus** `cg_polarization_threshold 0.9` | R (genuinely incoherent per-block) | exchange only; ADAPTIVE | `accepted_transitions == 0`; final `resolution_counts coarse=0`; `final_state` all nonzero (atomistic) | **is itself the negative control** for `parity_adaptive_cpu` (identical thresholds, aligned `Initmag 3` state, which *does* coarsen) — a true differential test | **rejection** (RCG-03 `POL-CANCELLATION` production regression). The state plausibly *does* move physically for 3 steps (random, unaligned), but no energy/field/trajectory evidence of that motion is recorded or asserted — the fixture proves the gate blocks coarsening, not that dynamics is correct |

### 3.5 Flagged: `initmag_spin_spiral` (`e2e/initmag_spin_spiral`)

`Initmag 8`, `initpropvec (1/6, 0, 0)`, `initrotvec (0,0,1)`, `initrotang 0`,
`Nstep 2`, exchange-only (`jfile`, no DMI/anisotropy/field),
`refine 0.80`/`coarsen 0.20`. This is a **planar** in-plane spin spiral
`m(x) = (cos(qx), sin(qx), 0)` with commensurate `q = 2π/6` matching
`ncell_x=6` (exactly one turn across the supercell). The harness only
asserts the setup banner (`"spin spiral modulation"`,
`"initial_state_source=initmag mode=8"`) and a loose aggregate bound
(`|direction_sum| < 20`, out of a maximum possible 48) — no torque, no
per-step trajectory, no comparison against any moving-state reference.

This is a **near-exact instance of the construction the RCG-04 prompt pack
explicitly warns against** (governing rules §3.4: "Do not use a plain
unperturbed planar Heisenberg spiral as the generic moving state. It can be
a stationary eigenconfiguration"; §6.4 gives the identical planar-spiral
form as the example to avoid). For an isotropic-exchange-only Hamiltonian
with symmetric (±) neighbour shells — which `Sym 1` crystal-symmetry
auto-generation makes likely, though this was not independently proven by
running a new diagnostic in this slice, which is out of RCG-04A's scope —
a planar spin spiral is an exact eigenstate: the local exchange field on
every site remains parallel to that site's own moment, giving zero LLG
torque identically, independent of `q`. **This fixture's actual
nontriviality is unverified and should be treated as suspect** pending an
independent torque check in RCG-04B/C. Category: **currently unsupported
claim** (labelled "spin-spiral"/inhomogeneous, but not demonstrated to be
moving).

The same reasoning applies, with less certainty (uniform + aligned
anisotropy easy axis is unambiguously stationary; a post-IP-relaxation state
in §3.2 is unverified either way), to every `U` and `U+K` row above: none of
them has ever had a torque or field-misalignment norm computed and
asserted, by this harness, at any step.

### 3.6 Static/parity CPU↔GPU family

| Fixture (path) | Init & steps | State | Terms & ownership | Observables asserted | Neg. control | Claim category |
| --- | --- | --- | --- | --- | --- | --- |
| `parity_static_cpu`/`parity_static_gpu` | `Initmag 1` (`mseed 92731`), `Nstep 3`, PROJECTED, STATIC, `mask.dat` mixed | R | exchange only; fine+coarse+interface | `|direction_sum| < 40`; named energies/field checksums/`final_state`/restart-state CPU vs GPU (rel 5e-4/8e-4, restart abs ≤8e-5) | none | **parity** (final-state + aggregate only; random-but-seeded state, not an independently characterized nonzero-torque construction — violates the spirit, though not the unseeded letter, of governing rule §3.4 "do not use an unseeded or implementation-dependent random initial state for claim-bearing parity") |
| `gpu_static_mixed`/`gpu_adaptive_mixed` | mirrors of `static_mixed`/`adaptive_mixed` on `do_gpu Y` | U | same as CPU counterparts | `"Gpu: AdaptiveCG initial active_atoms="`; `0 < static_active < 48`; adaptive `initial`/`final active_atoms` decreasing | none | **setup smoke** / **zero-torque smoke** (GPU mirrors) |
| `gpu_fft_static_mixed` | `Initmag 1`, `Nstep 1`, STATIC (reuses `parity_static_gpu`'s mask), uniform FFT dipole enabled | R | exchange + uniform dipole (FFT/Ewald path); STATIC mixed | `"EWALD3D_FFT production field/energy operator ready"`; `coarse_dipole != 0`; `fft > 0`; `final_state == parity_cpu["static"]`'s state | none | **setup smoke** (dipole-path exercised; no dynamics claim — `Nstep 1`) |

### 3.7 Setup-time rejection cases (not a dynamics claim)

`invalid_partial_block`, `invalid_mask`, `unsupported_temperature`,
`unsupported_initial_phase_x`, `missing_alat` all set `Nstep 0` or `1` and
are asserted only for nonzero return code plus a specific diagnostic
substring (`"block_size_x/y/z"`, `"duplicate block id"`,
`"Temp/do_qhb/do_3tm"`, `"ip_mode"`, `"explicit positive SI lattice
parameter in metres"`). Claim category: **rejection**. No state or motion
claim is or should be made here; these are correctly out of RCG-04's
moving-dynamics scope and require no reclassification.

### 3.8 Untracked/orphaned directory found during this audit

`tests/coarse_graining/e2e/unsupported_initial_phase_mc/` exists on disk but
is **empty** (no files) and is not tracked by git (`git ls-files` returns
nothing under it — an empty directory cannot be tracked) and is not
referenced by `fixture_dependencies.py`, `run_production_e2e.py`, or
`run_setup_rejection_matrix.py`. It is a harmless leftover scratch
directory, not a fixture; no source or test change is needed, but a later
slice's directory listing should not assume every `e2e/*` subdirectory is a
live fixture. Likewise `examples/AdaptiveCoarseGraining/static/` (distinct
from `static_mixed`) is a tracked, git-committed example directory that is
**not** referenced by `EXAMPLE_CASES` or any test — an untested
documentation example, out of this audit's harness scope but noted for
completeness.

---

## 4. Assertions-to-claims summary

| What the suite currently proves | What it does **not** prove |
| --- | --- |
| Setup/CLI plumbing accepts valid inputs and rejects specific invalid ones with the right diagnostic (§3.7). | That any accepted configuration produces correct *dynamics*. |
| Ownership bookkeeping (`active_atom_updates` ordering fine>mixed>coarse; nonzero `reduced_atom_updates`) is self-consistent (§3.1). | That the coarse operator reproduces atomistic *physics* under motion. |
| The RCG-03 polarization gate distinguishes an aligned state (coarsens) from a genuinely incoherent one (does not) under identical thresholds (§3.4). | Anything about the *trajectory* of either state — only the final transition count and resolution vector are checked. |
| A closed-form static uniaxial anisotropy energy matches an independent hand-derived value to 2e-12 (§3.3, `anisotropy_uniform_*`). | Any DMI or anisotropy *dynamical* or *chirality* claim — `dmi_anisotropy_mixed` only checks that three named terms are nonzero. |
| CPU and GPU production runs agree on final restart state, final resolution state, and aggregate energy/field checksums within loose tolerances (5e-4–8e-4 relative; restart 8e-5 absolute) for four fixture pairs (§3.2, §3.6). | Per-step agreement, or agreement on any *moving* state — every CPU/GPU comparison in the suite starts from either a uniform (zero-torque) or seeded-random-but-uncharacterized state. |
| The off/all-fine restart-state identity (`feature_off` vs `static_all_fine`, 1e-10 absolute) is a genuine equality check. | That this equality survives a *nonstationary* state — the compared state is the exact uniform exchange ground state, so both paths trivially agree. |
| `accepted_transitions > 0` fires for several "adaptive"-named fixtures. | That any transition was caused by real spatial motion — every such fixture uses either a saturated threshold pair or a uniform/zero-misalignment state (§3.2, §3.5), so the transition is a threshold/initialization artifact, not moving-dynamics evidence. |

**No current fixture in `run_production_e2e.py` supports a "moving
dynamics" claim category as defined by the RCG-04 prompt pack.** Every case
falls into setup smoke, zero-torque smoke, rejection, weak (final-state/
aggregate-only) parity, or currently-unsupported claim. This confirms the
prompt pack's own baseline assessment (§1): "Existing feature-off, all-fine,
all-coarse, static-mixed, and adaptive-mixed fixtures remain dominated by
aligned, stationary, very short, or final-state evidence... they cannot
carry moving-dynamics claims."

Per the instruction to avoid broad renaming churn, no fixture directory,
CTest name, or `inpsd.dat` file was renamed in this slice. The
reclassification above is recorded here, in the evidence document, as the
authoritative labelling; §7 also proposes light, additive
docstring/comment-only relabelling for RCG-04B to apply alongside the real
fixture work, rather than as a separate churn-only patch.

---

## 5. Canonical moving-evidence record (schema for RCG-04B onward)

Later slices must produce a machine-readable evidence record per accepted
moving fixture (JSON/YAML sidecar or an equivalently structured stdout
block) containing at minimum:

1. **Provenance:** source commit, `git status --short`, generator
   name/version and parameters (RCG-04B), input content hash, configure/
   build/run commands, compiler, backend, device, precision.
2. **Per-spin trajectory or justified equivalent:** either the complete
   `(step, atom, ensemble, direction, moment)` series — reusing
   `do_tottraj`/`moment.<simid>.out` per the governing blueprint's "prefer
   existing ordinary production output" rule (§3.3) wherever it satisfies
   the contract — or an equivalently robust representation (e.g. a
   collision-resistant per-step hash keyed by `(step, atom, ensemble)` with
   an explicit, tested collision/ordering argument). A same-step aggregate
   sum/norm, as used throughout the current suite, does **not** qualify.
3. **Named energy and field records:** the existing eight `last_energy_j`
   terms and four `last_field_checksums_t` sums, but recorded **per sampled
   step**, not only at the final step.
4. **Restart-state comparison:** the existing `restart.<simid>.out`
   final-state comparison, retained as one (not the only) check.
5. **Transition history where applicable:** the full per-event
   `print_transition_events` record (step, block, old/new state, accepted,
   reason, before/after/jump energy, polarization ratio), not only the
   aggregate `accepted_transitions` count.
6. **Nonzero-torque and nonzero-displacement evidence:** an independently
   computed (not read back from the production diagnostic under test)
   initial effective-field/torque norm, and a maximum/final per-spin
   displacement from the initial state, each checked against a documented,
   non-zero floor — directly answering this slice's central finding that no
   current fixture does this.
7. **Physical mode-specific observables:** conical-mode amplitude/phase/
   frequency, signed chirality, domain-wall centre/displacement/crossing,
   as applicable to the specific slice (RCG-04D–H).
8. **Precision:** which of fp32/fp64, and (RCG-04I) the backend-specific
   tolerance and its justification — explicitly not decided in this slice.

This schema is a specification for RCG-04C to implement, not new code
introduced here.

---

## 6. Required raw-evidence fields (evidence policy, restated for RCG-04)

Every later RCG-04 slice's claim-bearing run must record: exact commit and
dirty/clean status; configure/build/run/test commands; compiler, backend,
device, precision; fixture/generator provenance; observable definitions and
units; raw/machine-readable result location; the tolerance used and why;
negative-control modification and its expected failure; restoration of the
unmodified source after any fault injection; and any missing backend or
deferred review — per
`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` §2.2 and the
prompt pack §3.3, both already governing and unchanged by this slice.

---

## 7. Later-slice responsibility map

| Slice | Responsibility, informed by this inventory |
| --- | --- |
| RCG-04B | Deterministic conical-spiral (not planar — see §3.5's flag on `initmag_spin_spiral`), domain-wall-pair, and ±q DMI generators, each with an explicit nonzero-torque rejection/label for degenerate parameter choices. |
| RCG-04C | Reuse `do_tottraj`/`moment.<simid>.out` (currently `N` everywhere, §2.1) and `print_transition_events` (currently unparsed, §2.1) as the primary trajectory/transition sources; implement the §5 schema's comparison metrics. |
| RCG-04D | Redo the `feature_off`/`static_all_fine` identity claim (§3.1) starting from a genuinely nonstationary state instead of the current uniform exchange ground state; retain the current pair as a labelled zero-torque smoke test per the prompt pack. |
| RCG-04E | All-coarse dynamics against an atomistic reference — `static_all_coarse` today (§3.1) only checks update-count bookkeeping, not physics. |
| RCG-04F | Static-mixed interface dynamics — `static_mixed` (§3.1) already has a genuine 3-way fine/coarse/interface split via `mask.dat`, which is reusable topology, but its state must become moving. |
| RCG-04G | Adaptive boundary-crossing — every current "adaptive" fixture's transitions are threshold/initialization artifacts (§3.2), not motion-driven; `polarization_gate_cpu/gpu` (§3.4) is the one fixture with a genuinely incoherent state and should be examined as a possible basis for a real moving/negative-control pair, since its motion, while plausible, is currently unrecorded. |
| RCG-04H | DMI/anisotropy chirality — `dmi_anisotropy_mixed[_gpu]` (§3.3) today only proves three terms are nonzero; no sign, oracle, or dynamics is checked, and its DMI file/anisotropy easy axis are not stated as displaced from the DMI minimum. |
| RCG-04I | CPU/GPU parity — the four existing parity pairs (§3.2, §3.6) are final-state/aggregate-only with 5e-4–8e-4 relative and 8e-5 absolute tolerances chosen without documented derivation; RCG-04I must derive and freeze real budgets from moving-state error scaling. |
| RCG-04J | Reconcile all of the above against the five required exit-evidence packages; this document's §4 table is the authoritative "before" baseline to cite. |

---

## 8. Fresh build/test evidence (this slice)

**Environment:** GNU Fortran 13.3.0, GNU C/C++ 12.4.0, CMake 3.28.3, CPU
backend, fp64 (default precision), Release build type.

**Fixture-dependency audit** (run first, no build required):

```text
$ python3 tests/coarse_graining/audit_fixture_dependencies.py
adaptive-CG fixture dependency audit: PASS (38 fixture directories, 60 input paths)
```

**Fresh out-of-tree configure/build:**

```text
$ cmake -S . -B /tmp/rcg04a-cpu-uZjr1I -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
...
-- Git tag found: VERSION="v6.0.2-448-ge382-dirty".   # dirty solely from the
                                                        # pre-existing, unrelated
                                                        # anomaly recorded above
-- Configuring done
$ cmake --build /tmp/rcg04a-cpu-uZjr1I -j$(nproc)
...
[100%] Built target polarization_gate_tests            # exit 0, no warnings
                                                          other than one
                                                          pre-existing
                                                          executable-stack
                                                          linker warning
                                                          unrelated to this
                                                          slice
```

**Test run (`cg13-cpu` label):**

```text
$ ctest --test-dir /tmp/rcg04a-cpu-uZjr1I -L cg13-cpu --output-on-failure
...
100% tests passed, 0 tests failed out of 12
 3: coarse-graining-block-topology .............. Passed
 4: coarse-graining-stiffness-material .......... Passed
 5: coarse-graining-dmi-dimer-energy ............ Passed
 6: coarse-graining-tensor-operator ............. Passed
 7: coarse-graining-multichannel-tensor-operator  Passed
 8: coarse-graining-smooth-projected-operator ... Passed
 9: coarse-graining-static-hybrid-operator ...... Passed
10: coarse-graining-block-selector .............. Passed
11: coarse-graining-adaptive-hybrid-solver ...... Passed
12: coarse-graining-polarization-gate ........... Passed
13: adaptive-cg-production-e2e .................. Passed (1.02 sec)
14: adaptive-cg-setup-rejection-matrix .......... Passed (2.43 sec)
```

No CUDA/HIP evidence was gathered in this slice (not required for an
inventory-only task; RCG-04D onward will need fresh backend evidence for
their own fixtures).

**Worktree check after the run:**

```text
$ git status --short --porcelain=v1 | grep -v '^??'
 M docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md   # pre-existing, see above; untouched by this slice
```

`adaptive-cg-production-e2e` runs the real binary against tracked
`examples/AdaptiveCoarseGraining/*` inputs, which receive their own
provenance-stamped `uppasd.*.yaml` output as a side effect of testing (the
same behaviour RCG-03's evidence already documented); those were restored
with `git checkout` after the run so no test-run byproduct remained staged
or modified. No other tracked or untracked file was created, modified, or
deleted by this slice other than this new evidence document.

---

## 9. RCG-04A checklist

- [x] Current e2e harness, fixtures, registrations, and diagnostics are inventoried. (§2, §3)
- [x] Each current fixture's initial state and actual motion evidence are identified. (§3.1–§3.8, with `initmag_spin_spiral` flagged in §3.5)
- [x] Existing assertions are mapped to the claims they can genuinely support. (§4)
- [x] Uniform/aligned fixed-point cases are labelled smoke or zero-torque only. (§3.1, §3.2, §3.3, §3.6 claim-category columns)
- [x] No current stationary case is presented as moving-dynamics evidence. (§4: "No current fixture... supports a 'moving dynamics' claim category")
- [x] The canonical moving-evidence observable schema is documented. (§5)
- [x] Required provenance and raw-evidence fields are documented. (§6)
- [x] Later responsibility is mapped to RCG-04B through RCG-04J. (§7)
- [x] No final tolerance has been fitted or accepted in this inventory slice. (§5 explicitly defers precision/tolerance to RCG-04I; no numeric tolerance was changed anywhere)
- [x] Fixture dependency audit still passes, or any pre-existing failure is isolated. (§8: PASS, 38/60)
- [x] Fresh configure/build/test commands and results are recorded where run. (§8: CPU configure/build/`cg13-cpu` 12/12)
- [x] Unrelated worktree changes remain untouched and unstaged. (§8; pre-existing `REMEDIATION_BLUEPRINT.md` anomaly identified, reported to the user, and left untouched by explicit user decision)

All twelve RCG-04A checklist items are complete and evidenced.

---

## 10. RCG-04B: deterministic moving-state generators

**Status: RCG-04B only. Generator and provenance layer — no fixture,
production diagnostic, moving-dynamics claim, or production-consuming test
run was added in this slice.**

**Base commit:** `3fe7d600c98f52c6238fe4558dae11ac24d5d5fa` ("RCG-04A: define
moving e2e evidence contract"), the accepted RCG-04A commit. `git status
--short` at session start showed no modified tracked files (the
`docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` anomaly recorded in
§0 above was fixed, per explicit user instruction, before the RCG-04A commit
and is no longer present).

### 10.1 Inspection of existing conventions (before implementation)

Before writing any generator code, the following were read from source:

- `source/Input/inputhandler_ext.f90:read_moments` — **`momfile` is indexed
  by `(isite, ichem)`, i.e. one row per *basis site* (`na`), not one row per
  global atom.** It is a per-basis-site template, broadcast identically to
  every unit cell during setup. This means an ordinary `momfile` alone
  cannot express a spatially varying state — matches and explains this
  document's §3.1/§3.2 finding that every `Initmag 3` fixture is uniform.
- `source/System/magnetizationinit.f90:402-524` (`initmag==8`) — UppASD's
  *existing* spin-spiral initializer rotates the `momfile` base vector for
  each atom by an angle `q.x + phase` about a caller-chosen axis
  (Rodrigues' rotation formula, lines 469-486). Critically, **the base
  vector need not be perpendicular to the rotation axis** — the existing
  `initmag_spin_spiral` fixture (flagged in §3.5) happens to choose one that
  is (giving a planar spiral), but the same mechanism, given a base vector
  at a genuine cone angle, produces an honest conical spiral with no source
  change at all.
- `source/System/magnetizationinit.f90:267-282` (`initmag==4`) and
  `source/System/restart.f90:read_mag_conf_std` — ordinary UppASD's restart
  loader reads a full per-atom configuration (7-line header, then
  `iter ens iatom |M| Mx My Mz` per atom, matching the writer
  `prn_mag_conf_iter`/`prn_mag_conf_time` exactly). This is the only
  existing mechanism general enough to express a domain wall.
- `source/CoarseGraining/adaptivecgproduction.f90:692-696`
  (`validate_configuration`) — **AdaptiveCG production accepts `Initmag` in
  `{1, 2, 3, 5, 8}` only; `Initmag=4` (restart) is explicitly rejected**
  with the diagnostic `"initmag=4 restart is unsupported until adaptive
  state serialization is implemented"`. This is a real, pre-existing
  capability boundary discovered while designing the domain-wall generator,
  not a defect introduced by this slice — see §10.4.
- `source/System/geometry.f90:445` (`coord(1,i)=I1*C1(1)+I2*C2(1)+I3*C3(1)+Bas(1,I0)`)
  and the atom-index formula `i=I0+I1*NA+I2*N1*NA+I3*N2*N1*NA` used
  throughout `magnetizationinit.f90` — the exact ordering/position
  conventions the generator must reproduce to be directly consumable.
- `tests/dipole_validation/generate_cases.py` and
  `tests/dipole_validation/test_open_fft_oracle.py` — the repository's
  existing convention for a deterministic Python generator/oracle module
  plus a `unittest`-based, directly-CTest-registered test file
  (`python3 <test>.py -v`). RCG-04B follows this convention rather than
  inventing a new one (no pytest dependency is introduced).

No isolated test-only momfile/atom-ordering convention was created; both
new generator mechanisms reuse an existing UppASD input path exactly.

### 10.2 What was implemented

`tests/coarse_graining/moving_state_generator.py` (new module):

- `Geometry` — a small dataclass reproducing the atom-index and
  Cartesian-position formulas above, restricted to the identity-`cell`
  case (`cell 1 0 0 / 0 1 0 / 0 0 1`), which is every geometry file under
  `tests/coarse_graining/e2e` today. This restriction is stated explicitly
  in the module docstring rather than silently assumed.
- `conical_spiral_state(geometry, *, cone_angle_deg, turns, axis,
  phase_deg, modulation_cell_axis, moment_magnitude, landeg,
  degeneracy_tolerance_deg)` — produces a `momfile` plus the
  `Initmag 8`/`initpropvec`/`initrotvec`/`initrotang` `inpsd.dat` records
  needed to drive UppASD's existing spin-spiral initializer into a genuine
  cone. Also returns the full per-atom direction dictionary (computed
  independently in Python, for self-testing and for later oracle reuse) and
  a `manifest` provenance record.
- `chiral_partner_pair(...)` — calls `conical_spiral_state` twice with
  `turns` and `-turns` (all other parameters identical); neither function
  accepts or reads a DMI file, so the `+q`/`-q` choice cannot be derived
  from whichever DMI input a later fixture tests against (checked
  structurally in `test_does_not_depend_on_any_dmi_file`).
- `domain_wall_pair_state(geometry, *, axis_cell_index, wall_centers_cells,
  width_cells, easy_axis, wall_type, chirality, cant_deg,
  moment_magnitude, simid, separation_margin_widths)` — produces a
  restart-format text block for a periodic **two**-kink domain-wall pair
  (`theta(u) = 2*atan(exp((u-c1)/w)) - 2*atan(exp((u-c2)/w))`, giving an
  up-down-up pattern with zero net winding, required for periodic boundary
  conditions), with `NEEL`/`BLOCH` wall-plane selection, a `chirality`
  sign, and a small deterministic `cant_deg` tilt (scaled by `sin(theta)`
  so it vanishes away from either wall core) that breaks exact wall-plane
  symmetry pinning without introducing any randomness.
- `content_hash`/`manifest_json` — SHA-256 content hash and canonical
  (`sort_keys=True`) JSON serialization for provenance/regeneration checks.

**Analytic justification for the conical-spiral degeneracy rejection**
(independently derived from the written Hamiltonian in the module
docstring, not inferred from any production diagnostic): for a Bravais
lattice with isotropic, centrosymmetric pairwise exchange, a conical-spiral
ansatz `m(x) = R_axis(q.x).s0` is an exact LLG fixed point iff `s0` is an
eigenvector of `M(q) = sum_delta J(delta) R_axis(q.delta)`. In the rotation
eigenbasis, `M(q)` is diagonal with eigenvalue `J0 = sum_delta J(delta)`
along the axis and `J(q) = sum_delta J(delta) cos(q.delta)` (real, by
centrosymmetry) in the perpendicular plane. A cone angle strictly between 0
and 180 degrees, and not equal to 90, mixes eigenvectors with different
eigenvalues whenever `J0 != J(q)` — generically true for any nonzero `q` in
the long-wavelength regime, where `J0-J(q) ~= q^2 * (exchange stiffness) /
2`. `conical_spiral_state` therefore rejects cone angles within
`degeneracy_tolerance_deg` (default 5 degrees) of 0/90/180 degrees, and
rejects `turns=0` outright, raising `DegenerateStateError` with the
reasoning inline. This is the same physical mechanism §3.5 of this
document flagged as a suspected defect in the existing
`initmag_spin_spiral` fixture (a planar, `theta=90`, construction);
RCG-04B's generator makes the non-degenerate region the default path
rather than special-casing the flagged fixture itself (no fixture files
were touched in this slice).

`domain_wall_pair_state` similarly rejects (`DegenerateStateError`) wall
placements whose minimum separation (from each other and from the periodic
boundary) is less than `separation_margin_widths` (default 4) times
`width_cells`, since the two-kink profile only closes periodically up to a
residual of order `exp(-separation/width)`.

### 10.3 Tests (`tests/coarse_graining/test_moving_state_generator.py`)

40 `unittest` cases, run directly (`python3 ... -v`) and via the new CTest
target `coarse-graining-moving-state-generator` (labels
`coarse-graining;cg13;cg13-cpu;reference`, and appended to `cg13-cuda`/
`cg13-hip` alongside the other host-only reference tests when those
backends are configured, since generation has no backend dependency),
covering every category the RCG-04B prompt requires:

- **atom ordering:** `iter_atoms()` is a permutation of `1..natom`; its
  nesting order matches `magnetizationinit.f90` exactly (basis fastest,
  then `i1`, `i2`, `i3`); `atom_index` matches the production formula by
  direct computation.
- **normalization:** every generated direction, for both generators, has
  unit norm to at least 1e-10.
- **periodic closure:** the conical spiral's azimuth advances linearly by
  exactly `turns*2*pi/n1` per cell step (12 assertions across one full
  supercell traversal); the domain-wall pair's polar angle at the two ends
  of the periodic axis agrees to within the analytically expected
  `exp(-separation/width)` residual (a loosened, not exact, tolerance —
  documented inline as an expected discretization effect, not evidence of
  a construction error).
- **wall count/placement:** exactly two sign changes of the easy-axis
  component occur across the periodic domain, located within 1.5 cells of
  the requested centres.
- **opposite chiral partners:** `+q`/`-q` partners share an identical cone
  (axis-component) profile but wind in opposite azimuthal directions; a
  structural test confirms neither generator function's signature accepts
  a DMI-file argument.
- **deterministic regeneration:** two independent calls with identical
  parameters produce byte-identical `momfile`/restart text and an identical
  recorded SHA-256 hash; changing one parameter changes the hash; the
  manifest's JSON serialization is itself byte-stable across calls.
- **malformed input:** mismatched basis/`na` length, non-positive
  geometry extents, out-of-range cell indices, invalid `modulation_cell_axis`/
  `axis_cell_index`, non-positive `moment_magnitude`/`width_cells`,
  unordered/out-of-range wall centres, invalid `wall_type`/`chirality`, and
  an `easy_axis` parallel to the wall propagation axis all raise a typed
  error (`MalformedGeneratorInputError`) with a specific message.

No RNG, `set()` iteration, or filesystem/locale dependency is used anywhere
in the generator, satisfying the determinism requirement without needing a
fixed seed (there is nothing stochastic to seed).

### 10.4 Production capability gap found, then fixed by explicit user direction

While designing the domain-wall-pair generator, AdaptiveCG production's
explicit rejection of `Initmag=4` (§10.1, `adaptivecgproduction.f90:692-693`)
was found to block **any** genuinely spatially-varying, non-helical initial
state — not just this generator's output — from ever being loaded into an
AdaptiveCG-enabled production run. This is a real, load-bearing scope
boundary for RCG-04G (adaptive wall motion), which depends on this
generator per the RCG-04 dependency graph. This was reported to the user
(not silently worked around) before any fix was attempted, per the
remediation blueprint's governing rule that a discovered production
capability gap "must be demonstrated and reported... before expanding the
slice."

**The user explicitly asked for the rejection to be fixed** as an
in-session follow-up, after reviewing this finding. The remainder of this
section records that fix, which is therefore in scope for this commit
despite being outside RCG-04B's originally drafted prompt.

#### 10.4.1 Why lifting the rejection is safe

The rejection message read: `"initmag=4 restart is unsupported until
adaptive state serialization is implemented"`. Before changing anything,
the call sequence in `source/uppasd.f90` was traced to check whether that
concern actually applies to this use case:

1. `magninit` (line 1245) populates `emom`/`emomM` — including, for
   `initmag==4`, via the ordinary `read_mag_conf`/`read_mag_conf_std`
   restart loader — **before** any AdaptiveCG code runs at all.
2. `preflight_adaptive_cg_production` (called later, line ~1375) only
   validates the configuration; it does not touch `emom`.
3. `setup_adaptive_cg_production` (called later still, after any initial
   phase, per its own comment: "Construct adaptive ownership only now,
   from the completed atomistic handoff state") always builds a **fresh**
   block/channel classification from whatever `emom` currently holds — the
   same construction runs identically regardless of which `Initmag` value
   produced that `emom` state.

There is no code path, for any `Initmag` value, that resumes a *previous
AdaptiveCG run's own* resolution/dwell/transition-history state — that
capability (the genuine meaning of "adaptive state serialization") simply
does not exist yet for anything. The rejection therefore did not protect
against an actual functional difference between `Initmag=4` and the five
already-accepted values; it blocked a case that is architecturally
identical to them (restart is just another way to seed the same cold-start
atomistic state that `Initmag` 1/2/3/5/8 already seed). `Initmag=4` was
added to the accepted set in `validate_configuration`
(`source/CoarseGraining/adaptivecgproduction.f90`), replacing the outright
rejection, with the reasoning above recorded inline as a source comment.

#### 10.4.2 Regression: the now-obsolete rejection-matrix negative control

`tests/coarse_graining/run_setup_rejection_matrix.py` had a `"restart"`
case asserting `initmag=4` rejection. Inspecting it surfaced a second,
independent, pre-existing fragility unrelated to this fix: the case relied
on `tests/coarse_graining/e2e/static_mixed/restart.cg105mix.out`, which
`git ls-files` shows is **not tracked** (it matches the repository's
`restart.*.out` `.gitignore` pattern) — a leftover runtime artifact from a
previous local test run, not a reproducible fixture input. On a genuinely
clean clone this case would already have failed differently (a generic
"restartfile does not exist" stop, never reaching AdaptiveCG's own check),
independent of today's fix. This case has been removed (not
reworked to assert something else), since there is no longer any rejection
to assert; its removal incidentally also removes that latent fragility.

#### 10.4.3 New positive regression: `initmag_restart_atomistic`

`tests/coarse_graining/e2e/initmag_restart_atomistic/` is a new production
e2e fixture (full provenance in that directory's `README.md`) proving the
capability now works end to end through the real `sd.f95` executable: the
standard 48-atom host geometry, `Initmag 4`, and a restart-format seed file
(`restart_seed.out`) generated deterministically by this slice's own
`domain_wall_pair_state` (RCG-04B generator consuming its own output,
closing part of the "no fixture yet consumes this generator" gap noted
below). `restart_seed.out` is deliberately **not** named
`restart.<simid>.out`, and `simid` was kept within UppASD's
`character(len=8)` field width — both found and fixed during development,
see the fixture's `README.md` for the exact silent-truncation and
self-overwrite hazards that motivated them. `run_production_e2e.py` now
asserts `returncode == 0`, `"AdaptiveCG: capability accepted"`, and
`"AdaptiveCG: initial_state_source=initmag mode=4"`. This is a
**setup/capability smoke test only** (`Nstep 1`, no trajectory, energy, or
field assertion beyond the standard capability-accepted banner) — it does
not claim the resulting dynamics is correct, consistent with RCG-04B's "no
moving-dynamics, parity, or accuracy claim" boundary.

`domain_wall_pair_state`'s output remains, additionally, directly
consumable by **ordinary** (`AdaptiveCG`-disabled) UppASD, useful for a
future feature-off/all-fine reference.

### 10.5 RCG-04B checklist

- [x] Existing input/generator conventions were inspected before implementation. (§10.1)
- [x] A deterministic periodic conical-spiral generator is implemented. (§10.2, `conical_spiral_state`)
- [x] Special stationary parameter choices are rejected or explicitly identified. (§10.2 analytic justification; `DegenerateStateError`)
- [x] A deterministic periodic domain-wall-pair generator is implemented. (§10.2, `domain_wall_pair_state`)
- [x] A deterministic `+q`/`-q` DMI-sensitive state generator is implemented. (§10.2, `chiral_partner_pair`)
- [x] Atom ordering, basis/material identity, and moment magnitudes are preserved. (§10.1, §10.3 atom-ordering tests)
- [x] Every generated spin direction is normalized within a documented budget. (§10.3, 1e-10)
- [x] Periodic closure and wall topology are tested. (§10.3)
- [x] Repeated generation is byte-stable and produces a recorded stable hash. (§10.3, `DeterministicRegenerationTests`)
- [x] Generator parameters and provenance are tracked in a manifest or equivalent. (§10.2 `manifest`/`manifest_json`; `GENERATOR_MANIFEST.json` in §10.4.3's fixture)
- [ ] CPU and GPU fixture paths can consume identical generated bytes. Partially evidenced on CPU only: §10.4.3's `initmag_restart_atomistic` fixture proves CPU production consumes this generator's restart-format output; no GPU-path or conical-spiral-path fixture exists yet, so this box remains open rather than fully ticked.
- [x] Malformed or physically incompatible generator requests fail clearly. (§10.3 malformed-input tests)
- [x] Tracked-fixture/package audit covers the generator and required inputs. `restartfile` was added to `audit_fixture_dependencies.py`'s tracked-input-keyword set, and §10.4.3's fixture (whose `restart_seed.out` is this generator's output) is now covered: `audit_fixture_dependencies.py` passes at 39 fixture directories / 62 input paths (§10.6).
- [x] No moving-dynamics, parity, or accuracy claim is made in this slice. (§10 status line and §10.4.3: the new fixture is explicitly a setup/capability smoke test, `Nstep 1`)
- [x] Unrelated worktree changes remain untouched and unstaged. (§10.6)

### 10.6 Fresh build/test evidence (this slice)

**Environment:** GNU Fortran 13.3.0, GNU C/C++ 12.4.0, CMake 3.28.3, CPU
backend, fp64, Release build type.

```text
$ python3 tests/coarse_graining/test_moving_state_generator.py -v
...
Ran 40 tests in 0.016s
OK

$ cmake -S . -B /tmp/rcg04b-fix-cpu-Nt2NEM -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
-- Git tag found: VERSION="v6.0.2-449-g3fe7-dirty".   # dirty from this
                                                        # slice's own
                                                        # uncommitted
                                                        # changes, expected
                                                        # mid-slice
$ cmake --build /tmp/rcg04b-fix-cpu-Nt2NEM -j$(nproc)
...
[100%] Built target polarization_gate_tests   # exit 0 (rebuild after the
                                                # adaptivecgproduction.f90 fix)

# Manual pre-CTest check that the new fixture works end to end and that its
# tracked restart_seed.out is not overwritten by the run:
$ cd tests/coarse_graining/e2e/initmag_restart_atomistic && \
  md5sum restart_seed.out > /tmp/before.md5 && \
  /tmp/rcg04b-fix-cpu-Nt2NEM/bin/sd.f95 > /tmp/out2.log 2>&1; echo "returncode=$?"; \
  md5sum -c /tmp/before.md5; grep -n "capability accepted\|initial_state_source" /tmp/out2.log
returncode=0
restart_seed.out: OK
45:AdaptiveCG: capability accepted: regular periodic single-FM deterministic Heun
52:AdaptiveCG: initial_state_source=initmag mode=4

$ ctest --test-dir /tmp/rcg04b-fix-cpu-Nt2NEM -L cg13-cpu --output-on-failure
...
13: coarse-graining-moving-state-generator ...... Passed (0.05 sec)
14: adaptive-cg-production-e2e .................. Passed (0.75 sec)
15: adaptive-cg-setup-rejection-matrix .......... Passed (2.34 sec)
100% tests passed, 0 tests failed out of 13

$ ctest --test-dir /tmp/rcg04b-fix-cpu-Nt2NEM -R '^asd-tests$' --output-on-failure
1/1 Test #2: asd-tests ........................   Passed    7.01 sec
100% tests passed, 0 tests failed out of 1    # legacy feature-off/non-CG
                                                # regression suite, confirming
                                                # the fix and new fixture
                                                # introduced no unrelated
                                                # regression

$ python3 tests/coarse_graining/audit_fixture_dependencies.py
adaptive-CG fixture dependency audit: PASS (39 fixture directories, 62 input paths)
```

**Worktree check after the run:** `adaptive-cg-production-e2e` again
touched the three tracked `examples/AdaptiveCoarseGraining/*/uppasd.adaptive.yaml`
provenance files as a side effect of testing (same behaviour documented in
§8 and in the RCG-03 evidence); they were restored with `git checkout`
after the run. Manually running the new fixture also produced its own
runtime byproducts (`restart.cg105ir4.out`, `inp.cg105ir4.json`,
`uppasd.cg105ir4.yaml`, all matching existing `.gitignore` patterns), which
were deleted (not committed) after confirming they were byproducts, not
inputs. `restart_seed.out` itself required `git add -f` to stage, since it
matches the repository's broad `*.out` ignore pattern despite being a
genuine tracked fixture input (consistent with `restart.*.out` already
being an explicit e2e-scoped ignore pattern for the same reason). After
restoration and staging, `git status --short` showed exactly this slice's
intended files: `CMakeLists.txt`,
`source/CoarseGraining/adaptivecgproduction.f90`, five modified
`tests/coarse_graining/*.py` harness files, two new generator/test `.py`
files, the four new `initmag_restart_atomistic` fixture files, and this
evidence document. `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
(fixed before the RCG-04A commit) remained untouched and clean throughout.

---

## Open items (carried forward, not blocking RCG-04A/B)

- ~~AdaptiveCG's `Initmag=4` rejection blocks RCG-04G~~ — **fixed** in this
  slice (§10.4): `validate_configuration` now accepts `Initmag=4`, and
  `initmag_restart_atomistic` (§10.4.3) proves it end to end on CPU. RCG-04G
  can now consume `domain_wall_pair_state`'s output in an AdaptiveCG-enabled
  run; this no longer blocks it. GPU-path evidence for `Initmag=4` remains
  unestablished (no GPU fixture was added) and CUDA/HIP builds were not
  exercised in this slice — a later slice (RCG-04G or RCG-04I) should
  confirm the GPU dispatch path also accepts and honours a restart-loaded
  state before relying on it there.
- Whether `initmag_spin_spiral` truly has zero initial torque (§3.5) is
  argued from Hamiltonian symmetry, not from a freshly run diagnostic
  (RCG-04A/B do not implement one). RCG-04B's own analytic derivation
  (§10.2) independently confirms the mechanism (planar, `theta=90`, cone
  angle under isotropic centrosymmetric exchange is exactly the rejected
  degenerate case) but this remains a source-level argument, not a
  production measurement; RCG-04C should verify it independently before
  deciding whether that fixture is reusable, needs a nonzero cant, or
  should be retired in favour of `conical_spiral_state`.
- No fixture yet consumes any RCG-04B generator's output through the
  UppASD executable (CPU or GPU); the "CPU and GPU fixture paths can
  consume identical generated bytes" and "tracked-fixture/package audit
  covers the generator" checklist items are explicitly left open above
  until RCG-04D or a later slice wires a generator into a real e2e case.
