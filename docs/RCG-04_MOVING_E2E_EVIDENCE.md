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

## Open items (carried forward, not blocking RCG-04A)

- The pre-existing, out-of-scope corruption in
  `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` (see top of this
  document) remains in the working tree by explicit user decision
  (2026-08-08) and should be reverted or repaired in a separate, unrelated
  change whenever the user chooses.
- Whether `initmag_spin_spiral` truly has zero initial torque (§3.5) is
  argued from Hamiltonian symmetry, not from a freshly run diagnostic in
  this slice (implementing such a diagnostic is explicitly out of RCG-04A's
  scope). RCG-04B/C must verify this independently before deciding whether
  the fixture is reusable, needs a nonzero cant, or should be retired in
  favour of the conical-spiral generator.
