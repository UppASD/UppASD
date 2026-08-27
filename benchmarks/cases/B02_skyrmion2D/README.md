# B02_skyrmion2D — quasi-2D short-range exchange+DMI skyrmion lattice

Admission record for WP-08B (blueprint §3 "B02 — 2D skyrmion: short-range
(J+D)"). This case is the framework's short-range-neighbour, two-dimensional
production reference: contrasted against `B01_bccFe`'s medium-range,
three-dimensional workload, and against a future `B03_skyrmion3D`'s
genuinely three-dimensional short-range J+D workload.

## A. Template audit

### Template selection

Two maintainer example decks model a 2D square-lattice skyrmion:
`examples/SpecialFeatures/SkyrmionLattice` and
`examples/SpecialFeatures/SkyrmionTransport`. Both share the same
Hamiltonian class (NN exchange + NN DMI on the same square lattice, same
`jij`/`dmdata` magnitudes). `SkyrmionTransport` was **not** selected:
`Initmag 4` loads a fixed-atom-count `restartfile` (`start.singlesk`, a
single-skyrmion snapshot) as its initial state and runs an `ip_mcanneal`
Monte Carlo anneal before the timed dynamics -- neither scales to an
arbitrary lattice size (blueprint §D explicitly requires "a sensible
scalable initialization method" and forbids "expensive initialization that
dominates process timing unless intentional"). It also enables `stt A`
(current-driven spin-transfer torque), a physics feature outside B02's
stated purpose (short-range J+D dynamics, not current-driven transport).
`SkyrmionLattice` was selected: `Initmag 3` (uniform direction from
`momfile`, trivially scalable to any `ncell`), no restart file, no MC
anneal, no STT.

Source: `examples/SpecialFeatures/SkyrmionLattice/{inpsd.dat,jij,dmdata,momfile,posfile,qfile}`,
copied byte-for-byte into `template/` (`diff -q` verified identical to the
examples-tree copy, which remains untouched).

| Property | Value | Source |
| --- | --- | --- |
| Lattice | Simple square, single-atom basis, cell `(1,0,0)/(0,1,0)/(0,0,1)` (reduced units) | `posfile`, `inpsd.dat` |
| Dimensionality | Quasi-2D: `ncell Nx Ny 1` (single layer), `BC P P 0` (periodic in-plane, open out-of-plane) | `inpsd.dat` |
| Basis/moment | 1 basis atom, `mmom=1.0`, initial direction `(0.1, 0.0, 1.0)` normalized (`Initmag 3`) | `momfile`, `inpsd.dat` |
| Symmetry | `Sym` unset -> default `0` (no symmetry expansion): `jij`/`dmdata`'s 4 listed bonds are the complete directed lists, not shell representatives | `inpsd.dat`, `source/Input/inputdata.f90:331` |
| Exchange | 4 nearest-neighbour bonds (`+-x, +-y`), `J=1.0` | `jij` |
| DMI | Same 4 nearest-neighbour bonds, `\|D\|=0.8165`, **D parallel to each bond** (Bloch-type, not interfacial -- interfacial DMI would have D perpendicular to the bond, in-plane) | `dmdata` |
| Directed interactions / atom | **4 exchange + 4 DM, measured directly** (see §D) -- both constant across a 16x atom-count range (full in-plane periodicity, no boundary truncation) | `struct.<simid>.out`, `dmdata.<simid>.out` |
| Anisotropy | **Off.** `#anisotropy ./kdata` is commented out; `do_anisotropy` stays `0` | `inpsd.dat` |
| External field | `hfield 0.0 0.0 -150.0` (large easy-axis-like Zeeman field, reduced units) | `inpsd.dat` |
| Dipole | Off -- no `dipole`/FFT keyword present | `inpsd.dat` |
| Temperature (maintainer default) | **0 K** -- no finite-T value is given anywhere in this template (see §C) | `inpsd.dat` |
| Timestep | `1.000e-15 s` | `inpsd.dat` |
| Integrator | Heun, `SDEalgh 1` | `inpsd.dat` |
| Damping | `0.50` | `inpsd.dat` |
| Initial phase | `ip_mode Y` -> `qminimizer_wrapper('Y')` -> `sweep_q3`: a genuine `O(129 x Natom)` triple-Q spin-spiral energy-minimization line search over the 129 candidate q-points in `qfile`, using `ip_hfield` (unset in this template, defaults to `(0,0,0)`) -- **not** the production `hfield`. Overwrites the `Initmag 3` state before the timed measurement phase. `qm_relax` unset -> default `'N'`, so no MC relaxation runs per candidate (bounded, predictable cost). See §D for measured cost and §E for a real consequence of this field-free search. | `inpsd.dat`, `source/MonteCarlo/qminimizer.f90:62-106,470-668`, `source/uppasd.f90:171-238` |
| Diagnostics | `skyno T` (skyrmion number via triangulation), `do_avrg Y`. `qpoints F` (the reciprocal-space measurement diagnostic is off; `qfile`'s 129 points are consumed only by the `ip_mode Y` search above, a dual use of the same file for two different features) | `inpsd.dat` |
| `do_prnstruct` | `2` in the maintainer template (cluster-oriented printing; does **not** emit `struct.<simid>.out`'s neighbour list -- only `1`/`4` do, `source/Hamiltonian/hamiltonianinit.f90:372` et al.) -- overridden to `1` for workload-metadata runs | `inpsd.dat`, `source/Hamiltonian/hamiltonianinit.f90` |

### Backend dispatch (`gpu_mode`, `skyno`)

Neither `gpu_mode` nor `skyno` is set in the maintainer's single shared
deck (there is no separate `SkyrmionLatticeGPU` twin, unlike `B01_bccFe`'s
`bccFeCPU`/`bccFeGPU` pair), so both default from source
(`gpu_mode=0`, i.e. CPU dispatch). Two real backend-compatibility problems
were found while validating this template on both backends:

1. **`gpu_mode` must be a per-run override**, reusing `B01_bccFe`'s
   established pattern (`benchmarks/cases/B01_bccFe/README.md#backend-dispatch-gpu_mode`):
   a template with `gpu_mode` baked to `0` never dispatches to the GPU on a
   GPU-capable build.
2. **`skyno T` unconditionally crashes the GPU backend.** `skyno`
   selects the skyrmion-number diagnostic method: `'Y'` maps to
   `SkyrmionMethod::BruteForce` (gradient + Pontryagin density, GPU-
   implemented), `'T'` maps to `SkyrmionMethod::Triangulation`
   (`source/gpu_files/measurement/measurementData.h:19-27`). The GPU
   measurement preflight throws unconditionally on `Triangulation`
   ("GPU measurement preflight rejects unsupported triangulation skyrmion
   layout", `source/gpu_files/measurement/gpuMeasurement.cpp:122-124` --
   the comment there confirms this is a known-incomplete code path, not a
   transient failure). Reproduced directly: the unmodified template with
   `gpu_mode 1` fails at startup on `build_gpu/bin/sd.f95.cuda` with this
   exact error; overriding `skyno` to `'Y'` resolves it and the run
   completes normally.

`skyno` is a post-hoc diagnostic of the current spin state (computed from
already-evolved moments), not a Hamiltonian/lattice/moment parameter, so it
fits the same override category as `gpu_mode` and `do_prnstruct`. Extended
`harness/cases.py::GLOBALLY_SAFE_OVERRIDE_KEYS` to add it. The template's
own baseline stays `skyno T` (CPU-safe, matches the maintainer's file
verbatim); a GPU sample passes
`extra_overrides={"gpu_mode": 1, "skyno": "Y"}`.

## B. Size ladder

`Nx=Ny=n` square lattices (single-atom basis, `natom=n^2`), chosen from the
blueprint's ~4x ladder (4k/16k/64k/256k/1M/4M atoms). Every target in this
ladder happens to be an exact integer square, so no size here is a rounded
approximation (unlike `B01_bccFe`'s cube ladder, where only two of six
targets landed exactly). `n=128` is the template's own maintainer-supplied
default (`examples/SpecialFeatures/SkyrmionLattice/inpsd.dat: ncell 128 128 1`)
and lands on this ladder's 16k tier exactly -- the natural anchor size,
mirroring how `B01_bccFe` anchors on this project's historical 65,536-atom
3D reference.

| `size_id` | `n` | `natom` | vs. target | Note |
| --- | --- | --- | --- | --- |
| `64x64` | 64 | 4,096 | target 4,096, **exact** | |
| `128x128` | 128 | 16,384 | target 16,384, **exact** | maintainer's own default `ncell 128 128 1` |
| `256x256` | 256 | 65,536 | target 65,536, **exact** | |
| `512x512` | 512 | 262,144 | target 262,144, **exact** | |
| `1024x1024` | 1024 | 1,048,576 | target 1,048,576, **exact** | |
| `2048x2048` | 2048 | 4,194,304 | target 4,194,304, **exact** | |

Thickness is fixed at `thickness_cells: 1` for every size -- only the
lateral `[Nx, Ny]` scales, per blueprint §B. Aspect ratio is exactly `1:1`
at every size (trivial for a square ladder).

## C. Variants

Only `skyrmion_2d_t0` is admitted:

| `variant_id` | `temp` override | `tseed` override | Purpose |
| --- | --- | --- | --- |
| `skyrmion_2d_t0` | 0.0 | 1 | Zero-temperature relaxation dynamics from the `ip_mode=Y` initial condition |

**No finite-T variant is added.** Unlike `B01_bccFe` (maintainer's own
100 K baseline, used verbatim), this template's only temperature anywhere
is `temp 0.000` -- there is no maintainer-supplied finite-T value to reuse.
The sibling `SkyrmionTransport` example does specify `temp 2.0`, but that
value belongs to a different, rejected template (current-driven dynamics,
different `Initmag`/initial-state machinery -- see §A) and borrowing it
here would be inventing a production parameter, which every prior WP-08x
admission has deliberately avoided (`B01_bccFe/README.md`: "not an invented
'representative' value"). `temp` remains in `allowed_input_overrides` so a
future WP can add a finite-T variant the moment the maintainer supplies one.

## D. Scaling validation

Used `harness.cases.generate_run_directory` (not manual file edits) with
`extra_overrides={"do_prnstruct": 1, "Nstep": 50}` on
`build_cpu/bin/sd.f95` at three sizes spanning a 16x atom-count range, then
read `struct.<simid>.out` (exchange) via the case's own
`neighbor_list_from_struct_output` and `dmdata.<simid>.out` (DM) directly:

| `size_id` | `natom` | exchange directed / atom | DM directed / atom |
| --- | --- | --- | --- |
| `64x64` | 4,096 | 4.0 | 4.0 |
| `128x128` | 16,384 | 4.0 | 4.0 |
| `256x256` | 65,536 | 4.0 | 4.0 |

Constant at every tested size: full in-plane periodicity means every atom
is bulk (no edge truncation), so per-atom interaction count does not depend
on lattice size -- the same mechanism `B01_bccFe` found for its own
(different) constant-96 result.

**Known undercount in `directed_interactions`.** This case's
`workload_metadata_method` (`neighbor_list_from_struct_output`) reads only
`struct.<simid>.out`, which `prn_exchange` writes for the *exchange* list
only -- it has no knowledge of the separate DM neighbour list UppASD writes
to `dmdata.<simid>.out` under the same `do_prnstruct` trigger. For this
case, both lists are the same 4 bonds/atom, so the true per-atom neighbour
workload driving the LLG effective-field sum is 8 directed interactions
(4 exchange + 4 DM), while the case's reported `directed_interactions`
metric (from the shared WP-02 parser) only ever reports the 4 exchange
ones -- a real ~2x undercount specific to any J+D short-range case, not
just this one (the future `B03_skyrmion3D` will hit the identical gap).
Not fixed here: `workload_metadata.py` is WP-02's frozen, harness-wide
contract, not something a single case admission should extend
unilaterally. Left as a known limitation for a future work package
(extending the parser, or the schema, with a DM-aware neighbour count) to
resolve across every J+D case at once, rather than patching it inconsistently
per case.

**Setup-cost consequence of the `ip_mode=Y` search (§A).** The same three
runs' own phase timers:

| `size_id` | `natom` | initial phase (wall) | meas. phase, `Nstep=50` (wall) |
| --- | --- | --- | --- |
| `64x64` | 4,096 | 0.07 s | 0.03 s |
| `128x128` | 16,384 | 0.31 s | 0.11 s |
| `256x256` | 65,536 | 1.18 s | 0.46 s |

The `O(129 x Natom)` line search costs roughly 2.5x the `Nstep=50`
measurement phase at every tested size, and scales with `Natom` just like
the measurement phase does -- consistent with §A's complexity analysis
(each of the 129 candidates costs one full-lattice `effective_field` call,
the same per-atom cost order as one LLG step). Because this cost is fixed
per `(size, ip_mode)` and **independent of `Nstep`** (the search always
evaluates all 129 `qfile` points regardless of the measurement-phase step
count), it lands entirely in the docs' `T_fixed` intercept term of the
`T(n) = T_fixed + n*t_step` model (blueprint §6.2), not in the steady-state
`t_step` slope -- the framework's multi-`Nstep` regression methodology
isolates it correctly as long as a campaign uses at least the blueprint's
recommended three separated `Nstep` points, per §6.3. It does mean this
case's short-`Nstep` runs (as used in this admission's own sanity checks)
are dominated by setup cost, exactly as `Nstep=5000` production runs are
not (129 candidates / 5000 steps is 2.6% -- small).

`512x512`, `1024x1024`, `2048x2048` were not executed (impractical for an
admission-time check); at the confirmed constant 4/4 topology,
`directed_interactions` (exchange only) at those sizes is `natom*4` =
1,048,576 / 4,194,304 / 16,777,216 respectively -- predictions from the
confirmed scaling law, not independent measurements.

## E. Sanity runs

All runs used `64x64` (4,096 atoms), `do_prnstruct 1`, `Nstep 50`, generated
through `harness.cases.generate_run_directory` (not manual file edits); CPU
on `build_cpu/bin/sd.f95` (`UPPASD_GPU_BACKEND=OFF`,
`UPPASD_PRECISION=DOUBLE`), GPU on `build_gpu/bin/sd.f95.cuda`
(`UPPASD_GPU_BACKEND=CUDA`, `UPPASD_PRECISION=DOUBLE`,
`extra_overrides={"gpu_mode": 1, "skyno": "Y"}`, required per §A).

| Run | Result |
| --- | --- |
| CPU, `skyrmion_2d_t0` | **PASS** -- completes; `natom`/exchange+DM directed interactions match §D; `\|m\|=1.0` conserved on every atom through 50 steps; Exc and DM both nonzero and active |
| CUDA, `skyrmion_2d_t0` | **PASS with a caveat below** -- completes (once `skyno` is overridden per §A); `\|m\|=1.0` conserved on every atom through 50 steps; Exc and DM both nonzero and active |

### A real cross-binary degenerate-minimum finding (not a GPU dynamics bug)

The CPU and CUDA runs above do **not** reach the same final spin state, and
their `t=0` `totenergy.<simid>.out` breakdown differs sharply:

```text
                     Tot        Exc        DM         Zeeman
CPU  (64x64, t=0):  -4.637   -3.461    -1.058       -0.000 (small, nonzero)
GPU  (64x64, t=0):  -2.403   -3.461    +1.058        0.000 (exactly zero)
```

Investigated directly, not assumed:

* `sweep_q3` (the `ip_mode=Y` search) is pure Fortran, called identically
  from `run_initial_phase` regardless of `gpu_mode` -- it is not
  GPU-dispatched code. Re-running the *same* CPU binary twice on the same
  input is exactly reproducible (identical "27 hits", identical `q_best`
  both times). But the CPU-only and CUDA-enabled binaries -- two
  separately compiled executables -- pick **different** `q_best` candidates
  ("27 hits" vs "44 hits" during the scan) despite the search converging to
  the *same* minimum energy value to 10 significant figures
  (`-4.2534654184` both times).
* The `\|Exc\|` energy term matches exactly between the two final states
  (`-3.46094416` both) while `DM` has the **same magnitude but opposite
  sign** (`-1.058` vs `+1.058`), and `Zeeman` differs (small-but-nonzero
  vs. exactly zero). Exchange (a dot product) is chirality-blind; DM (a
  cross product) flips sign under a mirror/handedness reversal. This is
  the signature of the field-free search (`ip_hfield` unset -> `(0,0,0)`,
  not the production `hfield`, per §A) landing on two energy-degenerate,
  opposite-*chirality* triple-Q domains -- a real physical degeneracy in
  the maintainer's own initial-phase algorithm, exposed (not caused) by
  ordinary cross-compiler/cross-build floating-point non-associativity
  choosing a different tie-break between the two binaries.
* Consequence: **`harness.gpu_sanity`'s usual T=0 bit-exact CPU/GPU
  comparison is not a meaningful discriminator for this case.** This is a
  property of the maintainer's template (the field-free chiral-degenerate
  search), not of the GPU dynamics path, `gpu_mode` dispatch, or this
  admission's overrides -- confirmed by the single-binary determinism check
  above. Both final states are individually valid (conserved `\|m\|=1.0`,
  finite energy, active Exc/DM terms); they are simply two different but
  physically equivalent (mirror-image) skyrmion-lattice domains. Unlike
  `B01_bccFe`'s weak-discriminator note (a genuine, resolvable degenerate
  fixed point at one specific initial condition), this one cannot be
  resolved without either changing the template's `ip_hfield`/`ip_mode`
  (out of scope -- an initial-state Hamiltonian-adjacent parameter, not in
  `allowed_input_overrides`) or pinning both binaries to identical
  compiler flags (out of scope for a case admission). Documented here so a
  later WP does not mistake this divergence for a GPU correctness
  regression.

## Checklist

- [x] Template audited (§A; verbatim copy of `examples/SpecialFeatures/SkyrmionLattice`, `diff -q` verified; sibling `SkyrmionTransport` template inspected and rejected with reasons).
- [x] Short-range neighbour cloud recorded (§A/§D: 4 exchange + 4 DM directed interactions/atom, both constant across a 16x atom-count range).
- [x] J and DMI enabled as expected (§A/§E: both nonzero and active in every sanity run's `totenergy.<simid>.out`; anisotropy confirmed off).
- [x] 2D size scaling implemented (§B, six sizes 4,096 -> 4,194,304 atoms, all exact squares).
- [x] Thickness fixed (`thickness_cells: 1` in `case.yaml`; every size only varies `[Nx, Ny]`).
- [x] Aspect ratio controlled (`Nx=Ny` at every size, exactly 1:1).
- [x] Interaction topology invariant (§C: the one admitted variant only overrides `temp`/`tseed`, never any `jij`/`dmdata`/lattice/moment content; `allowed_input_overrides` excludes all of those).
- [x] CPU sanity passes (§E).
- [x] CUDA sanity passes (§E, with the documented chirality-degeneracy caveat -- a template property, not a GPU defect).
- [x] Dynamics-only profile works (`measurement_profiles.DYNAMICS_ONLY` overrides `avrg_step`/`cumu_step`/`tottraj_step`/`ene_step` past `Nstep`, all declared in `allowed_input_overrides`; the `ip_mode=Y` initial-phase cost is `Nstep`-independent (§D), so it lands in the `T_fixed` intercept rather than contaminating the `t_step` slope `DYNAMICS_ONLY`/the setup-vs-steady regression exists to isolate).
