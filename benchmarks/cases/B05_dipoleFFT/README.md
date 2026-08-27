# B05_dipoleFFT — open-boundary simple-cubic dipole/FFT workload

Admission record for WP-08E (blueprint §3 "B05 — dipole/FFT workload"). This
case is the framework's algorithmically different production benchmark:
unlike `B01_bccFe`-`B04_dhcpNd` (all neighbour-list workloads), its
dominant per-step cost is a global FFT-based dipole/magnetostatic
calculation, not a local neighbour sum. Per the blueprint, this case
deliberately does **not** use `directed_interaction_visits/s` as its
principal metric (§H below).

## A. Template audit

### No maintainer example deck exists; this project's own dipole-validation base is the closest approved source

`examples/` has exactly one deck naming `do_dip`
(`PhaseDiagrams/Ising-Tsweep/Base`), and it sets `do_dip 0` (dipole
disabled) on a Monte Carlo (`mode I`) system -- not usable.

This project's own **internal** dipole development/validation pipeline
(`tests/dipole_validation/generate_cases.py`) builds every one of its own
dipole regression and WP10 acceptance cases from two fixtures,
`tests/SC2d_dipole` and `tests/SC3d_dipole` (identical `jij`/`dmdata`/
`kdata`/`momfile`/`posfile`, differing only in `inpsd.dat`'s `ncell`/`BC`)
-- confirmed by reading `generate_cases.py` directly:
`sc3d = ROOT / "tests/SC3d_dipole"`, `sc2d = ROOT / "tests/SC2d_dipole"`.
These fixtures are untracked in this branch's git history (`git ls-files`
returns nothing for either path) but are the real, load-bearing base every
tracked dipole acceptance script in `tests/dipole_validation/` depends on --
this admission uses the same base, not a fresh invention.

**Two disclosed, minimal changes from that base**, both documented in full
rather than left implicit:

| File | Base fixture value | This case's value | Rationale |
| --- | --- | --- | --- |
| `jij` | 4 in-plane bonds (`+-x`, `+-y`), `J=-0.00000` (literally zero -- confirmed by reading the raw file bytes) | 6 isotropic simple-cubic bonds (`+-x, +-y, +-z`), `J=1.0` | The base fixture's `J=0` is deliberate *for its own purpose* -- isolating the dipole field from exchange in a validation test (`tests/dipole_validation/README.md`: independent point-dipole analytic references, not ASD dynamics). Blueprint B05 section A requires "ordinary exchange remains present enough for valid ASD dynamics," so this admission restores a real `J`, reusing the exact magnitude convention `B02_skyrmion2D`/`B03_skyrmion3D` already established for the same bond topology (`J=1.0`). The missing `+-z` bonds in the un-modified base (present even in `SC3d_dipole`'s own "3D" `ncell 16 16 16` variant) would leave every `z`-layer of a genuinely 3D shape exchange-decoupled from its neighbours -- not physically sensible for this case's 3D-shaped ladder tiers, so `+-z` bonds are added, matching `B03_skyrmion3D`'s own 6-neighbour isotropic convention exactly. |
| `momfile` | direction `(0, 0, 1)` | direction `(0.1, 0.0, 1.0)` | Reused verbatim from `B03_skyrmion3D`'s own established non-trivial-direction convention (`B03_skyrmion3D/README.md`), needed here for an unrelated but analogous reason: see §B. |

`dmdata`/`kdata` (present in the source fixtures but unreferenced by
`inpsd.dat` there too) were not copied -- this case's Hamiltonian is
exchange + dipole only, matching blueprint B05 section A's "unrelated
complicated Hamiltonian terms do not dominate."

### `do_dip` (CPU) and `gpu_dipole_mode` (GPU) are two independent, unrelated keywords

`source/Input/inputhandler.f90` parses them separately:

```text
do_dip            0=off, 1=finite direct sum, 2=finite macrocell, 3=finite zero-padded FFT (legacy CPU, dipolemanager.f90)
gpu_dipole_mode   OFF / EWALD3D_FFT / OPEN_FFT (WP10 CV6, GPU-only, source/gpu_files/...)
```

Selecting a GPU mode requires `do_dip=0` -- confirmed directly: a run with
`do_dip=1` and `gpu_dipole_mode=OPEN_FFT` together
(`extra_overrides={"gpu_dipole_mode": "OPEN_FFT"}` without also overriding
`do_dip` back to `0`) fails hard at startup:

```text
ERROR STOP OPEN_FFT requires do_dip=0; legacy and GPU dipoles cannot be combined
```

exactly the WP10 blueprint's own documented startup validation gate
(section 8: "unsupported ... irrelevant Ewald inputs fail before
allocation"), reproduced here for real on this admission's own generated
run.

### Legacy CPU `do_dip=3` (the CPU FFT dipole) is confirmed broken; this case never uses it

`docs/WP10_OPEN_FFT_BLUEPRINT.md` section 9.2 ("Confirmed defects") records
a real, measured defect in the legacy CPU FFTW dipole path: wrong kernel
temporary dimensions and an `I1`/`I3` index-convention mismatch between
packing and unpacking, producing a **30.59% maximum field error** and
**-8.88% dipole energy error** versus the validated `do_dip=1` reference on
a real `3x2x1` test case (dated 2026-07-28). WP11 (the CPU dipole repair
project) has not started -- its checklist in the same document is entirely
`[ ]`. Per `tests/dipole_validation/README.md`'s own validation-status
section, only `do_dip=1` (finite direct point-dipole sum) has passed
independent analytic-reference validation; `do_dip=2`/`do_dip=3` explicitly
have not.

**Consequence for this case**: the CPU dipole arm uses `do_dip=1`
(validated, correct, but `O(Natom^2)` -- no FFT), never `do_dip=3`. This
means this case's CPU/GPU comparison is a genuine **algorithmic** crossover
(`O(Natom^2)` direct sum vs. `O(NFFT log NFFT)` FFT convolution), not merely
a same-algorithm hardware comparison -- there is no correct
same-algorithm CPU counterpart to compare against today. This is disclosed
explicitly rather than silently presented as an apples-to-apples FFT-vs-FFT
comparison.

### The GPU `OPEN_FFT` path is real, accepted, production dispatch

Unlike the CPU FFT path, GPU `gpu_dipole_mode OPEN_FFT` is accepted
production code: `docs/luna_wp10_final_acceptance.md` records "GO -- WP10.9
is complete for the accepted CUDA scope" (evaluated revision `07d14cc`);
commit `378010c3` ("Close WP10 OPEN_FFT acceptance") is confirmed an
ancestor of this branch's current `HEAD`
(`git merge-base --is-ancestor 378010c3 HEAD`). The accepted scope is BC
`0 0 0`, block `(1,1,1)`, `NA=1`, CUDA fp64 (accepted) and fp32
(accepted under a documented `5e-5` mixed-precision budget); HIP is absent
and non-gating. This case's template (single-atom basis, open boundary,
default block size) sits entirely inside that accepted scope.

### A required-together override trap, found and disclosed

A GPU dipole sample needs **three** overrides set together:
`gpu_mode=1` (backend dispatch), `do_dip=0` (mandatory whenever
`gpu_dipole_mode` is not `OFF`), and `gpu_dipole_mode=OPEN_FFT`. Verified
directly what happens if `gpu_mode` is left at this template's CPU-dispatch
baseline (`0`) while the other two are set correctly
(`extra_overrides={"do_dip": 0, "gpu_dipole_mode": "OPEN_FFT"}`, no
`gpu_mode` override): the run **completes successfully** (`exit=0`, no
error, no warning, no `OPEN_FFT` startup diagnostic line at all) but
silently reports `Dip=0.00000000E+00` -- the dipole term was never
constructed because `do_gpu` never became `'Y'`, so the GPU driver (and
therefore `gpu_dipole_mode`) was never reached; `do_dip=0` disabled the CPU
path too. This is exactly blueprint B05 section F's warning ("the benchmark
does not accidentally fall back to another backend") made concrete: any
future campaign generating GPU dipole samples for this case must pass all
three overrides together, or it will silently collect a real, schema-valid,
*wrong* record (dipole physically absent, but not flagged `FAILED`).

## B. A real, shape-dependent dipole-cancellation finding (not a bug)

An initial sanity probe at an **isotropic cube** (`4x4x4`, `Initmag 3`,
uniform direction) produced `Dip=1.15e-22` -- indistinguishable from
floating-point noise, next to `Exc=-4.5`. This is a real physical effect,
not a defect: a simple-cubic lattice of parallel point dipoles has a net
local dipolar field that cancels by cubic symmetry (the classical
Lorentz-cavity result for a cubic Bravais lattice) -- a symmetric finite
cube approximates this cancellation closely. Re-running the identical
system at a **thin-film** shape (`16x16x1`) instead gives
`Dip=1.40498391E-05` -- fourteen orders of magnitude larger, because a thin
film's demagnetizing tensor is strongly shape-anisotropic (real, textbook
magnetostatics -- shape anisotropy is the dominant magnetic-thin-film
effect). This is the same "uniform/symmetric state is a degenerate special
point" theme `B01_bccFe`/`B03_skyrmion3D` already documented for their own
Hamiltonians, here driven by lattice/sample *shape* rather than initial
spin texture.

**Consequence**: this case's size ladder (§C) is deliberately thin-film-
and racetrack-shaped, not isotropic cubes -- both because that is the
physically sensible geometry for a magnetostatics-dominated benchmark
(blueprint B05 section A) and because a cube would make the dipole term's
own sanity check ("output is numerically sane," blueprint B05 section F)
indistinguishable from a bug. A small number of cube/near-cube tiers are
kept specifically to document the contrast, not as the ladder's primary
shape.

## C. Size ladder

`[Nx, Ny, Nz]`, single-atom basis (`natom = Nx*Ny*Nz`), `NFFT` = the exact
zero-padded embedding every finite/open FFT dipole path (legacy CPU
`do_dip=3`'s intended design, and the accepted GPU `OPEN_FFT`) uses:
`NFFT_i = 2*Nx_i - 1` per axis (WP10 blueprint section 8: "`P_i >= 2G_i-1`
and cleared padding"; `harness/workload_metadata.py::fft_grid_from_replication`
implements exactly this formula and needed no change for this case).
Verified directly against this case's own GPU startup diagnostic at every
tested size (§D): predicted and reported `NFFT` match exactly.

| `size_id` | `[Nx,Ny,Nz]` | `natom` | `NFFT` (padded) | Shape class |
| --- | --- | --- | --- | --- |
| `thin_16x16x1` | `[16,16,1]` | 256 | `31x31x1` = 961 | thin film |
| `thin_32x32x1` | `[32,32,1]` | 1,024 | `63x63x1` = 3,969 | thin film |
| `racetrack_256x16x1` | `[256,16,1]` | 4,096 | `511x31x1` = 15,841 | racetrack (matches `docs/terra_wp10_8_open_performance.md`'s own racetrack sample) |
| `thin_64x64x1` | `[64,64,1]` | 4,096 | `127x127x1` = 16,129 | thin film (same `natom` as the racetrack above -- same atom count, very different `NFFT`, see §H) |
| `thin_128x64x1` | `[128,64,1]` | 8,192 | `255x127x1` = 32,385 | thin film (matches `docs/terra_wp10_8_open_performance.md`'s own thin-film sample) |
| `racetrack_512x32x1` | `[512,32,1]` | 16,384 | `1023x63x1` = 64,449 | racetrack |
| `cube_8x8x8` | `[8,8,8]` | 512 | `15x15x15` = 3,375 | cube (§B contrast case) |
| `cube_16x16x16` | `[16,16,16]` | 4,096 | `31x31x31` = 29,791 | cube (§B contrast case; same `natom` as `thin_64x64x1`/`racetrack_256x16x1`) |
| `cube_32x24x16` | `[32,24,16]` | 12,288 | `63x47x31` = 91,791 | near-cube (matches `docs/terra_wp10_8_open_performance.md`'s own "fully 3D" sample) |

Three tiers deliberately match shapes already measured at the standalone
FFT-kernel level in WP10.8's own performance report, so this admission's
production-executable timings are directly comparable (not identical --
different measurement methodology, see §D) to that prior, narrower
measurement. No memory-limit classification (`gpu_memory.py`/
`host_memory.py`, as `B04_dhcpNd` needed) is required anywhere on this
ladder: real measured GPU dipole-phase peak allocation tops out at 22.5 MB
(§D), utterly negligible against this host's 16 GB GPUs.

## D. FFT/dipole metadata and real measurements

All numbers below are real: `harness.cases.generate_run_directory`
(never manual file edits), CPU on `build_cpu/bin/sd.f95`, GPU on
`build_gpu_wp10_cuda/bin/sd.f95.cuda` (the WP10-accepted CUDA fp64 build;
`build_gpu_fp32/bin/sd.f95.cuda` for the fp32 check), `Nstep=50`,
`dipole_on` variant (`extra_overrides={"gpu_mode": 1, "do_dip": 0,
"gpu_dipole_mode": "OPEN_FFT"}` for GPU samples).

| `size_id` | `natom` | `NFFT` | GPU startup diagnostic | Dip energy (fp64, mRy/atom) |
| --- | --- | --- | --- | --- |
| `thin_16x16x1` | 256 | 961 | `16 x 16 x 1 active cells, 1 basis channel; padded FFT grid 31 x 31 x 1, 496 half-spectrum points` | `1.40498391E-05` |
| `racetrack_256x16x1` | 4,096 | 15,841 | `256 x 16 x 1 active cells ...; padded FFT grid 511 x 31 x 1, 7936 half-spectrum points` | `1.53331072E-05` |
| `thin_128x64x1` | 8,192 | 32,385 | `128 x 64 x 1 active cells ...; padded FFT grid 255 x 127 x 1, 16256 half-spectrum points` | `1.65936711E-05` |
| `cube_32x24x16` | 12,288 | 91,791 | `32 x 24 x 16 active cells ...; padded FFT grid 63 x 47 x 31, 46624 half-spectrum points` | `3.10078585E-06` (smaller -- residual shape anisotropy only, §B) |

`natom` and `NFFT` (both active and padded) match this case's own §C
prediction exactly at every tested size -- confirmed, not assumed. Every
diagnostic line also names the dipole implementation identifier this
admission's blueprint section G asks for: `OPEN_FFT`, `1 basis channel`,
`block 1 x 1 x 1`, `1 ensemble` -- printed by
`source/gpu_files/...` at GPU dipole setup, not inferred.

Real measured GPU memory (bytes, from this case's own runs -- the
production binary's own live tensor-allocation tracker, not an estimate):

| `size_id` | persistent | construction | workspace | dipole-phase peak |
| --- | --- | --- | --- | --- |
| `thin_16x16x1` | 165,320 | 69,192 | 0 | 240,656 |
| `racetrack_256x16x1` | 2,665,160 | 1,140,552 | 2,537,040 | 6,441,056 |
| `thin_128x64x1` | 5,456,072 | 2,331,720 | 5,182,544 | 13,166,944 |
| `cube_32x24x16` | 15,595,880 | 6,608,952 | 0 | 22,499,744 |

### Difference from WP10.8's own performance numbers

`docs/terra_wp10_8_open_performance.md` measured these same three shapes
(`thin_128x64x1`, `racetrack_256x16x1`/`512x32x1`, `cube_32x24x16`... its
own `32x24x16`) using a **standalone kernel microbenchmark**
(`dipole_gpu_fft_benchmark --open --grid ... --fft-grid ... --warmup 10
--iterations 100`), deliberately isolating `pack/clear -> forward FFT ->
contract -> inverse FFT -> scatter -> energy` from the rest of an ASD
simulation. This admission instead times the complete production
executable end to end (this project's core principle since WP-01/WP-02:
"never a kernel microbenchmark"), so this case's own future performance
campaign numbers will include ordinary-ASD overhead (moment
update/normalization, measurement I/O, the -- deliberately minimal, §A --
exchange sum) on top of the dipole kernel WP10.8 isolated. The shapes are
reused specifically so a future analysis *can* compare the two
measurement layers directly; this admission does not claim they will
match.

## E. Setup versus steady-state (blueprint B05 section D)

This case introduces no new timing mechanism: `harness.steady_state`'s
existing `T(n) = T_fixed + n*t_step` multi-`Nstep` regression
(`FitResult.setup_or_fixed_intercept`, "deliberately named for what it is
-- a fitted intercept -- not assumed to be pure startup time") already
applies unchanged.

**What the intercept contains for this case, demonstrated, not assumed**:
UppASD's own phase report (printed by every run, §D and earlier sections)
separates `Time for initialization` from `Time for meas. phase` --
initialization (before the `Nstep` loop starts) is where Hamiltonian setup,
GPU device-tensor allocation, and this case's `OPEN_FFT` plan
creation/kernel upload all happen (confirmed: the `Gpu: OPEN_FFT ...
staged`/`... allocation ...`/`... operator ready` diagnostic lines in §D
are printed during this phase, before the first measured iteration).
Measured initialization time on this case's own runs: `0.01 s` at
`natom=256`, growing to `0.05 s` at `natom=12,288` (§D's four sizes,
`Nstep=50` sanity runs). A multi-`Nstep` regression's fitted intercept will
therefore contain real FFT-plan/kernel-upload cost **plus** ordinary
process/Hamiltonian startup cost that has nothing to do with FFT planning
(reading `inpsd.dat`, allocating moment arrays, building the trivial
6-neighbour exchange list) -- consistent with `steady_state.py`'s own
"not assumed to be pure startup time" caution, and this case is the first
admitted case where "FFT planning" is even a plausible intercept
component worth naming explicitly. A future performance campaign wanting
to isolate the FFT-plan-specific fraction of the intercept would need to
difference this case's `dipole_on` vs. `dipole_off` intercepts at matched
`natom` (§F) -- not attempted here, this is a case-admission record, not a
campaign.

## F. Variants

| `variant_id` | CPU override | GPU overrides (`extra_overrides`) | Purpose |
| --- | --- | --- | --- |
| `dipole_on` | `do_dip=1` | `do_dip=0`, `gpu_dipole_mode=OPEN_FFT`, `gpu_mode=1` | Dipole Hamiltonian active |
| `dipole_off` | `do_dip=0` | `gpu_dipole_mode=OFF` (already the template default), `gpu_mode=1` | Same exchange-only system, dipole disabled -- a physically legitimate control (blueprint B05 section E), not a synthetic kernel-only stand-in: it is the same lattice, same exchange, same integrator, only the dipole term removed |

Both variants were run and validated (§G). `dipole_off` is admitted because
disabling the dipole flag genuinely leaves an otherwise valid, directly
comparable ASD benchmark -- exactly blueprint B05 section E's admission bar
-- letting a future campaign estimate dipole's incremental production cost
directly (`dipole_on` minus `dipole_off` at matched size).

## G. CPU/GPU sanity

All runs used `thin_16x16x1` (256 atoms), `Nstep=50` (`dipole_on`) or
`Nstep=10` (rejection/trap checks, §A), generated through
`harness.cases.generate_run_directory`.

| Run | Result |
| --- | --- |
| CPU, `dipole_on` (`do_dip=1`) | **PASS** -- completes; `natom=256`; `Dip=1.40498391E-05` mRy/atom |
| CPU, `dipole_off` (`do_dip=0`) | **PASS** -- completes; `Dip=0.00000000E+00` exactly, `Tot=Exc` exactly |
| CUDA fp64, `dipole_on` (`OPEN_FFT`) | **PASS** -- completes; `Dip=1.40498391E-05` mRy/atom -- **bit-identical to the CPU `do_dip=1` value to all 9 printed significant figures**, the strongest possible confirmation that the GPU FFT-convolution dipole and the CPU direct-sum dipole agree, despite being algorithmically unrelated implementations (§A) |
| CUDA fp64, `dipole_off` | **PASS** -- completes; `Dip=0.00000000E+00`, no `OPEN_FFT` diagnostic line printed (confirms the control genuinely disables dipole on GPU, not merely on CPU) |
| CUDA fp32, `dipole_on` (`OPEN_FFT`) | **PASS** -- completes; `Dip=1.40498369E-05` -- matches fp64/CPU to 6 significant figures (relative deviation `1.6e-7`), consistent with WP10.8's own accepted CUDA fp32 error budget (`<5e-5`) |
| CUDA fp64, misconfigured (`do_dip=1` + `gpu_dipole_mode=OPEN_FFT`, no override back to `do_dip=0`) | **Correctly rejected** -- `ERROR STOP OPEN_FFT requires do_dip=0; legacy and GPU dipoles cannot be combined`, exit code 1 (§A) |
| CUDA fp64, misconfigured (`gpu_mode` left at CPU-dispatch baseline, `do_dip=0` + `gpu_dipole_mode=OPEN_FFT` set) | **Silent no-op, not a crash** -- exit 0, `Dip=0.00000000E+00`, no `OPEN_FFT` diagnostic printed at all (§A's disclosed override-ordering trap) |

`Mean-field estimate of Tc: 210.5 K` (printed by every run using this
case's `J=1.0` exchange) confirms the restored exchange term is genuinely
active, not inert -- a direct numerical-sanity check, not an assumption.

## H. Analysis expectation (blueprint B05 section H)

This case must not use `directed_interaction_visits/s` as its principal
normalized metric -- its Hamiltonian's neighbour-list term is a
deliberately minimal, non-dominant 6-bond-per-atom exchange sum (§A), and
its dominant per-step cost is the FFT dipole convolution, not neighbour
traversal. Principal metrics for a future performance campaign on this
case should be `atom-steps/s` (`natom * Nstep / steady_step_seconds`,
already how every prior case's throughput is computed) and
**`FFT-grid-points/step/s`** (`NFFT * Nstep / steady_step_seconds`, using
this case's own `fft_grid_points` workload metadata, §C/§D) -- the second
metric has no precedent in `B01_bccFe`-`B04_dhcpNd`, all of which are pure
neighbour-list cases with no `NFFT` concept. `thin_64x64x1` and
`racetrack_256x16x1` (identical `natom=4,096`, `NFFT` 16,129 vs. 15,841 --
close by design) versus `cube_16x16x16` (also `natom=4,096`, but
`NFFT=29,791`, nearly double) is this ladder's built-in demonstration that
`natom` alone does not predict FFT cost, and that shape -- not just size --
determines it (blueprint B05 section C: "The analysis must be able to plot
runtime versus both").

## Checklist

- [x] Production dipole template audited (§A: no maintainer `examples/`
  deck exists; this project's own internal dipole-validation base
  (`tests/SC2d_dipole`/`SC3d_dipole`, referenced by
  `tests/dipole_validation/generate_cases.py`) used instead, two disclosed
  minimal changes, every other file/value verbatim).
- [x] Dipole implementation identified (§A/§D: CPU `do_dip=1` finite direct
  sum, validated; GPU `gpu_dipole_mode OPEN_FFT`, WP10-accepted; legacy CPU
  `do_dip=3` explicitly identified as confirmed-broken and excluded).
- [x] FFT grid metadata captured (§D: `fft_grid`/`fft_grid_padded`/
  `fft_grid_points` via the existing `fft_grid_from_replication` parser,
  cross-checked against this case's own real GPU startup diagnostics at
  every tested size).
- [x] Geometry-preserving size ladder exists (§C: nine sizes, thin-film/
  racetrack-dominant by design, with cube contrast tiers; three tiers match
  WP10.8's own already-measured shapes).
- [x] Natom recorded (§C/§D, every size).
- [x] NFFT recorded (§C/§D, every size, both active and zero-padded).
- [x] Setup/intercept separated from steady-state timing (§E: reuses
  `harness.steady_state` unchanged; documents what this case's intercept
  demonstrably contains, per blueprint B05 section D's explicit caution).
- [x] `DIPOLE_ON` validated (§F/§G).
- [x] Optional `DIPOLE_OFF` control added -- verified physically legitimate
  (§F: identical exchange system, dipole term removed, not a synthetic
  kernel-only stand-in) and passes on both backends (§G).
- [x] CPU sanity passes (§G).
- [x] CUDA sanity passes (§G; fp64 and fp32; plus two real misconfiguration
  probes -- one correctly rejected, one a disclosed silent-no-op trap).
- [x] No neighbour-count throughput metric misapplied (§H: `workload_class:
  fft_dipole` -> `FFT_DIPOLE`, not `NEIGHBOR_LIST`;
  `directed_interaction_visits/s` explicitly named as the metric this case
  must not use).
- [x] FFT-specific analysis metadata integrated (§D/§H: `fft_grid_points`
  and the `atom-steps/s`/`FFT-grid-points/step/s` metric pair).
