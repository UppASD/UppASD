# RCG-10: release evidence reconciliation and parent handback

**Task:** RCG-10, "Reconcile release evidence with the parent blueprint"
**Parent specification:** `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`
**Remediation overlay:** `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md` §9
**Date:** 2026-08-15
**Status: CLOSED (2026-08-15, Human decision: Anders Bergman).** The branch is
accepted as a **correctness/reference implementation, not a
production-performance release**. Four Human decisions on 2026-08-15 closed the
items this task could not settle on its own (§9), the last of them approving the
release wording and the §14 item 9 amendment as written. This completes the
remediation program and hands authority back to the parent blueprint. One gap is
carried forward unmet rather than reworded: HIP release evidence, for want of
hardware.

---

## 1. Summary

RCG-10 is the only authorized handback to the parent blueprint (remediation
§14). It does not implement physics. It re-validates the candidate commit from
a clean clone, audits every parent acceptance box that the remediation
provisionally reopened, completes the parent's §12 cross-task review
checklists from current evidence, and removes stale or superseded claims.

Three things changed materially in the parent blueprint as a result:

1. **The `1.30x` speedup claim is withdrawn.** It compared the adaptive runtime
   at a low fine fraction against the adaptive runtime at 100% fine — two modes
   of the same object, never a comparison with UppASD's atomistic GPU path
   (finding F-10). Parent CG-10's "A real active-DOF crossover is measured" box
   is unchecked, and the GPU backend is described as a correctness prototype.
2. **Parent §14 item 9 is explicitly amended,** not quietly reinterpreted. The
   remediation blueprint's §13 provides for exactly this outcome: items 1–7 and
   10–12 pass, item 9 fails, and the release is a correctness/reference
   implementation.
3. **The finite-temperature boundary is made explicit in both directions**
   (§9.1). The audit found that nothing gated `Temp > 0` against coarse blocks;
   the Human owner resolved it by accepting the combination — finite
   temperature at the atomistic level, coarse blocks at T=0 — and this session
   paid the resulting evidence debt with the `finite_temperature_mixed`
   fixture and its executed negative control. Parent §14 item 10's
   finite-temperature clause is amended accordingly.

The audit also closed two gaps that earlier tasks left open rather than
inferring them shut: the RCG-09A.4 staged ASD parity harness was executed on a
real CUDA device for the first time (§4.3), and the three RCG-09A negative
controls that had only ever been compiled were executed and each failed at the
stage its own table predicted (§4.4).

---

## 2. Candidate and clean-clone provenance

Validation began from a fresh `git clone` of the branch, not from the working
tree, and not from any existing build directory (remediation §2.2).

```text
commit        abda6534ee6c60e438f2c72e6f63863a23dace62
subject       RCG-09A: complete adaptive all-fine parity coverage
branch        gpu_hip_cu_ab_cg
describe      v6.0.2-488-gabda
clone         git clone --no-local --single-branch --branch gpu_hip_cu_ab_cg
git status    0 entries (clean, no untracked and no modified files)
```

The source worktree was also confirmed clean of *tracked* modifications before
the clone. Its untracked entries are local build directories, run outputs, and
example artifacts; none is referenced by any harness. `adaptive-cg-fixture-
dependencies` (the `PKG-TRACKED-E2E` audit) passes on the clean clone, which is
the executable form of that claim: every input, interaction file, mask, and
restart the e2e harnesses consume comes from tracked source.

**Hygiene note, not a blocker.** Running the CPU CG-13 suite modifies three
*tracked* files in place — `examples/AdaptiveCoarseGraining/{adaptive,
initial_phase_texture,static_mixed}/uppasd.adaptive.yaml` — because those
committed run records are overwritten with the current date and `git_revision`.
Acceptance does not depend on them, but a validation run leaves the tree dirty.

### Environment

```text
OS          Linux 6.8.0-137-generic x86_64
CPU         11th Gen Intel Core i9-11900 @ 2.50GHz (2 cores available)
Fortran     GNU Fortran 13.3.0
C/C++       GCC 12.4.0
CUDA        nvcc 13.3, V13.3.73; CMAKE_CUDA_ARCHITECTURES=native
GPU         2x NVIDIA RTX A4000, driver 610.57.04
HIP/ROCm    absent: no hipcc on PATH, no /opt/rocm
Build type  Release
```

---

## 3. Release matrix

| Backend / precision | Configure + build | Suite | Result |
|---|---|---|---|
| CPU fp64 | clean out-of-tree, succeeded | `ctest -L cg13-cpu` | **29/29 passed**, 38.28 s |
| CPU fp64 | same build | full `ctest` (38 tests) | 37/38; one **pre-existing** failure, §3.2 |
| CPU fp64 | same build | `regression-test`, `asd-tests` | both **passed** |
| CUDA fp64 | clean out-of-tree, succeeded | `ctest -L cg13-cuda` | **32/32 passed**, 152.98 s |
| CUDA fp32 | clean out-of-tree, succeeded | `ctest -L cg13-cuda` | **31/32 passed**, 153.04 s; one known pre-existing failure, §3.3 |
| HIP any | **not attempted** | — | **UNAVAILABLE**: no toolchain, no device. Not inferred from CUDA. |
| Legacy GNU make | **fails** | — | **Not a supported build** for current GPU builds (Human decision 2026-08-15). Diagnosis in §3.4; not a release gate. |

Commands, verbatim:

The matrix below was re-run after this task added the `finite_temperature_mixed`
fixture (§9.1); the totals include it.

```sh
cmake -S clone -B build_cpu       -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=OFF  -DUPPASD_PRECISION=DOUBLE
cmake -S clone -B build_cuda      -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=DOUBLE
cmake -S clone -B build_cuda_fp32 -DCMAKE_BUILD_TYPE=Release -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=SINGLE
ctest --test-dir build_cpu       -L cg13-cpu  --output-on-failure
ctest --test-dir build_cpu                    --output-on-failure
ctest --test-dir build_cuda      -L cg13-cuda --output-on-failure
ctest --test-dir build_cuda_fp32 -L cg13-cuda --output-on-failure
```

### 3.1 Backend availability is stated, never inferred

HIP is unavailable in every environment used by this remediation program, from
RCG-02 onward. `cg13-hip` labels are added only by a HIP backend build, so a
CUDA build cannot report HIP validation. No HIP result is claimed, projected,
or inherited from CUDA anywhere in the reconciled documents. This is tracked as
the still-open RCG-04-FU1 / RCG-05-FU1 / RCG-06-FU1 / RCG-08-FU1 / RCG-09-FU4
deferral, all of which are the same deferral and may be picked up together.

### 3.2 The one CPU failure is a baseline failure, and is proven so

`dipole-open-fft-oracle` fails one of its nine sub-tests:

```text
FAIL: test_energy_finite_difference_derivative
AssertionError: -1.155179279520103e-25 != -1.1551792784777866e-25
  within 2.3103585569555732e-35 delta (1.0423163510038918e-34 difference)
```

The same test was executed at `f8e27f1e` — the last commit before the
remediation program began on 2026-07-30 — and fails there with **bit-identical
numbers**. It is an over-tight finite-difference tolerance (`2e-10` relative
against a ~9e-10 relative difference) in the WP10.1 open-boundary dipole
oracle, which is outside the adaptive coarse-graining package and carries none
of the `cg13-*` labels. It is a baseline failure, distinguished from introduced
failures as parent §12's pull-request review requires, and is not a
coarse-graining release blocker.

### 3.3 fp32 is a separate result, not an inherited one

fp32 does not inherit fp64 acceptance. The executed matrix is in §5. Of the two
previously documented fp32 failures, one is now fixed and one still reproduces;
neither is a newly introduced defect.

### 3.4 The legacy GNU build does not build — diagnosis, not a gate

**Not a supported build path.** The Human owner decided on 2026-08-15 that the
legacy GNU make path is not supported for current GPU builds; CMake is the
supported build system, and its CPU/CUDA matrices above are the release
evidence. What follows is retained as diagnosis for whoever revisits that build
system, not as an unmet requirement.

`make gfortran` fails before compiling anything:

```text
Makefile:184: make/dependencies.make: No such file or directory
make[2]: *** No rule to make target 'make/dependencies.make'.  Stop.
```

`make/dependencies.make` is never written because `make deps`
(`source/make/generateDependencies.py`) aborts:

```text
ERROR: C/CUDA/C++ depedencies not written
       Add <array> to list of excluded libs
```

Diagnosis, established by experiment rather than assumed:

* The generator's `external_CCCC` exclusion list does not know **fifteen**
  headers the current GPU sources include: `<array>`, `<algorithm>`,
  `<limits>`, `<sstream>`, `<stdexcept>`, `<complex>`, `<cassert>`,
  `<iomanip>`, `<memory>`, `<type_traits>`, `<utility>`,
  `<hip/hip_runtime.h>`, `<hiprand/hiprand.h>`,
  `<hip/hip_cooperative_groups.h>`, `<cooperative_groups.h>`.
* With all fifteen added locally, the generator gets further and then fails
  again on a second, independent defect: it resolves quoted includes only
  against `source/gpu_files/`, so a header in a subdirectory
  (`source/gpu_files/measurement/memoryMeasurement.h`) is reported missing.
* The same failure reproduces at `f8e27f1e`, i.e. it **predates the entire
  remediation program**, and does *not* occur on `master`. It is a
  GPU-branch build-tooling defect, unrelated to coarse graining. RCG-01 added
  the ordered `source/CoarseGraining` object set to the legacy Makefiles
  (closing F-02's object-list half); the object list is correct and the
  generator is what is broken.

No fix was applied. The local exclusion-list edits made to diagnose it were
reverted and the clone verified pristine afterwards. Under the 2026-08-15
decision the RCG-10 box "Legacy supported builds pass" is not applicable and
RCG-01's "Supported legacy GNU build succeeds" is withdrawn as a requirement;
RCG-10-FU2 is withdrawn with them. RCG-01's own half of F-02 — the ordered
`source/CoarseGraining` object set in the legacy Makefiles — was completed and
is unaffected.

An initial attempt to build also failed with `make: Permission denied`. That
was an environment artifact — an empty entry in `PATH` makes `make` resolve to
the `source/make/` *directory* — not a repository defect. It is recorded so a
future session does not misdiagnose it.

---

## 4. Evidence added by this task

### 4.1 What the audit is entitled to claim

Remediation §2.1 forbids a claim advancing from reviewer-observed to closed
without a reproducer. Two RCG-09A items were in exactly that state: the CUDA
runtime execution of the parity harness, and three of the five negative
controls. Both were resolved by execution on this host, not by argument.

### 4.2 Device state during the GPU runs

Both RTX A4000s in this host are shared with another user, so device
cleanliness is recorded rather than assumed (the RCG-08-FU3 condition). Before
and during the CUDA runs: utilization 0% on both devices, temperatures 48 °C
and 40 °C, 102 MiB and 15 MiB resident (graphics contexts only, no foreign
compute process). The CUDA results in §3 are pass/fail correctness suites, not
timing measurements, so contention would not change their verdict in any case.

### 4.3 The RCG-09A.4 parity harness executed on a device — first time

`docs/RCG-09A_ASD_PARITY_EVIDENCE.md` records that the direct parity command
returned the project's unavailable-backend code 77
(`ADAPTIVE-ASD-PARITY unavailable: no backend device`) on the validation
container, and that CTest consequently skipped the acceptance test through its
`SKIP_RETURN_CODE`. On this host the test **ran and passed**:

```text
ctest --test-dir build_cuda -R '^coarse-graining-adaptive-asd-parity$' -V
  → Passed, 0.43 s;  205 stages reported result=PASS, 0 result=FAIL
```

Representative stage results, all fixtures (exchange, DMI, uniaxial, combined,
finite-temperature) and both the one- and two-ensemble trajectory cases:

| Stage | Result |
|---|---|
| term-resolved energies | `rel=0.000000e+00`, exact |
| nonzero initial torque | `max_torque=1.561154e+00`, PASS |
| Hamiltonian before predictor | `max_rel≈2.2e-16` (fp64 roundoff) |
| thermal field | **bitwise identical**, tolerance `0.0` |
| Depondt predictor | **bitwise identical** |
| second thermal-field reuse | **bitwise identical** |
| Depondt corrector / final spin | **bitwise identical** |
| six-step finite-T trajectories, M=1 and M=2 | PASS, nonzero displacement at every step |

This is the discriminating evidence parent §12's "do signs and factors match
UppASD Hamiltonian conventions?" needs for the adaptive atomistic kernel, and
it is the reason RCG-09A is a hard dependency of RCG-10.

### 4.4 Three never-executed negative controls, executed

RCG-09A declared five fault injections and executed two (`DMI_SIGN`,
`NO_TRANSPOSE`). The other three compiled but never ran, and were explicitly
**not** reported as passing by inference. All three were built and run here as
separate fault-injection builds. Each fails, returns nonzero, and fails first
at exactly the stage RCG-09A's own table predicts:

| Mutation | Predicted first failure | Observed first failure | Magnitude |
|---|---|---|---|
| `RCG09A_NEGATIVE_THERMAL_AMPLITUDE` | thermal field | `stage=thermal-field` | `max_rel=1.180340e-01` |
| `RCG09A_NEGATIVE_DAMPING` | T=0 Depondt predictor | `stage=depondt-predictor` (T=0 fixture) | `max_rel=1.958363e-08` |
| `RCG09A_NEGATIVE_RNG_DISPLACEMENT` | finite-T thermal field | `stage=thermal-field` (finite-T fixture) | `max_rel=1.462537e+00` |

All five RCG-09A negative controls are now executed evidence. No mutation is
enabled in the production source or in any acceptance target; each is an
`OFF`-by-default CMake option.

---

## 5. fp32 precision matrix

Recorded separately from fp64, per the parent's CG-10 requirement that fp32
error budgets are documented separately and CG-13's rule that fp32 "is not
allowed to inherit fp64 claims".

```sh
cmake -S clone -B build_cuda_fp32 -DCMAKE_BUILD_TYPE=Release \
      -DUPPASD_GPU_BACKEND=CUDA -DUPPASD_PRECISION=SINGLE
ctest --test-dir build_cuda_fp32 -L cg13-cuda --output-on-failure
```

**Result: 31/32 passed, 153.04 s.** The whole coarse-graining operator,
selector, transition, reconstruction, ownership, moving-parity, setup-rejection
and GPU adaptive runtime set passes at fp32, including the RCG-09A.4 ASD parity
target. One test fails:

| Test | Status | Disposition |
|---|---|---|
| `adaptive-cg-production-e2e` | **fails** at `assert abs(float_metric(fft.stdout, "coarse_dipole")) > 0.0` | **RCG-04-FU5**, pre-existing and already tracked. RCG-04J established that this assertion tests a quantity indistinguishable from floating-point noise at *both* precisions; it is a fragile harness assertion, not a demonstrated fp32 precision floor. Not a coarse-graining regression. |

**RCG-05-FU4 is now closed by this run.** The `testSelectorPolicyDescriptorLayout`
float-literal comparison bug — which used to abort the whole
`coarse-graining-gpu-adaptive-runtime` binary under `UPPASD_PRECISION=SINGLE`
and thereby mask RCG-08's fp32 evidence — was fixed in `b235687` (RCG-08) by
comparing `copied.*Threshold` against the pre-copy `policy.*Threshold` value
instead of a hardcoded `double` literal. That test passes here at fp32
(`coarse-graining-gpu-adaptive-runtime`, 0.32 s). The follow-up list is updated
accordingly in §10.

fp32 remains a separate backend/precision result with its own recorded
tolerances, and no fp64 claim is extended to it.

---

## 6. Reopened parent items: current evidence

Remediation §4 provisionally reopened the following. Each row states whether
current evidence restores the acceptance, and what that evidence is.

| Parent area | Reopened acceptance | Verdict | Current evidence |
|---|---|---|---|
| CG-00 | CPU/CUDA/HIP clean builds and legacy build completeness | **Restored for the supported build system** | CPU fp64, CUDA fp64, CUDA fp32 clean out-of-tree CMake builds succeed (§3). HIP unavailable. Legacy GNU make is no longer a supported build (Human decision 2026-08-15), so its failure is not a gate; see §3.4. |
| CG-01 | DMI analytic fixture and convention approval | **Restored** | `coarse-graining-dmi-dimer-energy` hand-derived dimer oracle; Human physics approval recorded 2026-08-08 (Anders Bergman); RCG-02 closed. |
| CG-03 | DMI sign/handedness, small-q atomistic agreement, human sign/prefactor approval | **Restored** | RCG-02 fixed CPU/GPU cross-product order *and* the previously untouched Monte Carlo path; negative control fails with the exact pre-fix sign; `asd-tests` 31/31 with hand-verified golden reconciliation; Human approval 2026-08-08. |
| CG-04 | Coarse DMI chirality and derivative/dispersion evidence | **Restored** | RCG-04-FU4: 24 periodic exchange-only and DMI-only Fourier cases across four block sizes and six wavelengths collapse onto the exact lattice symbols; residuals `8.881e-16` / `9.992e-16`; fitted factors unity; enforced by `coarse-graining-dispersion`. |
| CG-06 | Buffer coverage, CPU/GPU ownership equivalence, non-overlap, moving interface fixtures | **Partially restored** | RCG-05 fixed the directional buffer-width scalarization (`GEO-ANISO-BUFFER`) and proved dipole/short-range/on-site ownership non-overlap; five `adaptive-cg-moving-*` fixtures pass. **Full CPU/GPU ownership-map identity is not met**: RCG-05-FU2 (GPU hard mask not sourced from `cg_static_mask_file`, 44/90 blocks differ on `ownership_aniso_buffer`). |
| CG-07 | Unsupported/unsafe blocks remain atomistic | **Restored** | RCG-03 polarization gate with its own `polarization-unsafe` diagnostic and ratio, on CPU and CUDA; `coarse-graining-polarization-gate`; setup-rejection matrix 31/31. |
| CG-08 | Compensated/low-polarization restriction safety, moving adaptive texture benchmarks | **Restored** | RCG-03 polarization gate; RCG-04 moving wall/skyrmion/DMI-chiral fixtures with recorded torque, displacement, and trajectory distance replacing the vacuous uniform-state oracle. |
| CG-09 | CPU/GPU descriptor equivalence, safe compact-list rebuild | **Partially restored** | Compact-list rebuild is safe: four Compute Sanitizer tools report zero errors (RCG-08), restriction and both Heun stages are bitwise identical to the serial reference, RCG-09C's compaction is deterministic and its host/device oracle disagreement was found and fixed. **Descriptor equivalence is not fully met** — same RCG-05-FU2 gap. |
| CG-10 | CPU/GPU parity, backend execution, performance isolation, crossover, precision budgets, production readiness | **Restored except performance** | fp64 parity, transitions, reconstruction, memory cleanup, and the RCG-09A.4 ASD parity harness all pass on CUDA (§4.3); fp32 budgets recorded separately (§5). **Crossover is not restored** — see §7. |
| CG-10.5 | Clean production dispatch, tracked fixtures, nontrivial trajectories, lifecycle, regression | **Restored** | `adaptive-cg-production-e2e` through the ordinary `sd` binary; `adaptive-cg-fixture-dependencies` audit; RCG-04 nonstationarity observables; `regression-test` and `asd-tests` pass. |
| CG-13 | Automated clean-clone validation, provenance, release wording, performance reporting | **Restored except Human approval** | This document is the clean-clone validation record; `docs/CG-13_RELEASE_VALIDATION.md` is updated with the final matrix and the withdrawn speedup claim. Human release approval is **not** recorded. |
| Section 12 | All review questions | **Completed** — see §8; four questions stay unchecked with named reasons. |
| Definition of success | Items 4–10 require renewed evidence | **Item 9 fails, item 10 fails** — see §7 and §9. |

---

## 7. Performance: what is claimed, and what is withdrawn

### 7.1 The withdrawn claim

Parent CG-10's evidence recorded "a 1.3016x speedup" and "the 1.3017x
active-DOF crossover". Finding F-10 established what that number actually
compared: the adaptive runtime at a low fine fraction against the adaptive
runtime at 100% fine, with a synthetic `featureOffAtomisticWork` FMA loop as
the control. It was never a comparison against UppASD's production atomistic
GPU path. It is withdrawn, and parent CG-10's crossover box is unchecked.

### 7.2 What was actually measured

Measured against `GpuHamiltonianCalculations` and `GpuDepondtIntegrator` — the
objects the feature-off timestep loop calls — and end to end through the
ordinary `sd` binary:

| Task | Production comparison | Outcome |
|---|---|---|
| RCG-09 | 16 384 atoms | adaptive **~3.4 × 10⁴ times slower**; cause isolated to a single-address FP64 CAS energy accumulator, `O(bonds^1.8..2.0)` |
| RCG-09B | five sizes, clean device | serialized accumulator removed at all eight call sites; adaptive still ~7.2–9.4x slower; **no crossover at any size or active fraction** |
| RCG-09C | 16 384 blocks, clean device, ABBA-interleaved | production step **276.6 µs** vs best adaptive step **310.2 µs**; **no production crossover** |

The trend across RCG-09 → 09B → 09C is a four-order-of-magnitude improvement,
and the remaining gap is now ~1.12x rather than ~3.4 × 10⁴. That is a
substantial and honestly-earned engineering result. It is still not a speedup,
and it is not reported as one.

### 7.3 Figures that are current, and figures that are not

* **Current:** RCG-09C's **11.37x** atomistic-phase ceiling at 16 384 blocks /
  65 536 atoms / 114 688 bonds (292.8 µs all-fine → 25.8 µs all-coarse), on a
  device meeting the RCG-08-FU3 cleanliness condition, ABBA-interleaved against
  a pristine RCG-09B build. This is the reduction the coarse-graining mechanism
  achieves *on its own atomistic phase* between its two limits. It is **not** a
  whole-step speedup and **not** a claim against production.
* **Superseded, do not cite:** the `6.0x` ceiling from RCG-09's
  negative-control sweep, and the `1.30x`/`1.3016x`/`1.3017x` active-DOF
  crossover.

### 7.4 Consequence under remediation §13

Items 1–7 and 11–12 pass. **Item 9 fails.** Item 10 is not satisfied (no
independent adversarial review; Human release approval not recorded). §13's
prescribed consequence is explicit and is adopted here: the result is accepted
only as a correctness/reference implementation, the parent blueprint leaves the
GPU production-performance gates open, and speedup wording is removed. Parent
§14 item 9 is amended in place rather than reinterpreted.

---

## 8. Parent section 12 cross-task review

Completed from current evidence. A box is checked only where a citable,
currently-passing artifact supports it. Thirty-three of the thirty-five boxes
are checked; the two that are not each carry a stated reason and are not
reworded to fit. The physics review's finite-temperature question is checked
against the boundary as amended on 2026-08-15, not against the original
wording — the amendment is recorded in parent §14.1 rather than applied
silently.

### 8.1 Physics review

| Question | Verdict | Evidence |
|---|---|---|
| Written energy for every field? | **Yes** | Exchange, DMI, uniaxial anisotropy, interface, and dipole each have a written energy and a derivative/oracle test: `coarse-graining-tensor-operator`, `coarse-graining-dmi-dimer-energy`, `coarse-graining-static-hybrid-operator`, `coarse-graining-adaptive-hamiltonian-contract`. RCG-09A records every `heisge` term as equivalent or setup-rejected, with no term left undetermined. |
| Signs and factors match UppASD conventions? | **Yes** | RCG-02 (one DMI convention across CPU, CUDA, extraction, coarse, hybrid, *and* Monte Carlo), Human-approved 2026-08-08; negative control fails with the exact pre-fix sign. §4.3's bitwise Depondt parity is the end-to-end confirmation. |
| Block-size dependence physically correct? | **Yes** | RCG-04-FU4: collapse versus `q·h` onto the exact exchange and DMI lattice symbols across four block sizes; fitted multiplicative factors unity; no additional coarse scale factor required. |
| DMI spin and spatial tensor indices unambiguous? | **Yes** | RCG-02's dimer oracle plus RCG-09A's `NO_TRANSPOSE` negative control, which fails at the initial Hamiltonian field/energy stage — the unique-pair scatter's transposed contraction is discriminated, not assumed. |
| Compensated systems as channels, not net macrospin? | **Yes, for the accepted scope** | The single-channel production path never forms a net macrospin from a compensated block: RCG-03's non-overridable polarization gate forces such blocks atomistic and logs `polarization-unsafe` with the triggering ratio. Genuine multi-channel dynamics remains CG-11 and is unaccepted; `coarse-graining-multichannel-tensor-operator` covers the operator only. |
| Gyromagnetic ratio and damping reductions justified? | **Yes** | Not reduced: the accepted model requires uniform Landé factor and damping, and setup rejects `Landeg`/`do_site_damping` heterogeneity (`adaptivecgproduction.f90:902`). Fine atoms use the production Depondt path with the production damping array, proven bitwise identical in §4.3. |
| FFT dipole approximation stated correctly? | **Yes** | `docs/CG-13_RELEASE_VALIDATION.md` "Error budget" states grid/padding, reciprocal truncation, basis deposition/interpolation, periodic-image convention, and macro-field approximation, and states that the adaptive mask does not refine the FFT grid. RCG-05F proves dipole ownership is independent of the short-range mask. |
| High-frequency/interface limitations acknowledged? | **Yes** | CG-13 "Error budget and known limitations" — block-scale texture variation, sharp interfaces, high-q spirals, projected ghosts, finite buffer width, block-centre quadrature, discretized ownership; explicitly not a promise of zero error for arbitrary textures. |
| Finite temperature explicitly outside the initial claim? | **Yes, once amended** | Finding RCG-10-F1 (§9.1), resolved by Human capability decision 2026-08-15. Finite temperature is *inside* the claim at the atomistic level (RCG-09A.4 bitwise thermal-field parity) and *outside* it for coarse blocks, which stay deterministic at T=0. The boundary is stated in CG-13 and held by the new `finite_temperature_mixed` fixture with an executed negative control. |

### 8.2 Numerical review

| Question | Verdict | Evidence |
|---|---|---|
| Every operator passes tangent energy derivatives? | **Yes** | Tangent/directional derivative fixtures in the tensor, smooth-projected, static-hybrid, and GPU adaptive runtime tests; RCG-04E/FU4's forward-gradient / discrete-transpose adjoint pair measured exact to `8.881e-16`. |
| Metric supports nonorthogonal cells? | **NO — unchecked** | The operator metric handles skew cells at unit level (RCG-05B generators), but **no skew fixture has ever run through the real executable**: `neighbourmap.f90` rejects the declared exchange bond with "no basis match" during setup. Pre-existing, unrelated to coarse graining. Tracked as RCG-05-FU3. |
| Boundaries and halos correct? | **Yes** | RCG-05D fixed directional buffer-width scalarization (F-07); `GEO-ANISO-BUFFER` is an executed anisotropic-cell negative control that detects the scalarized width. Buffer width is derived from the physical interaction radius. |
| Energy ownership non-overlapping? | **Yes** | RCG-05F ownership invariants; dipole field and energy included exactly once and independent of `onsite_owner`; `adaptive-cg-transition-ownership-invariants`. |
| Transitions synchronized with the integrator? | **Yes** | Complete-step publication with per-block energy-gate rollback; hysteresis and dwell; `adaptive-cg-transition-ownership-invariants` and the adaptive solver fixtures. RCG-09B additionally fixed a real `emom`/`emom2` state-handoff defect found by an existing `8.0e-5` restart budget that was **not** widened. |
| RNG results reproducible? | **Yes** | §4.3: thermal fields bitwise identical, tolerance exactly `0.0`, including the second-draw reuse. RCG-09C: RNG identity independent of active-list position. RCG-06D measured the reconstruction generator's spatial correlation; the zero-cone production default is unaffected by construction. |
| Error thresholds and precision budgets stated? | **Yes** | CG-13 "Precision and diagnostics": `1e-12` relative analytic CPU, `5e-4`/`8e-4` GPU fp64 accumulated-kernel budgets, fp32 recorded separately (§5). |
| Error decreases under spatial refinement? | **Yes** | RCG-04-FU4's `q·h` collapse across four block sizes and six wavelengths; refinement checks in the static-hybrid interface fixtures. |

### 8.3 Software review

| Question | Verdict | Evidence |
|---|---|---|
| Feature-off path unchanged? | **Yes** | RCG-09's paired feature-off delta `+0.023%` with zero inventory change; `regression-test` and `asd-tests` pass on the clean CPU build; feature-off e2e case in the production suite. |
| Immutable topology and mutable state separate? | **Yes** | `GpuAdaptiveDeviceTopology` versus `GpuAdaptiveDeviceRuntime`; the CPU analogue separates `topology` from `runtime` in `adaptive_cg_state`. |
| Spatial, basis, FFT, and dynamical channel counts separate? | **Yes** | CG-02's regular block topology; `cg_channel_mode` gates the dynamical channel count independently of basis and FFT mappings; `do_ralloy`/`Nchmax`/`do_lsf`/`ind_mom_flag` rejected at setup. |
| Unsupported geometries rejected early? | **Yes** | `adaptive-cg-setup-rejection-matrix`, 31/31 cases, each run in an isolated temporary directory requiring a nonzero exit *before* the measurement phase. Explicit-device geometry and non-`P P P` boundaries are rejected. |
| Public interfaces narrow and testable? | **Yes** | Every operator has a standalone test executable; the RCG-09A.1 contract probes are public precisely so the regression test can reach them without widening the runtime interface. |
| CPU and GPU descriptors identical in meaning? | **NO — unchecked** | RCG-05-FU2: the GPU `hardAtomisticBlockMask` is sourced only from the polarization gate, not from `cg_static_mask_file` as CPU does. 44/90 blocks still differ on `ownership_aniso_buffer` after RCG-05D's buffer-width fix. This is a real, independently-scoped defect, not a wording problem. |
| All allocations tracked and released? | **Yes** | RCG-06A memory preflight covers construction and runtime allocations including the scan buffers; the GPU adaptive runtime test asserts exact memory cleanup; `adaptive-cg-mem-large-host` covers the F-13 automatic-array hazard. |
| Compact GPU work lists rebuilt safely? | **Yes** | RCG-08: four Compute Sanitizer tools (memcheck, racecheck, synccheck, initcheck) report zero errors; restriction and both Heun stages bitwise identical to the serial reference; F-12's `O(N log N)` scan replaced with a linear-work stable scan carrying an executed negative control. RCG-09C's ghost-shell flag uses `atomicOr` specifically so a benign write-write race is not merely undetected. |
| Unrelated workflows untouched? | **Yes** | `asd-tests` and `regression-test` pass; the one failing test in the full CPU suite is proven pre-existing (§3.2). |
| μASD separation preserved? | **Yes** | CG-00's constellation removal was the only approved μASD-adjacent change; no remediation task touched `source/Multiscale`. |

### 8.4 Pull-request review

| Question | Verdict | Evidence |
|---|---|---|
| One blueprint task or declared subset per PR? | **Yes** | One `RCG-*` task per session and per commit, RCG-00 through RCG-09D, in the order §10 prescribes. RCG-09A/09B were deliberately split so a possible sign correction and a kernel rewrite did not face the same review (§5.2). |
| Prompt, dependencies, acceptance evidence linked? | **Yes** | Each task in the remediation blueprint carries its prompt, dependencies, checklist, exit evidence, and a linked evidence document. |
| New tests fail without the change and pass with it? | **Yes** | Every physics remediation carries an executed negative control: the DMI dimer oracle fails with the exact pre-fix sign; `GEO-ANISO-BUFFER` detects the scalarized buffer width; the F-12 scan fix carries an executed control; and as of §4.4 **all five** RCG-09A fault injections are executed rather than compiled-only. |
| Baseline failures distinguished from introduced? | **Yes** | §3.2 (dipole oracle, proven bit-identical at the pre-remediation commit), §3.4 (legacy build, proven pre-existing), §5 (the one remaining fp32 failure is pre-existing and tracked). |
| No unrelated formatting hides physics changes? | **Yes** | Reviewed across the RCG series; the sign/convention changes in RCG-02 are isolated to the Hamiltonian routines and their goldens, with the goldens reconciled by hand-verified sign compensation rather than silently regenerated. |
| Sign/unit changes receive explicit physics review? | **Yes** | RCG-02 and RCG-03 both carry recorded Human physics approval (Anders Bergman, 2026-08-08). RCG-09A introduced **no** convention change — it proved the existing convention equivalent to `heisge` — so its "any convention change has independent physics review" box is vacuous rather than unmet. |
| Performance claims include methodology and raw measurements? | **Yes** | RCG-09/09B/09C report medians, MADs, repetitions, raw per-point samples, phase decomposition, launch counts, host waits, device/driver/clock/thermal state, ABBA interleaving, and named the one asymmetry rather than normalizing it. Raw artifacts under `docs/rcg09/` and `docs/rcg08/`. |
| Follow-up limitations recorded, not hidden in TODOs? | **Yes** | Named, independently pickupable follow-up tasks throughout (RCG-04-FU1..5, RCG-05-FU1..4, RCG-06-FU1..2, RCG-08-FU1..4, RCG-09-FU1..5, RCG-10-FU1..3), each with an explicit open/closed/withdrawn state in §10. The coarse-graining sources contain exactly **one** TODO comment (`macrocells.f90:91`, a generalization note), and no limitation is recorded only as a code comment. |

---

## 9. Blocking findings

### 9.1 RCG-10-F1 — the finite-temperature boundary, found ungated and now set

**Resolved by Human capability decision (Anders Bergman, 2026-08-15): finite
temperature is allowed at the atomistic level; coarse blocks stay at T=0.**

#### What the audit found

Finite temperature in the **fine** region was already a deliberate, tested
capability — the runtime prints the boundary at setup (`AdaptiveCG: capability
accepted: production Depondt fine spins; deterministic coarse blocks`), the
tracked e2e case `finite_temperature` runs at `temp 10.0`, and RCG-09A.4 proves
bitwise thermal-field and trajectory parity against production ASD.

But that accepted case is **all-fine by construction** (`initial_coarse=0`), and
nothing exercised or gated the combination of `Temp > 0` with coarse blocks
present. Reproduced by changing one line of the tracked `static_mixed` fixture,
`temp 0.0` to `temp 10.0`:

```text
AdaptiveCG: capability accepted: production Depondt fine spins; deterministic coarse blocks
AdaptiveCG: initial_coarse=1 active_atoms=40 active_blocks=1
AdaptiveCG: resolution_counts fine=1 interface=4 coarse=1
exit status 0
```

Raw record: `docs/rcg10/rcg10_f1_finite_temperature.txt`.

The mechanism: `coarsetensoroperator.f90:167` and
`multichannelcoarsetensoroperator.f90:102` do reject a nonzero
`options%temperature_k` or any stochastic field, but
`adaptivecgproduction.f90:399` sets that field to zero unconditionally, and the
production preflight — 40 `call reject(...)` sites covering explicit-device
geometry, multi-channel, `do_qhb`/`do_3tm`, legacy dipoles, external fields —
never inspects the run's own temperature.

#### The decision, and why the mechanism is now correct rather than bypassed

The owner's decision accepts the combination rather than rejecting it. Under
that decision `adaptivecgproduction.f90:399` is **doing the right thing**: the
coarse operator is *meant* to be handed T=0 while the run itself is warm,
because the continuum blocks carry no thermal fluctuations by design. What was
missing was not a rejection but a statement of the boundary and a fixture
holding it.

#### Evidence debt paid in this session

- **New fixture** `tests/coarse_graining/e2e/finite_temperature_mixed`:
  `static_mixed` at `temp 10.0` with a fixed `tseed 4711`, giving
  `fine=1 interface=4 coarse=1`. `run_production_e2e.py` asserts the capability
  is accepted, that at least one block really is coarse, and that the
  trajectory checksum separates from the T=0 `static_mixed` reference
  (`48.051222` against `48.000000`).
- **Negative control, executed.** Neutralizing the thermal field collapses the
  warm trajectory onto the cold one and the fixture fails with
  `mixed-resolution thermal field did not move the fine atoms: warm=48.0
  cold=48.0`. The check therefore discriminates a thermal field that stops
  reaching fine atoms once a block coarsens — which is the silent failure this
  capability is exposed to.
- **Audit gap closed.** `all_e2e_cases()` in `fixture_dependencies.py` did not
  list `FINITE_TEMPERATURE_CASE`, so the existing finite-temperature fixture's
  inputs were never covered by the `PKG-TRACKED-E2E` audit. Both
  finite-temperature cases are now registered; the audit passes at 61 fixture
  directories and 126 input paths.
- **Documentation.** CG-13 gained a "Temperature" section stating the supported
  and unsupported halves, and an error-budget paragraph: the hybrid is a
  two-temperature-resolution model, magnon populations and thermal broadening
  live only in the atomistic region, and the warm-fine/cold-coarse interface is
  not detailed-balance preserving. No finite-temperature *coarse* physics claim
  is made.

Parent §14 item 10's finite-temperature clause is amended to match. RCG-10-FU1
is closed.

### 9.2 Independent adversarial review — performed

**Performed by the Human owner (Anders Bergman).** That review is what opened
RCG-09A, RCG-09B, RCG-09C and RCG-09D: the adaptive atomistic kernel's
equivalence to `heisge`, the single-address FP64 energy accumulation, the
live-bond compaction, and the quadratic setup algorithms. It is the reason the
performance verdict in §7 exists at all — the `1.30x` claim did not fall to this
reconciliation, it fell to that review and the four tasks it opened.

Every blocker it raised is closed and evidenced. The items it left open are
recorded as unmet rather than argued away: HIP execution evidence (no hardware)
and the performance result itself, which is reported as a negative finding.

This session's own additions to that review — the three previously unexecuted
negative controls (§4.4) and the finite-temperature boundary (§9.1) — are the
kind of thing an independent pass is for, and both are now closed rather than
left as observations.

### 9.3 Human release approval — recorded

Human physics approval exists for the DMI convention and the
polarization/anisotropy thresholds (2026-08-08), for each RCG-0x closure, and
for the capability decisions of 2026-08-15 recorded above.

**Release wording approved as written (Anders Bergman, 2026-08-15)**, including
the §14 item 9 amendment: the speedup claim is withdrawn and the branch is
accepted as a correctness/reference implementation, with the parent's GPU
production-performance gates left open.

The wording this approves, in short: the GPU implementation matches the CPU
reference and reproduces production atomistic dynamics exactly where the model
is all-fine; its performance is characterized, not claimed; no crossover against
UppASD's atomistic GPU path exists at any measured system size or active
fraction. The branch must not be described as a production-performance backend.

### 9.4 What remains unmet

**HIP release evidence.** No HIP toolchain or device exists in any environment
used by this remediation program. Nothing is inferred from CUDA. This is the
single deferral that appears five times under different names (RCG-04-FU1,
RCG-05-FU1, RCG-06-FU1, RCG-08-FU1, RCG-09-FU4) and is one piece of work
awaiting hardware.

## 10. Follow-ups

### Opened and closed by this task

- **RCG-10-F1 / RCG-10-FU1 — finite-temperature boundary. CLOSED**
  (Human capability decision 2026-08-15; fixture, negative control, and
  documentation delivered in this session — §9.1).
- **RCG-10-FU2 — legacy GNU build repair. WITHDRAWN** (Human decision
  2026-08-15: the legacy GNU make path is not a supported build for current GPU
  builds; CMake is the supported build system). The §3.4 diagnosis is retained
  for whoever revisits that build system.

### Opened by this task and still open

- **RCG-10-FU3 — validation-run tree hygiene (non-blocking).** The CPU CG-13
  suite overwrites three tracked `uppasd.adaptive.yaml` run records under
  `examples/AdaptiveCoarseGraining/`, so a validation run dirties the tree.
  Either untrack them or write run records into the build tree.

### Inherited and still open

RCG-04-FU1/2/3/5, RCG-05-FU1/2/3, RCG-06-FU1/2, RCG-08-FU1/4, RCG-09-FU3/4/5.

The single HIP deferral appears five times under different names (RCG-04-FU1,
RCG-05-FU1, RCG-06-FU1, RCG-08-FU1, RCG-09-FU4); it is one piece of work
awaiting hardware, not five.

### Inherited and closed

RCG-04-FU4 (dispersion scaling); RCG-05-FU4 (fp32 descriptor-layout test bug,
fixed in `b235687`, confirmed passing by this task's clean fp32 build);
RCG-09-FU1 (promoted to RCG-09B); RCG-09-FU2 (addressed by RCG-09D).

## 11. Raw artifacts

| File | Contents |
|---|---|
| `docs/rcg10/ctest_cg13_cpu.txt` | CPU fp64 `cg13-cpu`, 29/29 |
| `docs/rcg10/ctest_cpu_full.txt` | full CPU suite, 37/38, with the pre-existing dipole failure |
| `docs/rcg10/ctest_cg13_cuda.txt` | CUDA fp64 `cg13-cuda`, 32/32 |
| `docs/rcg10/ctest_cg13_cuda_fp32.txt` | CUDA fp32 `cg13-cuda` |
| `docs/rcg10/cuda_asd_parity.txt` | RCG-09A.4 staged parity, verbose, 205 PASS / 0 FAIL |
| `docs/rcg10/production_e2e_finite_t_mixed.txt` | production e2e including the new `finite_temperature_mixed` case |
| `docs/rcg10/neg_THERMAL_AMPLITUDE.txt` | fault injection, fails at `thermal-field` |
| `docs/rcg10/neg_DAMPING.txt` | fault injection, fails at `depondt-predictor` |
| `docs/rcg10/neg_RNG_DISPLACEMENT.txt` | fault injection, fails at finite-T `thermal-field` |
| `docs/rcg10/legacy_build.txt` | legacy GNU `make gfortran` failure |
| `docs/rcg10/rcg10_f1_finite_temperature.txt` | RCG-10-F1 reproduction: mixed resolution accepted at `temp 10.0` |
| `docs/rcg10/environment.txt` | commit, `git status`, toolchain, device, and driver state |
