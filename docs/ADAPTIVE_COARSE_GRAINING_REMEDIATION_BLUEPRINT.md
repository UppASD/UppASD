# Adaptive two-scale coarse graining: review remediation

## Physics, validation, build, numerical-efficiency, and release-recovery blueprint

**Status:** Working companion blueprint for sequential remediation sessions  
**Date:** 2026-07-30  
**Parent specification:** `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md`  
**Primary scope:** Close the independent-review findings against the
single-channel adaptive coarse-graining implementation, restore trustworthy
acceptance evidence, and return to the parent blueprint for final review  
**Explicitly out of scope:** New finite-temperature, explicit-device, general
multi-channel, or μASD functionality

---

## 1. Executive decision

The independent review has material bearing on the release state of the
adaptive coarse-graining implementation. It found both ordinary engineering
defects and silent-wrong-physics risks. Several central findings are directly
confirmed by the current source, and the clean CPU build failure has been
reproduced independently.

The implementation architecture is not rejected. The adjoint projection,
exchange prefactor, mixed-derivative symmetrization, anisotropy derivative,
coarse DMI derivative, and interaction-radius-derived buffer are valuable
assets that must be preserved. The remediation therefore proceeds as a
gated correction program rather than a redesign.

This document is a temporary execution overlay on the parent blueprint:

- the parent blueprint remains the authoritative physical and architectural
  specification;
- this document records review findings, reproduction status, reopened gates,
  corrective tasks, and acceptance fixtures;
- checked parent items contradicted by the review are considered
  **provisionally reopened**, even before their Markdown boxes are edited;
- no release, parity, chirality, or speedup claim is restored merely because a
  code patch compiles;
- after every remediation gate is accepted, the parent blueprint is updated
  with new evidence and becomes the final checklist again.

The work should be performed in sequential, reviewable sessions. Each session
should normally implement one `RCG-*` task, attach its evidence, and stop at
its exit gate. Physics-convention changes, GPU performance work, and final
release reconciliation require independent review.

---

## 2. Evidence policy

### 2.1 Finding states

Every review claim is assigned one of four states:

1. **Confirmed:** reproduced by a clean command or directly demonstrated by
   source with no material interpretive ambiguity.
2. **Corroborated:** the source structure supports the finding, but its
   numerical or performance consequence still requires a controlled fixture.
3. **Reviewer-observed:** reported from an environment or device not yet
   independently reproduced in this worktree.
4. **Closed:** corrected and accepted with the evidence required by this
   blueprint.

No claim advances directly from reviewer-observed to closed. It must first
receive a reproducer or a documented reason why reproduction is impossible.

### 2.2 Clean evidence

Build and test acceptance must come from tracked source in a fresh out-of-tree
build directory. Reusing an existing build directory is diagnostic only.

This requirement is already justified by the branch:

- an incremental `build_cpu` build completed because it retained
  `adaptivetimestepping.mod` and its object;
- a fresh out-of-tree CPU build failed at
  `source/MonteCarlo/montecarlo.f90:58` because the source still imports the
  deleted `AdaptiveTimeStepping` module.

Every acceptance record must therefore include:

- the exact commit;
- `git status --short`;
- configure command and relevant options;
- build command and result;
- test command and complete result summary;
- compiler/backend/precision;
- device and driver for hardware-dependent evidence;
- raw benchmark output for performance claims.

An acceptance run is invalid if it depends on untracked fixtures, stale module
files, manually copied binaries, or a prior build tree.

### 2.3 Negative controls

New tests must demonstrate that they distinguish the defect from the fix.
For each physics or validation remediation, record one of:

- the new test fails on the pre-fix commit and passes on the fix;
- a minimal fault injection fails while the accepted implementation passes;
- an analytic oracle detects the old result with a nonzero margin.

Merely executing a code path, producing a nonzero energy, or matching a
uniform fixed point is not sufficient evidence.

### 2.4 Claim discipline

Until the final gate:

- the GPU backend is described as a correctness prototype, not a production
  performance backend;
- DMI production support is unaccepted;
- current end-to-end parity is unaccepted;
- the reported `1.30x` active-DOF crossover is not a comparison with UppASD's
  real atomistic GPU path;
- the parent section 12 review boxes remain unchecked.

---

## 3. Review finding register

### 3.1 Confirmed findings

| ID | Finding | Evidence in the current branch | Consequence |
| --- | --- | --- | --- |
| F-01 | A clean CPU build fails | `montecarlo.f90` imports deleted `AdaptiveTimeStepping`; fresh build fails on the missing module | All later binary and release evidence must be regenerated |
| F-02 | Legacy Makefiles omit the new package | `source/Makefile.gfortran`, `source/Makefile.legacy`, and `source/Makefile.legacy.gfortran` contain neither the removed-module repair nor the new coarse-graining object set | Supported non-CMake builds cannot represent the production dependency graph |
| F-03 | Five referenced e2e directories are not tracked | The harness references `missing_alat`, two anisotropy cases, and two DMI/anisotropy cases; their local directories are absent from `git ls-files` | A clean clone cannot execute both advertised e2e suites |
| F-04 | CPU and GPU atomistic DMI cross-product orders differ | `HamiltonianActions:dzyaloshinskii_moriya_field` evaluates the component order of `m_j × D_ij`; `hamdev::dm_field` evaluates `D_ij × m_j` | CPU/GPU chirality and CPU hybrid interface physics cannot both match without a convention decision |
| F-05 | The small-\(q\) DMI validator follows only one of those signs | `validate_coarse_material_small_q` adds the documented positive sine term while the CPU atomistic field uses the opposite cross-product order | Internal coarse self-consistency does not establish CPU atomistic consistency |
| F-06 | Production anisotropy is taken from one central cell | `build_production_anisotropy` constructs atom data, then copies anisotropy from `countstart+basis` to every block | Spatially varying anisotropy can be silently homogenized |
| F-07 | GPU buffer dilation loses directional widths | Production staging sends `maxval(buffer_width_blocks)` as one scalar | Anisotropic cells produce different conservative CPU/GPU ownership masks |
| F-08 | CPU buffer dilation is quadratic in block count | `rebuild_static_hybrid_ownership` scans every block for every fine seed | Adaptive rebuild cost can become prohibitive |
| F-09 | Six adaptive GPU kernels are explicitly single-threaded | Restriction, publication, atomistic evaluation, energy finalization, and both Heun stages reject every thread except `(0,0)` | The hot adaptive path is serialized |
| F-10 | The feature-off benchmark is synthetic | `featureOffAtomisticWork` is an FMA loop, not UppASD's normal GPU Hamiltonian and integrator | It tests inactive overhead, not production speedup |
| F-11 | GPU energy terms use `real` storage | `energyTerms_` is `GpuTensor<real,1>` | FP32 global energy accumulation lacks an FP64 accumulator |
| F-12 | The GPU buffer scan is Hillis--Steele style | `scanAdaptiveWorkStep` is launched over doubling offsets | Compaction performs `O(N log N)` element work and multiple launches |
| F-13 | Large automatic arrays depend on `n_atoms` | Static hybrid and smooth projected routines declare multiple automatic full-system arrays | Large runs can overflow the host stack |
| F-14 | No polarization acceptance gate is present | Resultant and moment sums exist, but no coarsening decision compares their ratio with a threshold | Nearly compensated single-channel blocks can acquire a noise-defined direction |

### 3.2 Corroborated findings requiring controlled measurement

| ID | Finding | Required reproducer before closure |
| --- | --- | --- |
| F-15 | Existing e2e parity is physically weak or vacuous | Record initial torque, maximum per-spin displacement, trajectory distance, and a deliberately non-equilibrium reference |
| F-16 | The CPU hybrid obtains little or no coarse-fraction scaling | Measure field and step wall time versus coarse fraction after separating setup, ownership rebuild, Hamiltonian, and integration |
| F-17 | Per-step full-system allocation churn is material | Count allocations and measure their wall-time/memory effect in a representative production case |
| F-18 | CPU timers exclude dominant work and use unsuitable clock semantics | Compare reported phase totals with external wall time and include field evaluation |
| F-19 | The pending-state dilation write/read pattern is a fragile race invariant | Run race/sanitizer checks and prove that only monotonic, stable states are read during the kernel |
| F-20 | The tuple-seeded MINSTD reconstruction can correlate nearby refinements | Add a spatial/statistical fixture before deciding whether replacement is required for the supported zero-cone default |
| F-21 | Dipole ownership is correct but insufficiently explicit | Add an ownership assertion and documentation proving why dipole ignores short-range `onsite_owner` |

### 3.3 Reviewer-observed hardware findings

The following measurements are useful review evidence but must be rerun after
correctness and implementation remediation:

- adaptive all-atomistic GPU execution was reported roughly three orders of
  magnitude slower than the synthetic feature-off control;
- the reported `1.30x` result compared two modes of the adaptive runtime, not
  adaptive CG with UppASD's production atomistic GPU path;
- no HIP device execution evidence was available;
- FP32 energy tolerances may be large enough to mask accumulation error.

Their precise numeric values are not frozen as acceptance references. New raw
measurements supersede them.

### 3.4 Positive findings to preserve

The review also establishes regression contracts:

- `physical_forward_gradient` and
  `add_physical_gradient_transpose` form an exact discrete adjoint pair;
- the `0.25` exchange-stiffness prefactor is consistent with UppASD's energy
  and field ownership convention;
- mixed derivatives are second order when the full symmetric exchange tensor
  contraction is used;
- the implemented anisotropy energy and analytic derivative are
  self-consistent;
- the coarse DMI energy derivative contains both variation terms with the
  correct relative weight;
- normalized prolongation uses the tangent projection and normalization
  Jacobian required by the chain rule;
- CPU buffer width is derived from the physical interaction radius.

Each affected refactor must retain the corresponding unit or derivative test.
These findings are not permission to keep the conflicting DMI convention or
other release defects.

---

## 4. Provisional reopening of the parent blueprint

The following parent tasks are reopened for evidence purposes. Existing
checked boxes may remain as historical records while remediation is active,
but they must not be cited as current acceptance.

| Parent area | Reopened acceptance |
| --- | --- |
| CG-00 | CPU/CUDA/HIP clean builds and legacy build completeness |
| CG-01 | DMI analytic fixture and convention approval |
| CG-03 | DMI sign/handedness, small-\(q\) atomistic agreement, and human sign/prefactor approval |
| CG-04 | Coarse DMI chirality and relevant derivative/dispersion evidence |
| CG-06 | Buffer coverage, CPU/GPU ownership equivalence, non-overlap, and moving interface fixtures |
| CG-07 | Unsupported/unsafe blocks remain atomistic |
| CG-08 | Compensated/low-polarization restriction safety and moving adaptive texture benchmarks |
| CG-09 | CPU/GPU descriptor equivalence and safe compact-list rebuild |
| CG-10 | CPU/GPU parity, backend execution, performance isolation, crossover, precision budgets, and production readiness |
| CG-10.5 | Clean production dispatch, tracked executable fixtures, nontrivial trajectories, lifecycle evidence, and regression results |
| CG-13 | Automated clean-clone validation, provenance, release wording, and performance reporting |
| Section 12 | All physics, numerical, software, and pull-request review questions remain open |
| Definition of success | Items 4 through 10 require renewed evidence |

If a remediation task changes the underlying contract rather than restoring
it, the parent text must be amended explicitly and receive human approval.

---

## 5. Remediation principles

### 5.1 Correctness before optimization

The order is mandatory:

1. make clean builds and tests honest;
2. resolve conventions and silent-wrong-answer behavior;
3. add discriminating executable fixtures;
4. repair memory and geometry parity;
5. optimize CPU and GPU paths;
6. measure against production baselines;
7. reconcile release claims.

Performance work performed before the DMI, polarization, and fixture gates may
be kept as exploratory work but cannot receive acceptance.

### 5.2 One defect class per patch

Build cleanup, DMI convention changes, capability changes, test-fixture work,
host performance, device performance, and release reconciliation should be
separate commits or pull requests. This preserves bisectability and allows
physics sign changes to receive focused review.

### 5.3 Preserve the feature-off path

Repairs must not route feature-disabled runs through adaptive setup, storage,
Hamiltonian, synchronization, or integration. The clean feature-off baseline
is captured before cross-cutting performance work and compared again at every
later gate.

### 5.4 Reject rather than homogenize

Where the implemented model cannot represent an input, setup must reject it
with a named diagnostic. It must not select a representative cell, silently
average incompatible descriptors, or continue with a noise-defined
macrospin.

### 5.5 Same semantics on CPU, CUDA, and HIP

Descriptors may have backend-specific storage, but not backend-specific
meaning. Directional buffer widths, DMI handedness, polarization state,
energy precision, work ownership, and diagnostics must be equivalent.

---

## 6. Fixture architecture

### 6.1 Fixture layers

The remediation suite uses four layers:

1. **Analytic unit fixtures:** dimers, chains, tensor derivatives, threshold
   boundaries, and ownership masks.
2. **Operator fixtures:** all-coarse, smooth projection, static hybrid, and
   adaptive solver calls with exact or high-accuracy references.
3. **Production executable fixtures:** ordinary `inpsd.dat` inputs running the
   normal UppASD binary without direct internal staging.
4. **Performance fixtures:** production CPU/GPU baselines and adaptive sweeps
   with hardware metadata and raw phase timings.

Passing a lower layer does not substitute for the next layer.

### 6.2 Required observables

Nontrivial production fixtures should record:

- initial and final per-spin state or a stable hash/checksum sensitive to
  permutation and sign;
- maximum and RMS spin displacement from the initial state;
- initial torque norm and at least one later torque norm;
- total and named per-term energies;
- atomistic, coarse, interface, and dipole field checksums;
- fine, buffer, and coarse block counts;
- accepted/rejected transitions;
- phase wall times where performance is in scope.

`direction_sum` and `direction_norm2` may remain diagnostics, but they are not
sufficient parity or dynamics oracles.

### 6.3 Mandatory fixture catalogue

| Fixture ID | Layer | Purpose | Required oracle |
| --- | --- | --- | --- |
| BLD-CLEAN-CPU | Build | Detect stale modules and missing dependencies | Fresh CPU build succeeds from tracked source |
| BLD-CLEAN-CUDA | Build | Validate CUDA dependency graph | Fresh CUDA build succeeds on an available toolchain |
| BLD-CLEAN-HIP | Build | Validate HIP dependency graph | Fresh HIP build succeeds on an available toolchain |
| BLD-LEGACY-GFORTRAN | Build | Preserve supported legacy build | Legacy binary compiles and links all production CG modules |
| PKG-TRACKED-E2E | Packaging | Prevent local-only fixtures | Harness dependency inventory is a subset of `git ls-files` |
| DMI-DIMER-ENERGY | Analytic | Freeze Hamiltonian sign | Hand-computed dimer energy and both spin fields match |
| DMI-DIMER-CPU-GPU | Analytic/backend | Freeze backend parity | CPU, CUDA, and HIP fields match the accepted convention |
| DMI-SPIRAL-Q | Operator | Validate spiralization sign | Fitted minimum has accepted sign and converges at small \(q\) |
| DMI-HYBRID-CROSSING | Production | Validate interface handedness | A chiral texture crosses the interface without chirality reversal or spurious pinning beyond budget |
| POL-THRESHOLD | Analytic | Freeze polarization inequality | Values below, at, and above the threshold select documented states |
| POL-CANCELLATION | Operator | Prevent noise macrospins | Low-resultant single-channel blocks remain fine; exact zero never normalizes |
| ANI-UNIFORM-TRANSLATED | Operator/e2e | Remove central-cell dependence | Translating the unit-cell origin leaves uniform results unchanged |
| ANI-NONUNIFORM-REJECT | Production | Enforce capability boundary | Unsupported spatial variation fails before integration with a named diagnostic |
| E2E-MOVING-OFF-FINE | Production | Validate feature-off/all-fine identity | Nonzero torque and displacement; complete trajectories match within fp64 budget |
| E2E-MOVING-ALL-COARSE | Production | Validate long-wave dynamics | Frequency, phase, energy, and trajectory error meet the discretization budget |
| E2E-MOVING-STATIC | Production | Validate mixed ownership | Nonzero atomistic/coarse/interface work and decreasing interface error |
| E2E-MOVING-ADAPTIVE | Production | Validate transitions | Texture moves, transitions occur, and energy jumps remain bounded |
| GEO-ANISO-BUFFER | Operator/backend | Preserve directional halo semantics | CPU/CUDA/HIP block-state maps match exactly on unequal widths |
| MEM-LARGE-HOST | Robustness | Detect stack/per-step allocation failure | Large representative case runs under a documented stack limit |
| ENERGY-FP32-ACCUM | Backend | Validate accumulator precision | Error scaling remains below a separately derived FP32 budget |
| PERF-ATOMISTIC-PROD | Performance | Establish real baseline | Normal UppASD Hamiltonian plus integrator is timed |
| PERF-CG-SWEEP | Performance | Establish useful crossover | Same physics/input is swept over coarse fractions and compared with PERF-ATOMISTIC-PROD |

### 6.4 Moving-state construction

At least one production parity fixture must begin away from a stationary
state. Preferred deterministic choices are:

- a spin spiral whose wavevector is deliberately offset from the DMI energy
  minimum;
- a translated or perturbed domain wall whose width differs from
  \(\sqrt{A/K}\);
- a chiral chain with a known nonzero initial torque.

The fixture must prove nontriviality by asserting:

- initial torque norm is above a documented floor;
- maximum spin displacement after integration is above a documented floor;
- the oracle changes if the DMI sign, coarse operator, or integration update
  is disabled.

Increasing `Nstep` without these assertions is not sufficient.

### 6.5 Fixture provenance

Every committed fixture directory must include either a short README or
machine-readable metadata recording:

- physical model and supported capability;
- source of the analytic/reference value;
- units and sign convention;
- why the state is nonstationary or otherwise discriminating;
- precision/backend tolerances;
- expected failure mode before the fix.

Generated restarts must be reproducible from a committed generator or
documented ordinary UppASD input. Opaque locally generated restarts are not
accepted as sole references.

---

## 7. Phased remediation and quality gates

### Gate A: Evidence integrity

Deliver clean builds, legacy build repair, tracked fixtures, and CI jobs that
cannot consume stale module files.

**Exit gate:** A fresh CPU build and the CPU coarse-graining suite run from
tracked source; CUDA/HIP build availability is reported honestly; all harness
fixture paths are tracked.

### Gate B: Physics safety

Deliver the accepted DMI convention, cross-backend dimer/spiral tests,
polarization gate, anisotropy capability handling, and tensor symmetry
assertion.

**Exit gate:** CPU atomistic, GPU atomistic, material extraction, all-coarse,
and hybrid DMI share one signed analytic oracle; unsafe coarse blocks and
unsupported anisotropy reject or stay atomistic as specified.

### Gate C: Discriminating validation

Deliver moving feature-off, all-fine, all-coarse, static, adaptive, DMI, and
anisotropic-geometry executable fixtures.

**Exit gate:** End-to-end tests detect disabled/broken coarse dynamics, show
nonzero evolution, and compare complete state-sensitive observables.

### Gate D: Runtime robustness

Deliver persistent host workspaces, controlled allocation lifetime, vector
buffer semantics, FP64 energy accumulation, safe device-state publication,
and corrected timings.

**Exit gate:** Large-system and precision fixtures pass; CPU/GPU masks match;
reported phases account for the adaptive step within a documented residual.

### Gate E: Production efficiency

Deliver local CPU dilation, threaded CPU work, parallel device atomistic and
integration kernels, and linear-work compaction where supported.

**Exit gate:** Scaling evidence shows the expected algorithmic behavior and
correctness fixtures remain unchanged.

### Gate F: Release recovery

Deliver production-baseline benchmarks, backend matrices, independent review,
and parent-blueprint reconciliation.

**Exit gate:** The original definition of success is either supported by new
evidence or explicitly narrowed by an approved scope amendment.

---

## 8. Delegation guide

The labels follow the parent blueprint and describe capabilities rather than
provider requirements:

- **Human:** owns Hamiltonian conventions, capability/scope changes, accepted
  model error, and final release wording.
- **Opus/Terra:** high-reasoning derivation, cross-module architecture,
  adversarial physics review, and performance-methodology review.
- **Sol:** substantial Fortran/C++ implementation, CPU/GPU integration,
  debugging, kernel engineering, and profiling.
- **Luna/Sonnet:** bounded fixes, clean-build automation, fixture construction,
  focused tests, documentation, and independent evidence audits.

Required separation:

- the implementer of the DMI sign change must not be its sole physics reviewer;
- the author of GPU kernel parallelization must not be the sole parity or
  performance reviewer;
- the author of an e2e oracle must demonstrate its negative control to a
  separate reviewer;
- final parent-checklist reconciliation requires Human approval.

Sequencing rules:

- accepted remediation proceeds sequentially from RCG-00 through RCG-10
  because later work depends on accepted conventions, fixtures, descriptors,
  and reference results;
- exploratory analysis may happen separately, but no later task is accepted
  before all of its declared dependencies;
- CPU and GPU buffer semantics must not be changed independently;
- no two agents should concurrently edit central production dispatch,
  `adaptivecgproduction.f90`, or the DMI convention sources;
- performance tuning must rebase on the accepted correctness commit rather
  than carrying an older private copy of the kernels.

Every session starts by reading both blueprints and the evidence from its
dependencies. It ends by updating only its own checklist and evidence block.

---

## 9. Task prompts and acceptance checklists

### Task RCG-00: Freeze claims and capture a clean baseline

**Dependencies:** None  
**Suggested primary:** Luna/Sonnet  
**Suggested review:** Sol for build graph; Human for claim freeze  
**Risk:** Low code risk, high evidence importance

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Establish a reproducible pre-remediation baseline from tracked source.
> Inventory every adaptive-CG build target, CTest label, executable fixture,
> backend, precision, and claimed benchmark. Configure a new out-of-tree CPU
> build with no reusable module or object files and capture the complete
> failure. Inspect incremental build directories for stale
> `AdaptiveTimeStepping` artifacts and record why they can conceal the defect.
>
> Add a machine-checkable fixture-dependency audit for
> `run_production_e2e.py` and `run_setup_rejection_matrix.py`: every referenced
> fixture directory and required file must be tracked. Do not fix physics or
> optimize kernels in this task. Mark DMI, e2e parity, and GPU speedup claims
> suspended in the remediation evidence.

#### Checklist

- [x] The baseline commit and dirty/untracked state are recorded.
- [x] A fresh CPU build reproduces or disproves F-01.
- [x] Stale module/object artifacts that concealed the result are identified.
- [x] CMake, legacy, CUDA, and HIP build paths are inventoried.
- [x] Every e2e harness path is checked against tracked files.
- [x] Current CTest labels and backend membership are recorded.
- [x] Existing benchmark claims are classified by their actual baseline.
- [x] No implementation or physics behavior changes are mixed into the patch.
- [x] Human acknowledges the temporary release-claim freeze (2026-07-30).

**Exit evidence:** `BLD-CLEAN-CPU` pre-fix log, tracked-fixture audit output,
and a baseline matrix attached to the task.

**RCG-00 evidence (2026-07-30):** The baseline matrix, clean-build failure,
stale-artifact inspection, fixture-audit output, and suspended-claim register
are recorded in `docs/ADAPTIVE_COARSE_GRAINING_BASELINE_20260730.md`.

---

### Task RCG-01: Restore build and packaging completeness

**Dependencies:** RCG-00  
**Suggested primary:** Luna/Sonnet or Sol  
**Suggested review:** Luna/Sonnet not responsible for the patch  
**Risk:** Medium mechanical breadth

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Remove the dangling `AdaptiveTimeStepping` import and repair every supported
> build description. Add the complete ordered `source/CoarseGraining` object
> set to the legacy Makefiles, including every module required by
> `uppasd.f90` and `sd_driver.f90`. Commit all e2e fixture files intentionally
> required by the harness, or remove a harness case only if the corresponding
> capability claim is explicitly withdrawn.
>
> Add CI coverage that configures and builds from a clean checkout. The job
> must not reuse Fortran module directories between commits. Preserve the
> feature-off path and avoid unrelated formatting or physics edits.

#### Checklist

- [x] `AdaptiveTimeStepping` has no remaining live import or build artifact.
- [x] Fresh CMake CPU configuration and build succeed.
- [ ] Supported legacy GNU build succeeds and links the production binary.
- [ ] CUDA and HIP clean builds succeed where toolchains are available.
- [ ] Backend absence is reported as unavailable, not silently passed.
- [x] Every e2e fixture referenced by either harness is tracked.
- [ ] Fixture inputs, interaction files, masks, and restarts are complete.
- [x] CI includes a genuinely clean build job.
- [x] Feature-off regression tests pass.
- [ ] μASD and unrelated workflow regressions are unchanged.
- [x] The pre-fix clean-build negative control is attached.

**Exit evidence:** `BLD-CLEAN-*`, `BLD-LEGACY-GFORTRAN`, and
`PKG-TRACKED-E2E`.

**RCG-01 evidence (2026-07-30):** The tracked source contains no live
`AdaptiveTimeStepping` import or build entry. A fresh out-of-tree CMake CPU
configuration and Release build succeeded, as did the fixture dependency audit
(`PKG-TRACKED-E2E`), the production e2e suite, and the feature-off case. The
new CI job configures and builds in a commit-specific runner-temporary
directory, then runs the tracked-fixture audit and production e2e test.

The legacy GNU build is not marked complete: its dependency-generation step
stops on a pre-existing C/C++ dependency-generator error before compilation.
CUDA/HIP clean-build and backend-unavailability evidence were not produced in
this local environment, and the broader setup-rejection label retains
pre-existing failures outside RCG-01. Generated restart/output files remain
runtime products; the audit covers every input, interaction file, and mask
consumed by the harness from tracked source.

---

### Task RCG-02: Freeze and repair DMI handedness end to end

**Dependencies:** RCG-01  
**Suggested primary:** Opus/Terra for derivation; Sol for implementation  
**Required independent review:** Human physics owner and an independent
Opus/Terra or Sol reviewer  
**Risk:** Critical silent physics change

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Starting from the written UppASD Hamiltonian, derive one indexed DMI energy
> and both site fields without relying on existing code as the oracle. State
> the directed-neighbour convention, pair-counting factor, antisymmetry
> assumptions, cross-product order, moment/field units, and expected
> right-/left-handed dimer response.
>
> Implement `DMI-DIMER-ENERGY` first. Run it against the CPU atomistic action,
> `applyhamiltonian`, the GPU `dm_field`, stiffness/spiralization extraction,
> the small-\(q\) validator, the coarse tensor operator, and hybrid production.
> Resolve every disagreement at the source of the convention; do not flip
> only the coarse sign to make one internal validator pass.
>
> Add a small-\(q\) chiral-chain or film fixture with the analytic
> \(q_{\min}\) sign and a mixed-resolution interface fixture. Record the
> pre-fix failures. Any change to normal atomistic DMI requires focused
> feature-off CPU/GPU regression evidence and Human approval.

#### Checklist

- [x] One written energy and field derivation fixes all indices and signs.
- [x] Directed versus unique-pair counting is explicit.
- [x] `DMI-DIMER-ENERGY` matches the hand calculation.
- [x] CPU atomistic field matches the accepted dimer field.
- [x] `applyhamiltonian` and all active CPU Hamiltonian paths match.
- [x] CUDA and HIP device fields match where available. (HIP unavailable in
      every environment used so far; CUDA passes.)
- [x] Material spiralization has the accepted sign.
- [x] The direct atomistic small-\(q\) energy uses the production convention.
- [x] The coarse tensor energy derivative still passes.
- [x] The analytic \(q_{\min}\) sign and magnitude converge at small \(q\).
- [x] CPU hybrid atomistic and coarse regions prefer the same chirality.
- [x] A sign-reversed negative control fails.
- [x] Feature-off DMI regressions are attached.
- [x] Human physics approval is recorded.

**Exit evidence:** `DMI-DIMER-ENERGY`, `DMI-DIMER-CPU-GPU`,
`DMI-SPIRAL-Q`, and the operator-level portion of
`DMI-HYBRID-CROSSING`.

**RCG-02 evidence (2026-08-08, CLOSED):**
`docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md` records the independent indexed
derivation, the DMI dimer's pre-fix CPU-action failure, and Human approval of
the source-level correction (2026-07-31).  `DMI-DIMER-ENERGY`, material
small-\(q\), tensor DMI derivative/chiral-chain, static-hybrid crossing, and
the ordinary feature-off/production CPU suite pass on the corrected source.
A fresh CUDA fp64 build on RTX A4000 hardware passes the device dimer, GPU
adaptive runtime, and production-e2e fixtures; HIP is unavailable.

Two further gaps were found and closed on 2026-08-08, independent of the
DMI fix's original author: (1) the legacy feature-off ASD suite's four
DMI-sensitive golden-output changes (Kagome, SCsurf) are reconciled by an
exact, hand-verified sign compensation in `tests/kagome/dmfile` and
`tests/SCsurf/dmdata` rather than by changing the golden references, and
`ctest -R '^asd-tests$'` now passes 31/31 including those four cases; (2)
Monte Carlo mode (`calculate_efield`, `calculate_energy` in
`source/MonteCarlo/montecarlo.f90`/`montecarlo_common.f90`, plus a related
`emom`/`emomM` mixing defect also present in `spinice.f90`) still used the
pre-fix DMI handedness — untouched by the original commit and uncovered by
this session's audit of every active CPU Hamiltonian path. No committed
regression exercised Monte Carlo with DMI enabled, so this carried no
silent-wrong-published-result risk, but it did fail the "one convention
governs every active path" gate. Both are now fixed, with a hand-derived
negative control (`tests/coarse_graining/test_dmi_dimer_energy.f90`) that
fails with the exact pre-fix sign when only the fix is reverted and passes
otherwise; `ctest -L cg13-cpu` (12/12) and `ctest -R '^asd-tests$'` (31/31)
both pass unchanged on a fresh out-of-tree CPU build after the fix.

Human physics approval of both 2026-08-08 findings (golden reconciliation
and the Monte Carlo fix) is recorded (Anders Bergman, 2026-08-08). Both are
committed (`5275f733`, `2f700fda`); a clean-commit rebuild of `2f700fda`
(`git describe` reports no `-dirty` suffix) reproduces every result above on
both CPU and CUDA.

**RCG-02 is closed (2026-08-08, Human decision: Anders Bergman).** Two
items are explicitly deferred by approved scope decision, not blocking:
HIP execution evidence, because no HIP toolchain or device exists in any
environment used so far (postponed to whenever HIP hardware becomes
available); and a separate independent Opus/Terra or Sol adversarial
physics review distinct from Human approval, deferred to a later stage of
the remediation program. Neither deferral reflects a physics disagreement
or an unresolved correctness question — every fixture this document
requires passes now.

---

### Task RCG-03: Enforce coarse-model capability safety

**Dependencies:** RCG-02  
**Suggested primary:** Opus/Terra for threshold/capability design; Sol for
production wiring  
**Suggested review:** Luna/Sonnet for rejection matrix; Human for threshold
and supported anisotropy scope  
**Risk:** High silent-wrong-answer prevention

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Implement the section 4.10 polarization gate using the already available
> resultant and scalar moment sums. Define
> \(p_b=\lVert\sum_a\mathbf M_a\rVert/\sum_a\lVert\mathbf M_a\rVert\),
> its numerical floor, the accepted threshold, equality behavior, selector
> reason code, diagnostics, and interaction with hysteresis/dwell. A
> single-channel block below the threshold must remain or become atomistic.
> Multi-channel cancellation is not enabled by this task.
>
> Audit production anisotropy construction. For the initial supported model,
> either construct block-local descriptors correctly or prove exact
> cell-periodicity and reject any spatial variation before integration. Do
> not use one central cell as an implicit representative. Add an exchange
> tensor symmetry assertion wherever second-order mixed-derivative accuracy
> relies on symmetry.

#### Checklist

- [x] Polarization is defined with units, floors, and equality semantics.
- [x] Below-threshold single-channel blocks cannot be coarsened.
- [x] Already-coarse unsafe blocks refine at an accepted synchronization point.
- [x] Hysteresis/dwell cannot override the safety interlock.
- [x] Zero and near-zero resultants never obtain noise-defined directions.
- [x] Selector diagnostics report polarization and the reason for refinement.
- [x] CPU and GPU use the same threshold and comparison.
- [x] Threshold tests cover below, equal, above, and roundoff-scale values.
- [x] Uniform anisotropy remains supported and translation invariant.
- [x] Unsupported spatially varying anisotropy rejects before integration.
- [x] No central-cell sampling remains as an unstated model.
- [x] Required exchange tensor symmetry is checked during setup.
- [x] Multi-channel behavior is unchanged and remains separately gated.
- [x] Human approves the polarization threshold and anisotropy boundary.

**Exit evidence:** `POL-THRESHOLD`, `POL-CANCELLATION`,
`ANI-UNIFORM-TRANSLATED`, and `ANI-NONUNIFORM-REJECT`.

**RCG-03 evidence (2026-08-06, exploratory, not accepted):**
`docs/RCG-03_POLARIZATION_ANISOTROPY_EVIDENCE.md` records a dynamic,
non-overridable polarization gate (new `cg_polarization_threshold` input,
default `0.9`) reusing the existing `hard_atomistic_mask` interlock, and a
cell-periodicity uniformity check in `build_production_anisotropy` plus a
matching exchange-tensor-symmetry assertion in
`MultiChannelCoarseTensorOperator`. `POL-THRESHOLD` and `POL-CANCELLATION`
pass with a defect-injection negative control at both the unit and
production-fixture layers, on CPU and now also on CUDA (the polarization
gate was ported to `GpuAdaptiveRuntime` on real RTX A4000/CUDA 13.3
hardware, since the GPU backend runs its own independent selector kernels
and never consulted the CPU mask; device-level and production-e2e negative
controls both pass; HIP remains untested by explicit instruction).
`ANI-UNIFORM-TRANSLATED` passes as a regression (this check needed no GPU
port — it runs once at Fortran-side setup before GPU dispatch);
`ANI-NONUNIFORM-REJECT` has no production fixture yet (do_cluster was found
to be transitively rejected by the pre-existing geometry check in the
common case) and rests on source-level fault-injection evidence only.

Eleven of the fourteen checklist items above are ticked because each is
individually true against delivered, tested evidence. The remaining three
are left open honestly rather than ticked on the strength of adjacent work:
"already-coarse unsafe blocks refine at an accepted synchronization point"
was not implemented or tested (the gate recomputes hard exclusion from
live atomistic state every step, but no fixture drives an already-coarse
block unsafe and confirms it refines); "selector diagnostics report
polarization and the reason for refinement" is not done (the gate reuses
the existing generic `'hard-atomistic-exclusion'` reason string, and no
polarization ratio is printed as its own diagnostic); and "unsupported
spatially varying anisotropy rejects before integration" has no production
trigger, matching the missing `ANI-NONUNIFORM-REJECT` exit evidence above.

This task was performed while RCG-02 remains open, per the blueprint's
exploratory-work allowance (section 8/10). Ticking individual items above
records what is honestly demonstrated; it is not a claim that RCG-03 is
closed. Closure additionally requires: RCG-02 to close, independent review
by someone other than the implementer, the three items above, HIP evidence,
and a production `ANI-NONUNIFORM-REJECT` fixture. This evidence is not
accepted until RCG-02 closes and the work is rebased and rerun.

**RCG-03 evidence (2026-08-08, CLOSED):** RCG-02 closed at `fae4c413`;
those commits sit chronologically on top of the RCG-03 patches above in a
single line of history, so no rebase was needed, only a rerun, which was
done first from fresh out-of-tree builds before any further edits
(`ctest -L cg13-cpu` 12/12, `ctest -L cg13-cuda` 15/15 on `d8b8c5ab`,
including the `mode:` rejection-matrix case that both RCG-01 and the
2026-08-06 RCG-03 evidence had flagged as a pre-existing unrelated
failure — it now passes cleanly and was not investigated further, being
out of scope). The three items left open above are now delivered:
already-coarse-block-refines-at-sync-point has an operator-level fixture
(`tests/coarse_graining/test_adaptive_hybrid_solver.f90`,
`test_polarization_forces_refine_of_coarse_block`) after finding that
production reconstruction structurally prevents an already-coarse block
from ever being observed polarization-unsafe through the ordinary step
loop (it rebuilds dormant atoms to exact full polarization every step,
before the gate runs); the polarization-forced-atomistic path now logs its
own `'polarization-unsafe'` reason (distinct from the static-mask-only
`'hard-atomistic-exclusion'`) with the triggering ratio, on both CPU and
CUDA (`evaluate_polarization_gate`/`evaluateAdaptivePolarizationGate`
gained a `block_ratio`/`polarizationRatio` diagnostic output, threaded
through the transition log and the GPU diagnostic snapshot/print); and
`ANI-NONUNIFORM-REJECT` now has a production fixture
(`tests/coarse_graining/e2e/ani_nonuniform_reject_cluster`, a `do_cluster`
embedding placed exactly on an existing host lattice site so
`clus_expand=0` and the geometry check the original candidate tripped
never fires), reached through the ordinary `sd.f95` executable and wired
into `run_setup_rejection_matrix.py` (31/31 cases) and the tracked-fixture
audit. Full detail, including two unrelated pre-existing `do_cluster`
setup-order quirks worked around entirely from the fixture's own input
files (no source changes), is in
`docs/RCG-03_POLARIZATION_ANISOTROPY_EVIDENCE.md`'s "RCG-03 closure
(2026-08-08)" section.

**RCG-03 is closed (2026-08-08, Human decision: Anders Bergman).** Two
items are explicitly deferred, matching RCG-02's precedent and for the
same reasons, not blocking: HIP execution evidence, because no HIP
toolchain or device exists in any environment used so far; and a separate
independent Opus/Terra or Sol adversarial review distinct from Human
approval, deferred to a later stage of the remediation program. Neither
deferral reflects a physics disagreement or an unresolved correctness
question — every fixture this document requires passes now on CPU and
CUDA.

---

### Task RCG-04: Replace vacuous e2e evidence with moving fixtures

**Dependencies:** RCG-02, RCG-03  
**Suggested primary:** Luna/Sonnet for harness and fixture construction; Sol
for production diagnostics  
**Suggested independent review:** Opus/Terra for physical oracle; Human for
accepted error budgets  
**Risk:** High validation breadth

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Redesign the production executable suite so its parity and coarse-dynamics
> claims depend on nonstationary states. Retain simple uniform fixtures only
> for zero-torque and setup smoke tests. Add deterministic moving fixtures
> for feature-off/all-fine, all-coarse, static mixed, adaptive mixed, and DMI
> interface behavior.
>
> Record state-sensitive observables rather than only direction sums and
> norms. Prove a nonzero initial torque and nonzero evolution. Compare full
> restart states or stable per-spin hashes, named fields and energies,
> transition histories, and analytic frequency/chirality data as applicable.
> Add negative controls that disable or sign-flip the component under test.
> Keep runtime modest enough for normal CI while separating longer validation
> cases if necessary.

#### Checklist

- [x] Uniform fixed-point cases are labelled as smoke/zero-torque tests only.
- [x] Feature-off/all-fine parity begins from the same moving state.
- [x] Initial torque exceeds a documented nontriviality floor.
- [x] Final displacement exceeds a documented nontriviality floor.
- [x] Complete state-sensitive trajectories agree within the fp64 budget.
- [x] All-coarse long-wave dynamics match an analytic or atomistic reference.
- [x] Static mixed work exercises atomistic, coarse, and interface ownership.
- [x] Adaptive mixed work performs accepted transitions during real motion.
- [x] DMI/anisotropy tests assert chirality and dynamics, not merely nonzero terms.
- [x] Interface error is measured under spatial refinement.
- [x] At least one wall or skyrmion crosses a static/adaptive boundary.
- [x] A broken coarse operator fails at least one fixture.
- [x] A reversed DMI sign fails the chiral fixture.
- [x] CPU/GPU parity uses identical initial data and observable definitions.
- [x] Every fixture and generator is tracked with provenance.
- [x] Precision-specific tolerances are justified from observed error scaling.

**Exit evidence:** `E2E-MOVING-OFF-FINE`,
`E2E-MOVING-ALL-COARSE`, `E2E-MOVING-STATIC`,
`E2E-MOVING-ADAPTIVE`, and `DMI-HYBRID-CROSSING`.

**RCG-04 evidence (2026-08-10, CLOSED):** `docs/RCG-04_MOVING_E2E_EVIDENCE.md`
records nine sliced sessions (RCG-04A-I: evidence contract, deterministic
moving-state generators, trajectory/observable infrastructure, then the five
exit-evidence packages and CPU/GPU precision-budget derivation) followed by
an independent RCG-04J closure audit. The closure audit confirmed the linear
RCG-04A-I ancestry, ran fresh out-of-tree configure/build/test on every
backend/precision available (CPU fp64, CUDA fp64, CUDA fp32 — HIP deferred,
no toolchain on any host used), and reconciled all sixteen checklist items
above against direct evidence rather than commit count; every fixture and
generator was independently confirmed git-tracked. All five required
exit-evidence packages pass fresh at every tested precision. Two real
production defects were found and fixed during RCG-04D/E (a missing physical
gyromagnetic-ratio factor in the AdaptiveCG atomistic step; a coarse
`channel_gamma` scaling gap), each with a defect-sensitivity negative control
proving the fix matters. The closure audit itself additionally found, and
disclosed rather than fixed (out of RCG-04's own scope): no CI workflow
currently runs any `cg13`/`moving-parity` test, and a pre-existing,
non-RCG-04 dipole-energy assertion (`gpu_fft_static_mixed`) fails
reproducibly at CUDA fp32.

**RCG-04 is closed (2026-08-10, Human decision: Anders Bergman).** Two items
requiring Human judgement were reviewed directly against their underlying
evidence and accepted: the RCG-04H DMI handedness/chirality convention
(`docs/RCG-04_MOVING_E2E_EVIDENCE.md` §16.1, consistent with the already-closed
RCG-02 derivation, with a reversed-sign negative control that correctly
fails) and the RCG-04I frozen fp64/fp32 precision budgets (§17.9, chosen with
headroom over observed error but well below each fixture's own physical
displacement, re-verified against negative controls after freezing). HIP
execution evidence remains deferred, matching RCG-02/RCG-03 precedent — no
HIP toolchain or device exists in any environment used so far. Rather than
leaving the remaining gaps as passive prose deferrals, five were promoted to
actively tracked, independently pickupable follow-up tasks (RCG-04-FU1
through RCG-04-FU5, defined immediately below); none blocks this closure or
RCG-05's start. No physics disagreement or unresolved correctness question
blocks this closure — every fixture this document requires passes now, on
every backend/precision available.

#### RCG-04 follow-up tasks (opened 2026-08-10, Human decision: Anders Bergman)

Non-blocking, independently pickupable in any order. Each may become its own
session; none is a dependency of RCG-05.

- **RCG-04-FU1 — HIP execution evidence.** Re-run the five RCG-04 exit-evidence
  packages and the RCG-04I backend-parity harness on real HIP hardware once a
  toolchain/device is available. Dependencies: HIP toolchain/hardware only.
- **RCG-04-FU2 — CI wiring.** Add a CI job (CPU-only is sufficient, matching
  `adaptive-cg-clean.yml`'s existing runner) that runs at minimum the five
  `adaptive-cg-moving-*` CTest targets on every push/PR; extend to
  `adaptive-cg-moving-backend-parity` if a CUDA-capable runner is ever
  available. Dependencies: none.
- **RCG-04-FU3 — RCG-04G refine-direction gap.** Find or construct a
  moving-wall geometry in which an accepted refine (coarse-to-atomistic)
  transition is genuinely demonstrated during motion, or determine whether
  `RECONSTRUCTION_CONE` succeeds where `RECONSTRUCTION_ALIGNED` was rejected
  in RCG-04G's own investigation. Dependencies: accepted RCG-04G.
- **RCG-04-FU4 — RCG-04E quantitative rate reconciliation.** Derive the
  coarse operator's precession-rate dependence on block size quantitatively
  and reconcile it against the atomistic reference, resolving whether the
  discrete-Laplacian dispersion argument (§13.6) fully explains the residual
  mismatch or whether an additional scale factor remains. Dependencies:
  accepted RCG-04E.
- **RCG-04-FU5 — fp32 `gpu_fft_static_mixed` failure.** Determine whether
  CUDA fp32 `EWALD3D_FFT` dipole coupling has a genuine precision-floor
  problem, or replace the fragile `coarse_dipole != 0` assertion (found by
  RCG-04J to test a quantity indistinguishable from floating-point noise at
  both precisions) with a physically meaningful check. Dependencies: none —
  pre-existing, non-RCG-04 harness (`run_production_e2e.py`).

---

### Task RCG-05: Restore CPU/GPU geometry and ownership equivalence

**Dependencies:** RCG-04  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra for ownership semantics; Luna/Sonnet for mask
and anisotropic-cell fixtures  
**Risk:** High cross-backend correctness

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Preserve the three directional interaction-derived buffer widths through
> the Fortran/C++ descriptor and CUDA/HIP kernels. Do not replace them with a
> scalar maximum unless the parent contract is explicitly changed and CPU
> uses the same conservative semantics. Add unequal-width orthogonal and skew
> fixtures with periodic wrapping and compare complete CPU/GPU fine, buffer,
> and coarse maps.
>
> Audit short-range, on-site, external, and uniform FFT dipole ownership.
> State why the uniform dipole field is all-grid and independent of
> `onsite_owner`, then add assertions that it is included exactly once.
> Document the monotonic pending-state invariant in dilation and remove any
> read/write race if sanitizer or reasoning cannot prove it safe.

#### Checklist

- [x] Buffer widths retain three directional components end to end.
- [x] CPU, CUDA, and HIP dilation semantics are identical.
- [ ] Non-cubic and skew-cell masks match exactly. **Not met as literally
      stated — see the RCG-05 closure decision below.** The buffer-width
      scalarization defect this task was chartered around is fixed and
      proven; full CPU/GPU map identity is still blocked by a separate,
      pre-existing GPU seed-mask-sourcing gap, and no skew fixture has run
      end to end (blocked by an unrelated `neighbourmap.f90` limitation).
      Left honestly unchecked rather than misstating the evidence.
- [x] Periodic wrapping is tested in every direction.
- [x] Every atomistic cross-interface bond is covered.
- [x] Dilation has no unproved read/write race.
- [x] Static and adaptive mask rebuilds preserve ownership invariants.
- [x] Short-range and on-site energies have non-overlapping owners.
- [x] Uniform FFT dipole ownership is documented independently of the mask.
- [x] Dipole field and energy are included exactly once.
- [x] An anisotropic-cell negative control detects scalarized buffer width.
- [x] CUDA/HIP descriptor layout checks cover the new vector data.

**Exit evidence:** `GEO-ANISO-BUFFER`, sanitizer/race evidence, and dipole
ownership derivative checks.

**RCG-05 evidence (2026-08-12, CLOSED):**
`docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md` records six sliced sessions
(RCG-05A-F: evidence contract and defect reproduction, skew/unequal-width
geometry generators, a CPU/GPU ownership-map comparator demonstrating the
buffer-width scalarization defect, the fix itself, a dilation-race
sanitizer audit and by-construction double-buffer fix, and dipole/
short-range/on-site ownership invariants) followed by an independent
RCG-05G closure audit. The closure audit confirmed the linear RCG-05A-F
ancestry, ran fresh out-of-tree configure/build/test on every
backend/precision available (CPU fp64, CUDA fp64, CUDA fp32 — HIP
deferred, no toolchain on any host used), and reconciled all twelve
checklist items above against direct evidence rather than commit count.
`GEO-ANISO-BUFFER`'s buffer-width-shape claim and the sanitizer/race
evidence are complete at every tested precision; the dipole ownership
derivative checks are complete for the tested (orthogonal, anisotropic)
fixture at every tested precision. The closure audit's own fresh fp32 run
additionally found, root-caused, and disclosed rather than fixed (out of
RCG-05's own scope): `testSelectorPolicyDescriptorLayout()` (RCG-05D) fails
under `UPPASD_PRECISION=SINGLE` CUDA builds — independently confirmed to be
a test-only float-literal-comparison bug, not a real descriptor-layout
defect — and reconfirmed the pre-existing, non-RCG-05 RCG-04-FU5
`gpu_fft_static_mixed` failure still reproduces at CUDA fp32.

**RCG-05 is closed (2026-08-12, Human decision: Anders Bergman).** One item
required Human judgement beyond what direct re-evidencing could resolve:
whether to accept `GEO-ANISO-BUFFER`'s narrower, fully-evidenced claim (the
buffer-width scalarization defect is fixed and proven, at every available
precision) as sufficient for this task's closure, given that the parent
checklist's literal "masks match exactly" claim is not met. Reviewed
directly against the underlying evidence
(`docs/RCG-05_GEOMETRY_OWNERSHIP_EVIDENCE.md` §13.4/§13.5/§13.7): accepted.
The two specific, independently-scoped gaps that block the literal claim —
a GPU seed-mask-sourcing defect unrelated to buffer-width dilation, and an
unrelated pre-existing `neighbourmap.f90` skew-cell limitation — are each
promoted to an active, independently tracked follow-up task rather than
left as passive deferral prose. Rather than leaving the remaining gaps as
passive prose deferrals, four were promoted to actively tracked,
independently pickupable follow-up tasks (RCG-05-FU1 through RCG-05-FU4,
defined immediately below); none blocks this closure or RCG-06's start. No
physics disagreement or unresolved correctness question blocks this
closure beyond the one item explicitly accepted above.

#### RCG-05 follow-up tasks (opened 2026-08-12, Human decision: Anders Bergman)

Non-blocking, independently pickupable in any order. Each may become its
own session; none is a dependency of RCG-06.

- **RCG-05-FU1 — HIP execution evidence.** Re-run `GEO-ANISO-BUFFER`, the
  sanitizer/race evidence, and the dipole ownership derivative checks on
  real HIP hardware once a toolchain/device is available. Dependencies:
  HIP toolchain/hardware only.
- **RCG-05-FU2 — GPU seed-mask-sourcing gap.** GPU's `hardAtomisticBlockMask`
  is sourced only from the polarization gate, not `cg_static_mask_file`,
  independently of and unaffected by the buffer-width fix. This is the
  direct blocker preventing full CPU/GPU ownership-map identity on
  `ownership_aniso_buffer` even after RCG-05D's fix (44/90 blocks still
  differ). Source GPU's hard mask from `cg_static_mask_file` as CPU already
  does, then re-run RCG-05C's comparator to confirm full map identity
  becomes achievable. Dependencies: accepted RCG-05C/D.
- **RCG-05-FU3 — skew-cell `neighbourmap.f90` gap.** No skew
  (non-orthogonal-lattice) fixture has ever run through the real
  executable's adaptive-mask or dipole path: `neighbourmap.f90` rejects the
  declared exchange bond with a "no basis match" error during setup, a
  pre-existing, unrelated Hamiltonian-setup limitation. Characterize and
  either fix or formally narrow `neighbourmap.f90`'s skew-cell
  exchange-shell tolerance/mapping so at least one skew fixture can run end
  to end. Dependencies: none — pre-existing, non-RCG-05 limitation.
- **RCG-05-FU4 — `testSelectorPolicyDescriptorLayout()` fp32 test bug.**
  This host-only CTest target fails reproducibly under
  `UPPASD_PRECISION=SINGLE` CUDA builds due to an exact-equality comparison
  against a hardcoded `double` literal that a `float`-valued `real` cannot
  round-trip to (`tests/coarse_graining/test_gpu_adaptive_runtime.cpp:637-640`).
  Independently confirmed not a real descriptor-layout defect. Compare
  `copied.*Threshold` against the pre-copy `policy.*Threshold` value (or use
  a numeric tolerance), not a hardcoded double literal. Dependencies: none
  — pre-existing test-only defect, introduced by RCG-05D, first exercised
  by RCG-05G.

---

### Task RCG-06: Remove memory, allocation, timing, and precision hazards

**Dependencies:** RCG-05  
**Suggested primary:** Sol  
**Suggested review:** Luna/Sonnet for allocation/lifecycle tests; Opus/Terra
for precision budget  
**Risk:** Medium-high implementation breadth

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Move full-system automatic arrays and per-timestep temporaries into
> persistent, explicitly owned workspace with allocation preflight, cleanup,
> and memory diagnostics. Preserve multi-ensemble indexing and avoid changing
> operator mathematics. Add a representative large-host fixture under a
> documented stack limit.
>
> Accumulate global energy terms in double precision regardless of the field
> storage precision, while retaining explicit FP32 field/output conversion.
> Correct phase timing so wall-clock measurements include both field
> evaluations, integration stages, restriction/reconstruction, selector,
> compaction, and waits. Report unaccounted wall time.
>
> Measure reconstruction RNG spatial statistics. Replace the generator only
> if the accepted nonzero-cone model needs stronger independence; otherwise
> document why zero-cone production is unaffected and leave nonzero-cone
> acceptance open.

#### Checklist

- [ ] No routine places `O(n_atoms)` adaptive workspace on the stack.
- [ ] Per-step full-system allocations are hoisted or proven negligible.
- [ ] Workspace allocation, resize, reuse, and cleanup are tested.
- [ ] Memory preflight includes all persistent CPU/GPU adaptive workspace.
- [ ] A large-host fixture passes with the documented stack limit.
- [ ] GPU global energy accumulation uses FP64 storage and arithmetic.
- [ ] FP32 field parity and FP64 energy budgets are distinct.
- [ ] Energy error scaling is measured over increasing system size.
- [ ] Timers use a suitable wall/device clock.
- [ ] Both Heun field evaluations are included.
- [ ] Phase totals plus unaccounted time reconcile with external wall time.
- [ ] Multiple ensembles retain correct indexing.
- [ ] RNG correlation evidence is recorded and its scope decision is explicit.
- [ ] Existing derivative and moving-state fixtures remain unchanged.

**Exit evidence:** `MEM-LARGE-HOST`, `ENERGY-FP32-ACCUM`, allocation
lifecycle tests, and reconciled timing output.

---

### Task RCG-07: Repair CPU algorithmic scaling

**Dependencies:** RCG-06  
**Suggested primary:** Sol  
**Suggested review:** Opus/Terra for algorithm; Luna/Sonnet for scaling harness  
**Risk:** High performance change to accepted CPU reference

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Replace all-block-by-all-seed dilation with the local directional stencil
> already represented by the interaction-derived widths. Preserve periodic
> semantics and exact ownership maps. Profile the full CPU adaptive step over
> increasing block counts and coarse fractions.
>
> Add OpenMP or an equivalent accepted CPU parallel strategy to independent
> restriction, prolongation, Hamiltonian, and integration work where
> deterministic reduction and energy ownership permit it. Reduce unnecessary
> full-system passes or compact atomistic work when doing so does not alter
> the accepted adjoint/interface model. Do not trade correctness for an
> asymptotic claim.

#### Checklist

- [ ] CPU dilation visits only the local stencil around relevant seeds.
- [ ] Dilation work scales linearly in blocks for fixed physical width.
- [ ] New dilation maps are bitwise identical to the accepted reference.
- [ ] Parallel loops are race-free and deterministic within stated budgets.
- [ ] Energy reductions use stable precision and ownership.
- [ ] One-thread results match the accepted pre-optimization reference.
- [ ] Multi-thread scaling is reported with thread affinity/configuration.
- [ ] Coarse-fraction sweeps separate field, integration, and rebuild costs.
- [ ] Remaining unavoidable `O(N)` passes are documented.
- [ ] Feature-off CPU performance remains unchanged within noise.
- [ ] Moving and derivative fixtures pass after optimization.

**Exit evidence:** CPU ownership equivalence, complexity sweep, thread-scaling
report, and feature-off control.

---

### Task RCG-08: Parallelize the adaptive GPU production path

**Dependencies:** RCG-07  
**Suggested primary:** Sol  
**Required independent review:** Sol or Opus/Terra not responsible for the
kernel patch; Luna/Sonnet for device parity harness  
**Risk:** Very high

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Replace the single-threaded restriction, publication, atomistic
> Hamiltonian, energy finalization, predictor, and corrector reference kernels
> with parallel CUDA/HIP implementations. Preserve unique-pair ownership,
> normalized-projection adjoints, deterministic transition publication,
> multi-ensemble indexing, and the accepted Heun update.
>
> Use compact active lists where they reduce real Hamiltonian and integration
> work. Replace the `O(N log N)` scan with an accepted linear-work backend
> primitive or an equivalent portable implementation, accounting for
> temporary storage in preflight. Keep a serial reference available in tests
> if useful, but do not dispatch it as the production hot path.
>
> Establish FP64 parity before FP32 tuning. Run CUDA/HIP sanitizers and report
> occupancy, launch counts, phase times, and remaining serial sections.

#### Checklist

- [ ] No production hot kernel is guarded to a single device thread.
- [ ] Restriction covers all active atoms/blocks in parallel.
- [ ] Atomistic unique-pair fields and energies are race-safe.
- [ ] Predictor and corrector integrate compact active work in parallel.
- [ ] Publication preserves deterministic accepted state transitions.
- [ ] Energy finalization uses the accepted FP64 reduction.
- [ ] Compaction performs linear work and uses tracked temporary storage.
- [ ] Launch and synchronization counts are reported per phase.
- [ ] FP64 fields, energies, trajectories, and decisions match CPU references.
- [ ] FP32 has a separate accepted error budget.
- [ ] Multiple ensembles and anisotropic buffer descriptors pass.
- [ ] CUDA sanitizer reports no memory or race errors where available.
- [ ] HIP execution and sanitizer evidence are recorded where available.
- [ ] The serial reference, if retained, is test-only and clearly labelled.
- [ ] Feature-off device inventory and timing remain unchanged within noise.

**Exit evidence:** device parity suite, sanitizer logs, kernel/launch profile,
and compaction complexity evidence.

---

### Task RCG-09: Establish an honest production performance result

**Dependencies:** RCG-07, RCG-08  
**Suggested primary:** Sol for measurement; Luna/Sonnet for harness  
**Required review:** Opus/Terra for methodology; Human for any production or
speedup wording  
**Risk:** Medium code risk, high claim risk

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Benchmark adaptive coarse graining against UppASD's actual feature-off
> atomistic Hamiltonian and integrator using identical physics, precision,
> geometry, ensemble count, timestep, device, and measurement scope. Retain
> the synthetic FMA control only as an inactive-runtime microbenchmark and
> rename it accordingly.
>
> Sweep fine fractions, block sizes, system sizes, and representative
> textures. Report active atoms, active blocks, interface work, short-range
> bonds, selector/compaction cost, restriction/reconstruction, both field
> evaluations, integration, FFT when enabled, host waits, memory, median,
> dispersion, and raw samples. Separate setup from steady state.
>
> A crossover is accepted only when adaptive total wall time becomes lower
> than the real atomistic production baseline by more than the uncertainty
> margin while all correctness fixtures pass.

#### Checklist

- [ ] The baseline calls the real production atomistic path.
- [ ] Baseline and adaptive cases use identical supported physics.
- [ ] Synthetic inactive-overhead control is clearly renamed and separated.
- [ ] Setup, warm-up, and steady-state measurement scopes are explicit.
- [ ] CPU and GPU use appropriate wall/device timing.
- [ ] Active atom, block, interface, and bond fractions are reported.
- [ ] Every relevant phase and host wait is reported.
- [ ] Median, MAD or equivalent dispersion, repetitions, and raw samples exist.
- [ ] Device, driver, compiler, flags, clocks, and precision are recorded.
- [ ] CUDA and HIP results are separate.
- [ ] CPU thread count and affinity are recorded.
- [ ] Crossover is computed against the production atomistic baseline.
- [ ] Speedup exceeds its uncertainty margin.
- [ ] No result is extrapolated to an unavailable backend.
- [ ] All physics/parity fixtures pass on the measured commit.
- [ ] Human approves any restored performance wording.

**Exit evidence:** `PERF-ATOMISTIC-PROD`, `PERF-CG-SWEEP`, raw benchmark
artifacts, and independent methodology review.

---

### Task RCG-10: Reconcile release evidence with the parent blueprint

**Dependencies:** RCG-00 through RCG-09  
**Suggested primary:** Luna/Sonnet for evidence audit and documentation  
**Suggested review:** Opus/Terra for adversarial cross-task review; Sol for
implementation audit; Human for final acceptance  
**Risk:** High governance importance

#### Prompt

> Before editing, read `docs/ADAPTIVE_COARSE_GRAINING_REMEDIATION_BLUEPRINT.md`
> in full and read the parent
> `docs/ADAPTIVE_COARSE_GRAINING_BLUEPRINT.md` in full. Treat their evidence
> policy, dependencies, scope boundaries, and acceptance gates as governing
> context for this task.
>
> Starting from a clean clone of the candidate commit, run the complete CPU,
> CUDA, and HIP matrices where available. Audit every parent blueprint box
> provisionally reopened by this document. For each item, link the exact
> remediation task, fixture, command, result, and reviewer. Remove or correct
> stale evidence paragraphs rather than stacking contradictory claims.
>
> Work through every cross-task review question in parent section 12. A box is
> checked only after the reviewer can cite current evidence. Distinguish
> unavailable backend evidence from passing evidence. If meaningful
> production speedup has not been demonstrated, amend the supported claim to
> a correctness prototype and leave the parent GPU production/success boxes
> open.
>
> Preserve this remediation blueprint as the audit trail. Do not delete it
> after reconciliation.

#### Checklist

- [ ] Candidate validation begins from a clean clone.
- [ ] Repository status contains no required untracked files.
- [ ] CPU release matrix passes.
- [ ] CUDA release matrix passes on available hardware.
- [ ] HIP release matrix passes on available hardware.
- [ ] Unavailable backends are explicitly marked unavailable.
- [ ] Legacy supported builds pass.
- [ ] Every reopened parent item has current linked evidence.
- [ ] Superseded or false evidence wording is corrected.
- [ ] Parent section 12 physics review is completed.
- [ ] Parent section 12 numerical review is completed.
- [ ] Parent section 12 software review is completed.
- [ ] Parent section 12 pull-request review is completed.
- [ ] Performance wording matches RCG-09 evidence.
- [ ] Model limitations and rejected combinations remain prominent.
- [ ] No finite-temperature, explicit-device, or general multi-channel scope leaked in.
- [ ] Independent adversarial review reports no unresolved blocker.
- [ ] Human approves physics conventions, model error, and release wording.
- [ ] The parent definition of success is honestly accepted or explicitly amended.

**Exit evidence:** Updated parent blueprint and release-validation document,
clean-clone validation matrix, independent review report, and Human sign-off.

---

## 10. Sequential session protocol

Each implementation session should use this handoff:

1. Select the first incomplete `RCG-*` task whose dependencies are closed.
2. Read the parent blueprint, this blueprint, dependency evidence, and current
   implementation before editing.
3. Reproduce the relevant finding on the current commit.
4. Add the discriminating fixture or negative control before or with the fix.
5. Keep the patch within the selected task.
6. Run focused tests, then the required clean gate.
7. Record commands, environment, raw results, and remaining limitations.
8. Obtain the specified independent review.
9. Check only evidence-backed boxes in this document.
10. Stop and hand off the next task; do not opportunistically tick parent
    release boxes.

Recommended session order:

```text
RCG-00 -> RCG-01 -> RCG-02 -> RCG-03 -> RCG-04 -> RCG-05
       -> RCG-06 -> RCG-07 -> RCG-08 -> RCG-09 -> RCG-10
```

Do not accept tasks out of order. If exploratory work for a later task is
performed early, rebase it onto the accepted dependency commit and rerun all
required fixtures during its own session.

### Session evidence template

```text
Task:
Commit:
Finding reproduced:
Files intentionally changed:
Negative control:
Focused test command/result:
Clean build command/result:
Required suite command/result:
Backend/precision/device:
Performance raw artifact:
Independent reviewer:
Checklist items accepted:
Open limitations:
Next permitted task:
```

---

## 11. Cross-remediation review checklists

### Physics review

- [ ] One DMI convention governs CPU, CUDA, HIP, extraction, coarse, and hybrid paths.
- [ ] The analytic dimer and small-\(q\) signs agree.
- [ ] Every field still has a written energy and derivative test.
- [ ] Polarization prevents unsupported single-channel cancellation.
- [ ] Anisotropy is represented locally or rejected, never sampled silently.
- [ ] Exchange tensor assumptions needed for accuracy are asserted.
- [ ] Dipole ownership remains exactly once and independent of short-range masks.
- [ ] Moving interface fixtures preserve the intended texture within budget.

### Numerical review

- [ ] End-to-end states are demonstrably nonstationary where dynamics is claimed.
- [ ] Complete state-sensitive observables replace aggregate-only parity.
- [ ] CPU/GPU directional buffer masks match.
- [ ] Large host runs do not depend on stack size.
- [ ] Global energy accumulation meets FP32 and FP64 budgets.
- [ ] Transition, publication, and compaction operations are race-free.
- [ ] Error decreases under spatial refinement.
- [ ] Timing phases reconcile with measured wall time.

### Software review

- [ ] Clean builds cannot consume stale Fortran modules.
- [ ] Every supported build system contains the complete dependency graph.
- [ ] Every harness dependency is tracked.
- [ ] Persistent workspaces have explicit ownership and cleanup.
- [ ] Feature-off execution and inventory remain isolated.
- [ ] CPU and GPU descriptors have identical semantics.
- [ ] Tests contain defect-sensitive negative controls.
- [ ] CI distinguishes unavailable backends from passes.

### Performance review

- [ ] No production hot kernel is intentionally single-threaded.
- [ ] CPU dilation is local rather than quadratic.
- [ ] Compaction uses linear work.
- [ ] The benchmark baseline is UppASD's real atomistic path.
- [ ] Same-physics comparisons include the complete timestep.
- [ ] Active work, overhead, uncertainty, and raw measurements are reported.
- [ ] Speedup wording does not exceed observed hardware evidence.

### Pull-request review

- [ ] The patch implements one declared `RCG-*` task.
- [ ] Dependencies were accepted before implementation.
- [ ] Finding reproduction and negative control are attached.
- [ ] Sign/unit changes have independent physics review.
- [ ] Performance changes retain correctness references.
- [ ] Baseline and introduced failures are separated.
- [ ] Untracked local files do not contribute to acceptance.
- [ ] Parent boxes are not ticked early.

---

## 12. Risk register

| Risk | Consequence | Mitigation |
| --- | --- | --- |
| A stale build appears to pass | False acceptance survives another review | Fresh out-of-tree builds and clean CI are mandatory |
| DMI is fixed only in the coarse path | CPU hybrid remains chirally inconsistent | One analytic dimer oracle covers every active path |
| A self-consistent validator repeats the implementation sign | Wrong physics passes internally | Derive the oracle independently from the written Hamiltonian |
| Low-polarization blocks are normalized | Plausible-looking noise-directed coarse states | Non-overridable polarization safety gate |
| Central-cell anisotropy remains implicit | Heterogeneous systems are silently homogenized | Block-local construction or setup-time rejection |
| Uniform e2e states remain the parity oracle | Broken operators continue to pass | Nonzero torque/displacement and defect-sensitive controls |
| CPU/GPU halo shapes differ | Backend-dependent ownership and trajectories | Preserve directional widths in one semantic descriptor |
| Workspace refactor changes adjoint mathematics | Correct projection is damaged during optimization | Retain derivative and adjoint fixtures throughout |
| GPU parallelization changes pair ownership | Double-counted or missing fields/energies | Unique-pair analytic and CPU/GPU parity fixtures |
| FP64 accumulation hides FP32 field error | Precision claim becomes ambiguous | Separate field and accumulator budgets |
| Synthetic control is called a production baseline | Misleading speedup claim | Rename it and compare with real UppASD execution |
| CUDA evidence is generalized to HIP | Unsupported portability claim | Backend-specific matrices and explicit unavailable state |
| Remediation document becomes a second permanent specification | Contracts diverge | Parent remains authoritative; RCG-10 reconciles and closes |
| Scope expands during repair | Review becomes unfinishable | Reject new physics and preserve original out-of-scope boundary |

---

## 13. Definition of remediation success

Remediation is complete only if:

1. tracked source builds cleanly without stale artifacts through every
   supported build path;
2. the complete referenced fixture set is committed and CI-executable;
3. DMI has one Human-approved convention across CPU, CUDA, HIP, extraction,
   coarse, and hybrid implementations;
4. low-polarization and unsupported anisotropy inputs cannot silently enter an
   invalid coarse model;
5. executable validation contains deterministic moving states and
   defect-sensitive oracles;
6. CPU/GPU geometry, ownership, fields, energies, transitions, and
   trajectories match within justified precision budgets;
7. large systems do not rely on oversized automatic arrays or per-step
   allocation churn;
8. the production CPU and GPU algorithms expose parallel work and avoid the
   identified quadratic or `O(N log N)` bottlenecks where required;
9. any speedup claim compares with the real UppASD atomistic production path
   under identical physics and exceeds measurement uncertainty;
10. an independent reviewer and Human owner approve the restored evidence;
11. the parent blueprint and release documentation contain no stale or
    contradictory acceptance claims;
12. the parent section 12 checklists are completed from current evidence.

If items 1 through 7 and 10 through 12 pass but item 9 does not, the result may
be accepted only as a correctness/reference implementation. The parent
blueprint must then leave GPU production-performance gates open and remove
speedup wording. That is an honest partial outcome, not a production release.

---

## 14. Final handback to the parent blueprint

RCG-10 performs the only authorized final handback:

1. update parent evidence paragraphs with remediation commits and fixtures;
2. retick provisionally reopened task boxes only where current evidence exists;
3. complete the parent cross-task review checklists;
4. update `docs/CG-13_RELEASE_VALIDATION.md` with the final backend matrix,
   fixture provenance, precision budgets, and benchmark method;
5. record Human physics and release approval;
6. retain this document as the review-to-resolution audit trail.

Until that handback is complete, the adaptive coarse-graining branch should
not be described as release-accepted, DMI-conformant, CPU/GPU-parity accepted,
or performance-production-ready.
