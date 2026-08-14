# Open questions on the adaptive atomistic kernel

**Raised by:** Anders Bergman, 2026-08-14, during RCG-09 review
**Status:** Promoted to **task RCG-09A** in the remediation blueprint on
2026-08-14 (Human decision, Anders Bergman). This note remains the background
record; the acceptance criteria live in the task block. **Not investigated,
adjudicated, or acted on by the RCG-09 session** — RCG-09 is a measurement task
and physics-kernel refactoring is outside its scope.
**Audience:** whoever runs RCG-09A

---

## 1. Where this came from

RCG-09 measured the adaptive GPU path against UppASD's production atomistic
path (`docs/RCG-09_PRODUCTION_PERFORMANCE_EVIDENCE.md`). In explaining why the
all-fine adaptive case is *not* simply "atomistic plus overhead", that document
notes two structural differences between the two implementations
(§3.5.1 there):

* `heisge` walks a per-atom neighbour list and evaluates each pair twice, once
  from each endpoint; the adaptive path walks a compact unique-bond list,
  evaluates each pair once, and scatters to both endpoints.
* The production step runs Depondt with `thermfield.randomize()` every step;
  the adaptive step runs a deterministic Heun update with no stochastic field.

Both observations were made to explain a *timing* difference. Anders raised
that both of them are really physics questions, and that they point at the same
conclusion: the adaptive atomistic part arguably ought to call `heisge`, or at
minimum follow a demonstrably equivalent workflow.

---

## 2. Question 1 — does unique-pair scatter hold for antisymmetric couplings?

**The observation.** Evaluating each pair once and scattering to both endpoints
is obviously sound for Heisenberg exchange, where the coupling is symmetric. It
is not obviously sound for DMI, where `D_ij = -D_ji`. Halving the pair
arithmetic is only legitimate if the two endpoints receive contributions that
respect that antisymmetry.

**Where to look.**

* `source/gpu_files/gpuAdaptiveRuntime.cpp`, `evaluateAdaptiveAtomisticBonds`
  (around line 778). The kernel forms two quantities per bond — one contracting
  the bond matrix with `sj`, one contracting it with `si` over the transposed
  index order — and scatters one to each endpoint.
* `source/CoarseGraining/adaptivecgproduction.f90`, `build_unique_bonds`
  (around line 1358). This is where the per-bond `3x3` matrix is assembled; the
  DMI branch builds a cross-product matrix and applies an orientation sign
  keyed on `atom < target`.
* `source/gpu_files/gpuHamiltonianCalculations.cpp`, `Heisge::each`, for the
  production convention it must agree with.

**What a reviewer would need to establish.** Whether the transpose contraction
plus the orientation convention together reproduce `D_ij = -D_ji` exactly, for
both endpoint orderings, and whether that is covered by an executable oracle
rather than by inspection. Note that `docs/RCG-02_DMI_HANDEDNESS_EVIDENCE.md`
already discusses a "unique-pair antisymmetric matrix", so there may be existing
evidence bearing directly on this — it should be read before any new work
starts. RCG-09 did not evaluate whether that evidence is sufficient.

**Note on RCG-09's own parity check.** RCG-09's benchmark verifies that the
adaptive and production paths agree on the field to `2.8e-16`. That check runs
on an **isotropic-exchange-only** fixture with `do_dm` disabled, so it exercises
exactly the symmetric case where unique-pair scatter is uncontroversial. It is
not evidence about DMI and should not be cited as such.

---

## 3. Question 2 — thermal fields in the adaptive atomistic part

**The observation.** The atomistic region of an adaptive run is ordinary
atomistic spin dynamics. At finite temperature it should carry the same
stochastic field that ordinary UppASD applies there. The adaptive path today
generates none.

**Current state, for context.** This is a declared capability boundary rather
than an omission: `adaptivecgproduction.f90` (around line 829) rejects setup
outright when `Temp /= 0`, alongside `do_qhb` and `do_3tm`, with the diagnostic
"adaptive coarse graining requires deterministic T=0 dynamics", and setup prints
`capability accepted: regular periodic single-FM deterministic Heun`. The parent
blueprint lists finite-temperature coarse dynamics as out of scope, and the
remediation blueprint repeats that exclusion.

**What is worth a reviewer's attention anyway.**

* Whether zero-temperature-only remains the intended long-term boundary, or
  whether thermal fields in the *atomistic* region specifically are a nearer-term
  requirement than full finite-temperature coarse dynamics. Those are different
  scopes and the current rejection collapses them into one.
* If thermal fields are added, RCG-09's recorded RNG asymmetry (§4.4 of the
  evidence document) disappears, and the performance comparison changes in the
  production baseline's favour. Any future crossover claim would need
  re-measuring under that change, not adjusting.
* Reconstruction of coarse blocks at finite temperature is a separate and harder
  problem, tracked at RCG-06-FU2 (nonzero-cone acceptance). The two should not
  be conflated.

---

## 4. The underlying suggestion

Anders' preferred resolution, recorded as stated rather than assessed:

> The adaptive atomistic part should use `heisge`, or at least a fully
> equivalent workflow. The latter would however be harder to maintain, so that
> is just a plan B.

Calling `heisge` directly would make both questions above moot by construction —
one convention, one implementation, no equivalence argument to maintain. It
would also change the performance picture: the pair arithmetic would no longer
be halved, and the compact-active-list structure the adaptive path relies on
would need to be reconciled with `heisge`'s per-atom neighbour-list traversal.

RCG-09 takes no position on the trade-off and did not attempt to estimate its
cost.

---

## 5. Relationship to open follow-ups

* **RCG-09B** (which supersedes RCG-09-FU1) rewrites the energy accumulation
  inside `evaluateAdaptiveAtomisticBonds` and the coarse kernels. It depends on
  RCG-09A precisely because this note may change what that kernel is: the two
  must not run concurrently.
* **RCG-02** owns the DMI convention and may already answer §2.
* **RCG-06-FU2** owns nonzero-cone reconstruction, adjacent to §3.
