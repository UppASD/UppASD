# Handover — V3 + M1 (graceful OOM and upfront memory budget)

Written 2026-07-20 at commit `5d46f02`, branch `gpu_hip_cu_ab`.

Everything below was checked against the tree, not inferred from commit titles or checklist
ticks. Line numbers drift — re-locate by symbol.

---

## Read this first: V3 is not done, and M1 is written on top of it

M1's task text opens: *"V3 delivers the catch + report. Extend with a pre-allocation
estimate."* **V3 does not exist.** It was ticked "verified landed 2026-07-20" by `28f4b5a`, the
same reconciliation pass that over-ticked PR1–PR7. I un-ticked it in `5d46f02` with the
evidence; the short version:

- `fort_helper.cpp:42-56` — the four `extern "C"` entry points (`gpusim_initiateconstants_`,
  `gpusim_initiatematrices_`, `gpusim_gpurunsimulation_`, `gpusim_release_`) are bare one-line
  forwarders. No `try`, no `catch`.
- `grep -rn "MemGetInfo" source/gpu_files/` returns **nothing**. Part 2 of V3 (the post-init
  free/total summary) is absent too.
- `git log -S "catch (const std::exception" -- source/gpu_files/fort_helper.cpp` returns no
  commit. It never landed and was never reverted — it was simply never written.

The premise is intact: `ASSERT_GPU` (`base.hpp:32-36`) does `throw std::runtime_error`, and
nothing between there and Fortran catches it. So today an allocation failure propagates through
`extern "C"` and terminates with a bare `terminate called after throwing...`.

**Do V3 first, as its own commit.** It is P0, small, and M1's diagnostics are worthless without
it. Then M1.

**General warning:** treat anything ticked by `28f4b5a` as unverified. Two of its ticks have now
been found wrong (PR1–PR7, V3). Re-check against the tree before building on one.

---

## Source documents

- **`uppasd_gpu_implementation_checklist.md`** — the authority. V3 is in Phase 0; M1 is in
  Phase 2 (Memory). Read both entries in full before starting.
- **`uppasd_gpu_audit.md`** — the original code audit. **M1 and V3 have no audit item.** The
  audit's sections are A1–A5 (architecture), C1–C10 (correctness), P1–P7 (performance), I1–I3
  (interop), B1–B2 (build), plus hygiene. Nothing in it covers OOM or memory budgeting. Do not
  go looking for one; if you cite an audit finding in a commit or changelog, cite one that
  exists. (The previous handover was burned by exactly this — the checklist cites a broken
  `MATCHES "^[__86+PTX]"` regex from audit B2 that is not in the tree.)
- **`handover_cm_pr5.md`** — the previous handover. Its Part 1 (CM1–CM6) is now largely done;
  its conventions section is superseded by the one below where they differ.

---

## Repo conventions

### Commit messages: subject line only

`area: imperative summary`, <= ~72 chars, **empty prose body**, plus the harness-required
trailer. Prefixes in use: `gpu:`, `cmake:`, `tests:`/`test:`, `docs:`, `fortran:`.

```
gpu: compute versine stably in Depondt rotation
cmake: single source of truth for GPU sources
docs: un-tick V3, it was never implemented
```

```
Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
```

Two commits (`0c4179d`, `28f4b5a`) have long prose bodies. That was a departure, not a
standard — do not copy them. Reasoning goes in the checklist entry or in a handover file.

### Other conventions

- **Put the reasoning in the checklist.** When you finish an item, tick it *with the evidence
  you actually gathered* — commands run, values observed — not a general impression. That
  discipline is what caught the V3 and PR1–PR7 over-ticks.
- **Report honestly when a change measures as nothing.** PR5(a)'s versine fix produced zero
  change in double precision; the checklist says so rather than claiming an accuracy win.
- Mixed CRLF/LF exists in the tree (`gpuSDSimulation.cpp` notably). Editing tools may normalize
  whole regions and produce a huge spurious diff. Run `git diff --stat` after every edit.
- Stage only files you intend to change. The `tests/*` edits that were floating uncommitted for
  most of the CM work were the user's, and they committed them themselves in `27d9382`; the
  tree is clean now, but the habit stands.
- `git diff --check` before committing, scoped to your files.

---

## Build and validation

**Presets are now the supported path** (CM5). The hand-configured `build_ab`, `build_fastcopy`,
`build_ptds` still exist and still work, but prefer presets for new trees.

```sh
cmake --preset cuda-debug     && cmake --build --preset cuda-debug     --parallel 8
cmake --preset cuda-fastcopy  && cmake --build --preset cuda-fastcopy  --parallel 8
cmake --preset cuda-ptds      && cmake --build --preset cuda-ptds      --parallel 8
cmake --list-presets          # full list
```

`cmake --preset cuda-relwithdebinfo` is the one to use for profiling (optimised, `-lineinfo`).
Build types are now standard CMake names — `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`,
`Coverage` — accepted case-insensitively; `TESTING` is a deprecated alias of `Coverage`.

```sh
python3 tests/gpu_regression/run.py --binary build/cuda-debug/bin/sd.f95.cuda
python3 tests/gpu_regression/sanitize.py --binary build/cuda-debug/bin/sd.f95.cuda
```

**Expected: 11 PASS, 1 FAIL in every configuration.** The failure is `sc_dm_uniaxial` — see
below, it is directly relevant to M1. `kagome_tensor` reports `abs=1.00e-11, rel=3.62e-09`
rather than exact zero; that is the sensitive case to watch. Everything else compares exactly
`0.00e+00`.

`ctest` covers the CPU side and passes 2/2 (`regression-test`, `asd-tests`). **The GPU cases in
ctest are broken on obsolete reference values** — per the maintainer. Do not use ctest to judge
GPU work; `tests/gpu_regression` is the trustworthy harness.

`cuda-ptds` is the regression net for stream-ordering bugs — it removes the implicit edge
between the legacy default stream and the blocking `workStream`. Keep it passing.

---

## Evidence gathered for M1 (all confirmed present, 2026-07-20)

### The allocation sites are in three places, not one

M1 says *"keep the estimator as one function so it stays in sync with allocations"*. That is in
tension with the current layout — the device allocations are spread across three files:

| Site | What it allocates |
|---|---|
| `gpuSimulation.cpp:279-358` (`initiateMatrices`) | Hamiltonian (`aHam`, `ncoup`/`j_tensor`, `nlist`, `nlistsize`, `dmvect`/`dmlist`/`dmlistsize`, aniso `kaniso`/`eaniso`/`taniso`/`sb`), lattice (`extfield`, `beff`, `b2eff`, `eneff`, `emomM`, `emom`, `emom2`, `mmom`, `mmom0`, `mmom2`, `mmomi`), `energyM`, `btorque` |
| `gpuLatticeConvolutionHamiltonian.cpp:377-387, 501, 561` | `field_real`, `spin_fft`, `field_fft`, `kernel_fft`, `kernel_real`, `cell_vectors`, `basis_positions`, `interaction_rij` |
| `correlations/gpuCorrelations.cpp:87-97` | `r_mid`, `q`, `coord`, `dt` (`sc_max_nstep`-sized) |

The estimator must cover all three to be meaningful — M1's task explicitly lists the
convolution descriptor and `nq, nw, sc_max_nstep`. Decide early whether the estimator is one
function that reaches into all three, or a small per-component contract (e.g. each component
exposes `static size_t estimateBytes(...)`) summed in one place. The second keeps it in sync
more honestly; the first matches the task text more literally. **Say which you chose and why.**

### The current check is post-hoc, which is the whole point of M1

`initiateMatrices` allocates everything and *then* checks `gpuHasNoData()` at
`gpuSimulation.cpp:373-379`, reporting via `GPU_GET_ERROR_STRING(GPU_GET_LAST_ERROR())`. By
then the allocation has already failed. M1 replaces this with an estimate computed *before* the
first `Allocate`.

Drive-by while you are in there: the message at `:377` reads `"Gpu: Failed to allocate
meFlags.mory: %s"` — a botched find/replace of `mFlags` into the word "memory". Fix it in the
same commit; it is in the exact code path M1 rewrites.

### `TensorMemoryTracker` has no query API — you will have to add one

`measurement/memoryMeasurement.h:14-39`. It exposes only `printResults()`, `saveToFile()`,
`add_host/add_device/remove_host/remove_device`. The peak is computed *inside* `printResults()`
as `*std::max_element(device_memory)` and never returned. M1's acceptance is "estimator matches
tracker within 5 %", which needs something like `static int64_t peak_device()`. Add the getter
alongside; it is a two-line change and keeps the acceptance testable rather than eyeballed.

Note the members are `std::vector<int64_t>` under a mutex, seeded with `{0}` — a running series,
not a scalar. Read it before assuming semantics.

### `do_meminfo` is a name collision — do not wire M1 to it

`do_meminfo` (`Tools/profiling.f90:11`, `Input/inputhandler.f90:1103`) is the **Fortran-side**
CPU allocation logger and has nothing to do with the GPU tracker, which is printed
unconditionally at `gpuSimulation.cpp:481`. `tests/gpu_regression/cases.json` sets
`do_meminfo: 0` for `sc_dm_uniaxial`. If M1 needs a runtime switch, introduce a clearly
distinct name rather than overloading this one.

### `sc_dm_uniaxial` fails *inside the function M1 rewrites*

The one persistently failing regression case segfaults in `gpusim_initiatematrices` ←
`sd_mphasegpu` (`sd_driver.f90:1185`), in **all** configurations, and predates all recent work.
That is `initiateMatrices` — precisely M1's target.

This is an opportunity, not a warning: **V3's catch handler may well turn that bare segfault
into a readable diagnosis**, and the budget check may reveal it as a memory problem (or rule
that out). Run this case deliberately once V3 lands, before writing the estimator. If it turns
out to be an OOM, M1 has a real acceptance case handed to it. If it is an unrelated bug, say so
and leave it — do not let it silently expand M1's scope.

---

## Ordering

```
V3 (catch + report)  ->  M1 (estimator + budget abort)
```

Within M1: tracker getter -> estimator function(s) -> budget comparison + abort table -> the
5 % agreement test.

---

## Risks

1. **Don't make the estimator a second source of truth that silently rots.** This is the
   failure mode the CM4 work just cleaned up in CMake (a duplicated source list that drifted
   until it referenced a file that did not exist). An estimator that lists allocations
   separately from the code that performs them will drift the same way. The 5 % agreement test
   is the guard — make it run somewhere automatic, not as a one-off you did by hand once.
2. **`MemGetInfo` free bytes are not what you can actually allocate.** Fragmentation, context
   overhead and other processes all eat into it. The ~90 % threshold in the task is a heuristic;
   treat a near-miss as advisory, and make sure the abort message says the estimate is
   approximate rather than implying precision it does not have.
3. **Aborting on a bad estimate is worse than not aborting.** A false positive stops a run that
   would have worked. Consider an escape hatch (an input flag or env var to bypass the check)
   and document it in the abort message itself.
4. **The catch handler must not swallow non-OOM errors.** V3 catches `std::exception` broadly;
   make sure the banner distinguishes "GPU out of memory" from "GPU error", since `ASSERT_GPU`
   throws for every failed GPU call, not just allocation.
5. **Exit path.** V3 specifies `std::exit(EXIT_FAILURE)` after reporting. That runs from inside
   an `extern "C"` frame called by Fortran; check it does not deadlock against the measurement
   thread (audit A4 covers that thread's legacy design) or lose buffered output — flush before
   exiting.
6. No AMD hardware available. Anything you touch in the HIP path is unvalidated. Say so.

---

## Suggested commits

```
gpu: report GPU errors and memory state at the Fortran boundary
gpu: expose peak device usage from the tensor memory tracker
gpu: estimate device memory before allocating
gpu: abort with a budget table when a run will not fit
docs: record V3 and M1 in the checklist
```

---

## Definition of done

- V3 lands: an oversized run produces a readable OOM banner with the tracker breakdown and
  free/total, instead of `terminate called after throwing...`. Demonstrate it with an actually
  oversized system and paste the output.
- M1 lands: the estimator agrees with `TensorMemoryTracker` peak device bytes **within 5 % on
  three cases with different flag combinations** (vary `do_jtensor`, `do_dm`, `do_aniso`,
  `do_ene`, convolution on/off, correlations on/off). Report the three measured pairs — actual
  numbers, not "verified".
- An oversized run aborts **before any kernel launch** with the top-consumers table and
  actionable hints.
- `sc_dm_uniaxial` has been run under the new handler and its failure is either explained or
  explicitly declared out of scope with a reason.
- 11 PASS / 1 FAIL preserved on `cuda-debug`, `cuda-fastcopy`, `cuda-ptds`.
- Checklist updated: V3 and M1 ticked **with evidence**, or left unticked and marked partial
  with what is missing. Never tick on impression.
