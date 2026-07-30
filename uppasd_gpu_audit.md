# UppASD GPU port (`gpu_hip_cu_ab`) — code audit

Scope: `source/gpu_files/**`, plus `source/chelper.f90` and `source/uppasd.f90` (and `sd_driver.f90` where needed to follow the call chain). Reviewed at commit `1e0ff85` ("Step 2 in codex review"). Findings are ordered by level: architecture, correctness, performance, interop, build/CI, hygiene, then a suggested roadmap. Line numbers are approximate but close.

---

## 1. Executive summary

The overall porting strategy is sound and worth keeping: persistent device state across the whole phase, Fortran demoted to setup/I/O, a functor-based kernel dispatch (`GpuParallelizationHelper` + `each()` operators), a single-source CUDA/HIP abstraction via `gpu_wrappers.h`, and an FFT convolution backend for Bravais lattices. The Fejes-era architecture has aged well and the newer additions (reduced Hamiltonian support, GPU-side measurements with two-phase reductions, the convolution backend) fit into it naturally.

The three things that most need attention, in order:

1. **The block-wide energy reduction in the Heisge kernels is broken in three independent ways** (divergent `__syncthreads`, stale shared memory, missing barrier between successive reductions). Energy measurement is unreliable for any `N` not a multiple of 256 and racy even when it is. (C1)
2. **The thermal-field RNG call can silently fail for odd `3·N·M`** because cuRAND/hipRAND normal generators require even counts and the status is never checked. This corrupts the thermodynamics silently. (C2)
3. **Maintainability debt is compounding**: six hand-copied ~150-line Heisge kernels, two ~1100-line correlation kernel files that are 90 % identical between CUDA and HIP, `old/` shipping dead code, heavy commented-out blocks, and no GPU compilation in CI (which is how a Fortran syntax error in `chelper.f90` is currently sitting on the branch). (A2, A3, B1, CI)

None of this argues for redoing the port differently. It argues for a consolidation pass before more physics is added.

---

## 2. Architecture-level assessment

### A1. Overall design — keep it
The execution model (upload once per phase, run the entire phase on device, download at phase boundaries, `initiate → run → release` per phase driven from `sd_driver.f90`/`chelper.f90`) is the right one for spin dynamics. The functor pattern (`class Foo : ParallelizationHelper::AtomSiteEnsemble { __device__ void each(...) }`) keeps kernels readable and gives you index-mapping in one place. The transposed (site-fastest) neighbor layout for `ncoup`/`nlist` is the correct coalescing choice.

### A2. CUDA/HIP unification is incomplete — finish it
You already have the right machinery (`gpu_wrappers.h` with `GPU_*` macros, `WARPSIZE`, `SHFL_DOWN`, `GPU_CREAL`, `MAKE_GPU_COMPLEX`), but two subsystems still carry duplicated sources:

- `cudaThermfield.cu` vs `hipThermfield.cpp`: near-identical (a ~44-line diff dominated by API names already covered by your macros).
- `correlations/correlation_kernels.cpp` vs `correlation_kernels.hip.cpp`: ~1100 lines each, ~90 % identical. The real deltas are (i) `cg::tiled_partition<32>` vs `<64>`, (ii) `thrust::complex` member access vs `GPU_CREAL`/`GPU_CIMAG`, (iii) cooperative-groups index idioms vs raw `blockIdx`.

All three deltas are expressible in the existing macro layer plus one compile-time constant, e.g.

```cpp
#if defined(__HIP_PLATFORM_AMD__)
  #define WAVE_SIZE __AMDGCN_WAVEFRONT_SIZE   // 64 on CDNA, 32 on RDNA
#else
  #define WAVE_SIZE 32
#endif
```

Note the current `tiled_partition<64>` hardcode is itself a portability bug: it is wrong on RDNA (wave32) consumer cards, which people do use for development. Unify to one source per subsystem, compiled as CUDA or HIP by the build system, and delete the `.hip.cpp` twins. This is the single largest maintenance win available.

### A3. Kernel combinatorics in the Hamiltonian — template instead of copy
`gpuHamiltonianCalculations.cpp` contains `HeisgeJij`, `HeisgeJijDM`, `HeisgeJijAniso`, `HeisgeJijDMAniso`, `HeisgeJijTensor`, `HeisgeJijTensorAniso` — six hand-maintained variants of the same kernel, and the anisotropy/energy code inside them has already drifted (some have `__restrict__`, some don't; the energy prefactor comments differ). Every new interaction (biquadratic, ring exchange, SA...) doubles this. Replace with one kernel:

```cpp
template<bool HasDM, bool HasAniso, bool HasTensor, bool Measure>
class Heisge : public ParallelizationHelper::AtomSiteEnsemble {
   __device__ void each(...) {
      ... exchange ...
      if constexpr (HasDM)    { ... }
      if constexpr (HasAniso) { ... }
      if constexpr (Measure)  { ... energy path ... }
   }
};
```

with a small dispatch table in `heisge()`. Dead branches are compiled out, so performance is identical to the hand-copied versions, and making `Measure` a template parameter also removes the energy path's register pressure from the 99 % of steps that don't measure (see P3). This is the highest-leverage refactor in the file and directly mirrors what you did in the `librsrec` batching work.

### A4. The measurement thread architecture is legacy — decide its future
The pthread `MeasurementQueue` + Fortran-callback design predates GPU-side measurements. Now that `GpuMeasurement` computes averages/cumulants/skyrmion number/energy on device, the thread mainly services the CPU fallback. Its concurrency model is fragile: the worker thread calls Fortran module procedures (`measure`, `flush_measurements`) which hold global state, while the main thread also calls into Fortran (`fortran_calc_simulation_status_variables` from `printMdStatus`, `copyToFortran` writing `emomM` that the empty-measurement path reads live). Fortran module state is not thread-safe in general. Either (a) complete GPU measurement coverage and retire the thread, or (b) constrain it: snapshot-only (never pass live `FortranData::` pointers to the worker) and never call Fortran from two threads in the same window.

### A5. `old/` and dead code
`gpu_files/old/` ships an entire previous generation (`cudaMatrix.hpp`, `mdSimulation.cpp`, ...). The active files also carry large commented-out regions ("TEMPORARY PRINTING", disabled `cpuRestMeasurement`, dead `null_energy` kernel, the disabled `fortran_print_measurables` twin at the bottom of `chelper.f90`). Git remembers; delete. It materially slows down anyone (human or model) auditing the tree.

---

## 3. Correctness findings

### C1 — `sum_warp_energy`: three bugs in one function (severity: high)
`gpuHamiltonianCalculations.cpp:15-49`, used by every Heisge kernel when `measure == true`.

**(a) Divergent `__syncthreads()`.** The kernels are launched via `atom_site_ensemble_kernel`, which calls `op.each(...)` only when `GridHelper::index2d` returns true, i.e. only for threads with `site < N`. In the tail block (whenever `N % 256 != 0`), out-of-range threads skip `each()` entirely and therefore never reach the `__syncthreads()` inside `sum_warp_energy`. That is a divergent barrier: undefined behavior, and on Volta+/CDNA with independent thread scheduling it can hang or silently desynchronize.

**(b) Stale shared memory.** `__shared__ real warp_sums[32]` is written only by lane 0 of warps that have in-range threads, but warp 0 reads `warp_sums[lane]` for all `lane < num_warps` where `num_warps` is computed from `blockDim` — including warps that were entirely out of range and never wrote their slot this launch. Garbage from previous launches gets summed into the energy.

**(c) Missing trailing barrier.** The kernels call `sum_warp_energy` three times back-to-back on the same shared array. Between the read of `warp_sums` by warp 0 in call *k* and the write of `warp_sums[w]` by warp *w* in call *k+1* there is no `__syncthreads()` — a straightforward shared-memory race.

**Fix.** The cleanest repair is to move the bounds check *inside* `each()` for the measuring variants so all threads enter the reduction (out-of-range threads contribute 0), zero-initialize `warp_sums` (or guard the read with the active-warp count), and add a `__syncthreads()` at the end of the function. Better still: you already have a correctly written two-phase block-reduction library in `measurement/kernels.cpp` (`block_reduce_sum_1d` + partial/finalize kernels). Reuse that machinery for the Hamiltonian energies instead of maintaining a second, buggy reduction.

### C2 — RNG length not checked; odd counts silently fail (severity: high)
`cudaThermfield.cu:144-147` / `hipThermfield.cpp:146-149`:

```cpp
curandGenerateNormalDouble(gen, field.data(), field.size(), 0.0, 1.0);
```

cuRAND (and hipRAND) require the count to be **even** for the normal/log-normal generators with XORWOW/MRG32k3a etc.; otherwise they return `CURAND_STATUS_LENGTH_NOT_MULTIPLE` and generate nothing. `field.size() = 3·N·M` is odd whenever both N and M are odd. The return status is discarded, so the thermal field is left stale/uninitialized — silently wrong Langevin dynamics for e.g. N=1001, M=1. Fix: pad the allocation to even length (generate `size+size%2` into a buffer with one slack element) and `assert`/check every `curand*`/`hiprand*` status. The same applies to `hiprandGenerateNormal`.

### C3 — `chelper.f90:20` — syntax error (severity: build-breaking under gfortran)
```fortran
   use prn_averages,  
   use prn_topology,     only : skyno, ...
```
A `USE` statement with a trailing comma and no rename/only list is invalid Fortran; gfortran rejects it. Since CI never compiles the GPU configuration (see §6), this sits undetected — it looks like an artifact of a recent edit (possibly the "codex review" pass). Restore the intended `only :` list.

### C4 — In-place transpose of Fortran-owned arrays (severity: medium, latent)
`gpuSimulation.cpp:494-500`: `copyFromFortran()` transposes `cpuHamiltonian.ncoup`/`nlist` **in place inside the arrays owned by the Fortran `HamiltonianData` module**, copies to device, then transposes back. This (i) temporarily corrupts the Fortran view — any concurrent reader (OpenMP regions, the measurement thread, a signal-triggered dump) sees garbage; (ii) leaves the arrays permanently transposed if anything exits between the two transposes; (iii) does two full O(mnn·NH) host-side transposes per phase for no reason. Since the device copy is persistent, transpose into a scratch buffer once (or upload raw and transpose on device with a trivial kernel). This dovetails with the SINGLE_PREC issue in I2 — a single "upload with layout/type transform" function fixes both.

### C5 — MC acceptance is wrong for anisotropy (severity: high if GPU-MC + aniso is used)
`gpuMetropolis.cpp:57-95` (`MCSweep`): the acceptance uses `ΔE ∝ eneff·(S_new − S_old)` with `eneff` frozen from the last `heisge`. For strictly bilinear terms (Heisenberg, DM, Zeeman) under a valid sublattice decomposition this is exact — spins in one sublattice don't couple, so the field at each site is constant during the sweep. But uniaxial/cubic anisotropy energies are quadratic/quartic in the *on-site* spin: the "energy field" `eneff` was built from the **old** spin, so the linearized ΔE is simply not the anisotropy energy difference. The Boltzmann weights are wrong whenever `do_aniso != 0` in the GPU-MC path. The on-site anisotropy ΔE is cheap to evaluate exactly inside `MCSweep` (both spins are in registers); do that, or hard-refuse the combination until it is done.

### C6 — Sublattice independence derived only from the exchange list (severity: medium)
`split_lattice()` builds the coloring from `nlist` alone. If `dmlist` ever contains a pair not present in `nlist` (different cutoffs are legal input), two spins in the same sublattice interact via DM and the parallel update violates detailed balance. Color on the union of all active neighbor lists.

### C7 — `USE_BIG_GRID` breaks the ensemble index in energy/MC accumulation (severity: high when enabled)
The Heisge energy path and `MCSweep` read `mInd = blockIdx.y` directly. In the default small-grid mapping `blockIdx.y` *is* the ensemble, but under `USE_BIG_GRID` (the documented escape hatch for very large N) `GridHelper::index2d` derives the ensemble as `xy / X`, and `blockIdx.y` becomes a grid-splitting index. Energies get accumulated into wrong (or out-of-range) `energyM` rows exactly in the large-system regime the flag exists for. Pass the ensemble that `each(atom, site, ensemble)` already receives instead of touching `blockIdx` inside functors — that's the abstraction's whole point.

### C8 — `fastCopy` measurement path has a host/DMA race (severity: medium; path off by default)
`fortranMeasurement.cpp:99-135` (`copyQueueFast`): the queue callback is enqueued on `workStream`, which was made to wait only on `copyDone` — an event recorded **after the D2D staging copies but before the D2H copies** to the pinned buffers on `copyStream` (intentionally so, to unblock compute early — that part is good design). But the callback then `memcpy`s from the pinned buffers with no ordering against the in-flight D2H transfers. Record a second event after the D2H copies and enqueue the callback behind it (or put the callback on `copyStream`). Worth fixing now: `USE_FAST_COPY` is exactly the switch someone will flip when measurement overhead shows up in profiles, and the failure will look like intermittently corrupted averages.

### C9 — CUDA runtime touched in static initialization (severity: low-medium)
`ParallelizationHelperInstance` is a global; its `GridHelper` member constructor calls `cudaGetDevice`/`cudaGetDeviceProperties` before `main()`, unchecked. On a node without a visible device (login nodes, `CUDA_VISIBLE_DEVICES` games, Slurm prolog quirks) `maxGridSize1` is left uninitialized and every later `dim1d` decision is garbage. Query the device inside `initiate()` and check the error.

### C10 — Smaller correctness notes
- `GpuDepondtIntegrator::rotate(emom, delta_t)` ignores `delta_t` and uses the member `timestep`. Harmless today, a trap tomorrow.
- `Rotate::each` divides by `norm` of the effective field without a zero guard; a zero local field (possible at T=0 with no external field in symmetric configurations) produces NaNs that propagate silently.
- `SetupNeighbourList` mutates device `nlist`/`ncoup` in place (Fortran→C index shift + padding). This is safe only under the invariant "exactly one `initiate()` per fresh upload", which currently holds because each phase re-uploads — but nothing enforces it. A second `initiate()` without a re-upload silently corrupts the neighbor lists. Cheap guard: a `listsPrepared` flag, or make the setup idempotent by writing to a device-owned copy.
- Energy measurement is taken on the *predictor* configuration (heisge after `evolveFirst`), so measured energies are O(dt) off the corrected trajectory. Probably fine; verify it matches the Fortran convention or note it.
- `heisge(..., measure=true)` bypasses the convolution backend and falls back to the sparse path. Correct, but it means energy-measuring steps exercise different code/rounding; also the sparse structures must always stay resident even in convolution mode (they do today). Document this coupling, and consider computing energies from the convolved field so the FFT path covers measuring steps too.

---

## 4. Performance findings

### P1 — MC path: full-system field recomputation per sublattice (largest MC win)
`gpuMCSimulation.cpp:131-138, 240-245`: each MC step does `for sub in num_subL { MCSweep(sub); heisge(FULL SYSTEM); }`. One Metropolis step costs `num_subL` complete Hamiltonian evaluations. The standard alternative removes field precomputation from the sweep entirely: compute ΔE locally inside `MCSweep` by gathering the trial spin's neighbors from `nlist`/`ncoup` (all resident on device). Cost per attempt becomes O(mnn) reads — the same as the field kernel's per-site cost — and the per-sublattice `heisge` disappears. That is roughly a `num_subL`× speedup of the MC path, and it simultaneously fixes C5 because you can then evaluate the exact on-site anisotropy difference. Recompute `beff` once per step only when a measurement needs it.

### P2 — Stream discipline and gratuitous device-wide syncs
`GPU_DEVICE_SYNCHRONIZE()` appears at the end of every `GpuMeasurement::measure()`, after every correlation sample (`gpuCorrelations.cpp`, several sites), and at phase ends. Meanwhile the correlation kernels launch on the **default stream** while the integrator/Hamiltonian run on `workStream` — today this is only correct because legacy default-stream semantics serialize everything, which is also why it's slow. Pass `workStream` to every launch (the correlation kernels take no stream argument at all right now), keep the copy/work event choreography you already have in `copyQueueFast`, and synchronize only where the host actually consumes results (`flush`, `printMdStatus`). This would also make the code safe under `--default-stream per-thread`, which someone will eventually enable.

### P3 — Energy path bloats every field kernel
Because `measure` is a runtime bool, the energy code (extra registers, shared array, atomics) is compiled into the hot field kernel used every step. Registers are allocated for the worst-case path, reducing occupancy even when `measure == false`. Making `Measure` a template parameter (A3) gives the lean kernel for production steps for free.

### P4 — `j_tensor` layout is fully uncoalesced
The tensor path indexes `tensor[i + 3j + 9k + 9·mnn·l]` (Fortran layout, site `l` slowest). Consecutive threads (consecutive sites) hit stride `9·mnn` — every load is its own cache line. You already transpose `ncoup`/`nlist` to site-fastest for the scalar path; do the same for the 3×3 tensors (layout `[site + NH·(component + 9·neigh)]`). Expect a large bandwidth win in tensor mode.

### P5 — `split_lattice` is O(N²·mnn) with debug prints
`gpuMetropolis.cpp:155-300`: the greedy coloring scans all remaining sites against a growing neighbor union, and then prints the entire `num_subL × max_spins` table to stdout (with `setbuf(stdout, nullptr)` in effect, so every integer is an unbuffered write). For large N this dominates startup. Greedy graph coloring on the adjacency list is O(N·mnn): for each site, mark colors of already-colored neighbors, take the smallest free color, then bucket sites by color. Delete the prints.

### P6 — Integrator kernel fusion
`evolveFirst` is five element-wise passes (`Add`, `BuildEffectiveField`, `AddTo` for STT, `Rotate`, `ScalarMult`), each streaming the full 3NM arrays. All are per-atom with no cross-atom dependencies — fuse local-field + damping-term + rotation (+ `emomM` update) into one kernel and halve the integrator's memory traffic. This starts to matter exactly when the convolution backend makes the Hamiltonian cheap.

### P7 — Precision details
- `GPU_NORMAL_DOUBLE`/`GPU_RAND_UNIFORM_DOUBLE` always generate doubles, even under `SINGLE_PREC`; use the float variants conditionally.
- `Tensor::operator()` contains `assert(valid_index(...))` executed on device; ensure production builds define `NDEBUG` (note: the `-DDEBUG` that CMake injects is unrelated to `assert` — see B2). Device asserts in the innermost index computation of every kernel are expensive.

### P8 — Convolution backend
The FFT path is well-structured (kernel built once, two C2C transforms per `heisge`). Two cheap follow-ups: (i) the spin field is real — R2C/C2R halves FFT work and spectral-multiply footprint; (ii) `heisge(measure=true)` currently abandons it (see C10) — computing energies from the convolved field would let long measured runs stay on the fast path.

---

## 5. Interop layer (Fortran ↔ C++)

### I1 — Mixed binding conventions
`FortranData_setCorrelations` and `FortranData_setGpuGeometry` have proper `bind(C)` interfaces in `chelper.f90`; the other six `FortranData_set*` calls (plus `gpusim_*`) rely on implicit externals matching hand-mangled `extern "C" void fortrandata_setconstants_(...)` symbols — the classic trailing-underscore convention, plus the assumption that hidden character-length arguments are appended at the end and ignored. This works with gfortran and current flags, but it is compiler-convention-dependent (Cray/NVFortran/`-fno-underscoring` all break it) and it bypasses any argument checking. Since you already started `bind(C)` interfaces, finish the job — it's mechanical and removes a whole class of silent ABI breakage. While at it: Fortran default `integer` is passed where the C side declares `unsigned int*`; harmless for the values in play, but the `bind(C)` rewrite should settle on `integer(c_int)`/`int` consistently.

### I2 — `SINGLE_PREC` is incompatible with the bridge
The Fortran side always passes `real(dblprec)` (double) arrays; the C side types everything `real*`, which is `float*` under `SINGLE_PREC`. A single-precision GPU build therefore reinterprets doubles as floats over the whole bridge. Either fail at configure time when `SINGLE_PREC` is set, or (better, and the real goal if fp32 throughput on consumer/RDNA cards is wanted) convert at the upload boundary: host arrays stay double, `copyFromFortran` converts into fp32 device buffers. Combined with C4, this argues for one dedicated staging/transform layer for all uploads.

### I3 — `fortran_print_measurables` passes `allocatable, intent(in)` dummies across the callback boundary
Called from C++ via the measurement writer with raw pointers on the actual-argument side; assumed-shape/allocatable dummies require the caller to pass a descriptor, which C++ cannot construct without `bind(C)` + `CFI_cdesc_t`. If this currently works it is because the call chain happens to go Fortran→Fortran; audit that path once when doing I1.

---

## 6. Build system and CI

### B1 — No GPU configuration is ever compiled in CI
`.github/workflows/{linux,mac,win}builds.yml` build CPU-only. That is how C3 (a syntax error in a GPU-only file) can live on the branch. You don't need GPU hardware to catch 95 % of regressions: add a compile-only job with the CUDA toolkit (nvcc compiles everything without a device) and one with ROCm/hipcc in a container (`rocm/dev-ubuntu` images work on plain runners). Given the regression infrastructure you built for the RS-LMTO `fable` branches, the same pattern applies here; even a `cmake --build` smoke job would have caught C3, and a T=0 deterministic comparison job (thermal field off, fixed seed, compare GPU vs Fortran trajectories for ~100 steps across the J / J+DM / J+aniso / tensor / convolution permutations) would catch layout/indexing regressions essentially for free on a self-hosted GPU runner when one is available.

### B2 — CMake details
- The CUDA arch workaround `if (CMAKE_CUDA_ARCHITECTURES MATCHES "^[__86+PTX]")` is a regex character class: it matches any arch string *starting with* one of `_`, `8`, `6`, `+`, `P`, `T`, `X` — so e.g. a detected `80` or `61` silently gets pinned to 86. Use `CMAKE_CUDA_ARCHITECTURES=native` (CMake ≥ 3.24) and delete the workaround.
- `-DDEBUG` is appended unconditionally to `CMAKE_CUDA_FLAGS`/`CMAKE_HIP_FLAGS`; nothing in the tree uses `DEBUG`, and it does not control `assert` (that's `NDEBUG` via build type). Decide what it's for or remove it, and make sure release builds actually define `NDEBUG` (P7).
- Setting `CMAKE_CXX_COMPILER hipcc` *inside* `CMakeLists.txt` after `project()` is unreliable (CMake caches the compiler at first configure); prefer a toolchain file or `CC=/CXX=` environment, which is also what LUMI's PrgEnv expects.
- `gpuParallelizationHelper.tpp` and headers are listed in `target_sources` (harmless), while `gpuMeasurement.cpp`/`measurementQueue.cpp` are commented out in the subdir list but compiled via the top-level `set_source_files_properties` list — the source lists have drifted; reconcile them.
- The `ON_LUMI` thrust-vs-hipComplex fork in `real_type.h` becomes unnecessary once A2 lands (pick one complex representation for HIP everywhere).

---

## 7. Hygiene (quick sweep)

`#pragma once` in translation units (`fortranData.cpp`, `gpuMetropolis.cpp`, `gpuCorrelations.cpp`, `gpuParallelizationHelper.cpp`); mixed CRLF/LF across the tree; `std::setbuf(stdout, nullptr)` in production phase functions (unbuffered stdout is slow and only needed when debugging crashes); the stray `real cv{};  // Specific heat` pasted inside the SDmphase step loop (`gpuSDSimulation.cpp:~305`); dead locals `mnn`, `l`, `NH` recomputed from `j_tensor` extents in both SD phases; typo'd identifiers worth fixing while APIs are young (`redNeibourCount`, "Depontd"); `energyM` column meanings (0 exch, 1 aniso, 2 dm, 3 tensor, 4 ext, 5 total) as magic numbers — make an enum; the "To Anders: sign????" comments in `HeisgeJijAniso` mark genuinely unresolved physics questions in the energy prefactors (uniaxial `ax_en` sign, cubic /2 and /3 factors) — these should be settled against the Fortran reference and the comments deleted, since they currently flag the energy output as untrusted; `GridHelper::dim1d`'s "are we exiting or returning?" — return false and propagate rather than `exit()` from a library layer; the misleading "multiple of 32" performance note on HIP wave64; `pthread_attr_setschedpolicy` without `PTHREAD_EXPLICIT_SCHED` is a no-op (and `std::thread`/`std::mutex` would be more portable anyway); `cudaStreamAddCallback` is deprecated — `cudaLaunchHostFunc`/`hipLaunchHostFunc` is the modern equivalent and has cleaner semantics for C8.

---

## 8. Suggested roadmap

Ordered so each step de-risks the next; (1)–(3) are small and unblock trust in the code, (4)–(6) are the structural investments.

1. **Correctness triage** — fix C1 (reuse the `measurement/kernels.cpp` reduction), C2 (even-length padding + status checks), C3 (the `use` statement), C7 (take `ensemble` from `each()`), and guard C10's `initiate` invariant. Add the `Rotate` zero-field guard.
2. **CI** — compile-only CUDA + HIP jobs; then a T=0 deterministic GPU-vs-Fortran regression across Hamiltonian permutations.
3. **MC path** — local ΔE inside `MCSweep` (fixes C5, delivers P1); union-of-lists coloring (C6); O(N·mnn) coloring without prints (P5).
4. **Single-source unification (A2)** — merge thermfield and correlation kernel twins onto `gpu_wrappers.h` + `WAVE_SIZE`; drop `ON_LUMI`.
5. **Heisge templating (A3/P3)** — one kernel, compile-time feature flags including `Measure`; transpose `j_tensor` (P4) while touching the layouts; introduce the staging/transform upload layer (C4 + I2) at the same time.
6. **Stream discipline (P2)** — stream parameter through all launches, remove blanket device syncs, fix the `fastCopy` event ordering (C8) and consider making it the default; then integrator fusion (P6) and convolution follow-ups (P8) as profiling dictates.

The port does not need to be redone — the architecture is the right one. It needs the energy/RNG bugs fixed, a CI net under it, and one determined consolidation pass so that the next round of physics (and the MC rework) lands on a single, templated, stream-correct code path instead of six diverging copies.
