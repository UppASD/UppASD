# U5 — Finalize the OAM refactor · Opus · commits `[OAM-U5a]` … `[OAM-U5g]`

Paste `docs/dev/oam/prompts/00_PREAMBLE.md` above this prompt. The archived round is `docs/dev/oam/prompts/rounds/r4/`.

**Rules for this prompt:**
- Do the steps in order, with one tagged commit per step.
- The maintainer has approved the contract text below; paste it as written.
- Bit-identity comparisons use `CMAKE_BUILD_TYPE=Debug` builds (Release adds `-Ofast`).
- **Final acceptance:**
  - `ctest -L oam` passes on both an FFTW and a non-FFTW build;
  - both harnesses ALL PASS, and `--selftest` passes;
  - regression suite 20/20;
  - `git status` is clean.

---

## U5a — Lock in the band-packet bridge (B5.5)

**Why.** The audit showed that the production kernel reproduces the LSWT Berry term on the angularly varying honeycomb band packet. That result currently lives only in `reference/b55_band_packet_reference.py`. Port it into `run_bridge_checks.py` behind a `--band` flag, and write it independently of `docs/`.

Audit values (D = 1, upper band, ring width 0.15, 90×90 cell, `oam_axis 0 0 1`, one frozen step):

| k₀ | l | spectral | exact gradient | LSWT l − 2⟨F⟩ | FEM |
|---|---|---|---|---|---|
| 1.0 | 0 | 0.000606 | 0.000606 | 0.000689 | 0.000525 |
| 1.0 | 1 | 1.000602 | 1.000606 | 1.000689 | 0.882 |
| 2.0 | 0 | 0.040822 | 0.040822 | 0.042170 | 0.0232 |
| 2.0 | 1 | 1.040821 | 1.040822 | 1.042170 | 0.602 |

**Acceptance** (spectral only; FEM is reported, not gated):
- \|spectral − exact\| ≤ 1e-5;
- the l=1 minus l=0 shift is within 1e-5 of 1;
- \|exact − LSWT\| ≤ 3e-3 (the gap is the ⟨F⟩ quadrature, not physics).

Register it as CTest `oam-bridge-band`, labels `oam-traj;oam-bridge`. Skip it with a notice when the binary lacks FFTW. Runtime is about 30 s.

## U5b — Default `oam_gradient auto` (amends C17)

Replace C17's first paragraph with:

> **C17 — Gradient method (amended 2026-09-29).** The optional key `oam_gradient auto|fem|spectral` selects the trajectory gradient method.
> - The default is `auto`. It resolves to `spectral` when the build defines `USE_FFTW`, `BC1 = BC2 = 'P'`, and the shortest-image search stays inside its bounds. Otherwise it resolves to `fem`.
> - The resolution is printed once to stdout and recorded in the header as `# oam_gradient = <method> (auto)`.
> - An explicit `spectral` keeps the strict refusal described below. An explicit `fem` always uses FEM.

**Tasks.**
- Implement the resolution in `oam_init`.
- The harness must pass `oam_gradient` explicitly in every existing case, so each keeps its current method and oracle. Existing C0–C13 outputs must stay byte-identical (Debug).
- Add C15, with three sub-cases:
  - `auto` on a periodic FFTW build → header says `spectral (auto)`, and the result matches the spectral oracle;
  - `auto` with open BC → `fem (auto)`;
  - `auto` on a non-FFTW build → `fem (auto)`, with no refusal.
- Update the header text, `OAM_IO.md` and the U2 validity sentence: "FEM is biased at short wavelength; `auto` selects spectral where possible."

## U5c — Ensemble aggregation (new clause C18)

Add to the contract:

> **C18 — Ensembles (approved 2026-09-29).** With `Mensemble > 1`, each sample aggregates over the ensembles `k`.
> - **λ columns** (origin, centroid and per-sublattice) are norm-weighted: `λ = Σ_k L_k / Σ_k n_k`, where `L_k = Σ_i ℓ_i w_i` and `n_k = Σ_i |ψ_i|² w_i`. The sums run over the ensembles for which that column is valid: the norm guard for origin columns; the norm and spread guards for centroid columns.
> - **`N_m` and `dSz_hbar`** are arithmetic means over all ensembles.
> - **`Lz_tot_hbar` and `balance`:** `Lz_tot_hbar = N_m × lambda_L_centroid` and `balance = dSz_hbar + Lz_tot_hbar`, both from the aggregated values, so the C5 relations hold.
> - **`R_x`, `R_y`, `sigma_psi`** are `n_k`-weighted means over the centroid-valid ensembles.
> - If any ensemble is excluded from a λ column at a sample, stdout gets a one-time warning with the count.
> - With `Mensemble = 1`, output is byte-identical to before.

**Tasks.**
- Restructure `oam_evaluate_ensemble` and `oam_sample` to return and aggregate `L_k` and `n_k`.
- Add the aggregation to the oracle as an authorised `evaluate_ensembles()`.
- Add harness check C14. Restart file with `Mensemble = 2`:
  - ensemble 1: the l=+1 vortex;
  - ensemble 2: an l=+2 vortex of different amplitude;
  - pass when the output matches the oracle to 1e-10.
- Add a second C14 case where ensemble 2 is delocalised (the C9 field). Pass when the centroid columns come from ensemble 1 only and the warning is printed.

## U5d — GPU path copies moments on OAM steps (bug)

`CudaMeasurement::measure` copies moments from the device only when `fortran_do_measurements(mstep)` is true. `do_measurements` in `measurements.f90` does not know about `do_oam_traj`. On an OAM sample step with no other copy-triggering measurement, `oam_sample` therefore runs on **stale** `FortranData::emom`, without any warning.

**Tasks.**
- Export a function from the OAM module, `logical function oam_sample_due(mstep)`, using the same `mod(mstep-oam_rstep-1, oam_step_traj)` rule as `oam_sample`.
- In `do_measurements`, set `do_copy = 1` when `do_oam_traj == 'Y' .and. oam_sample_due(mstep)`. Note that the OAM rule uses `mstep`, not the log-sampled `sstep`.
- Check the `gpu_mode = 2` (C++) path and the other drivers (`ms_driver`, `sld_driver`, `sx_driver`) for the same pattern, and fix any you find.
- If a CUDA toolchain is available, compare a 100-step GPU run against the CPU run on the C12b fixture, with `oam_step 10` and every other measurement off. Otherwise, state in the commit message that this was verified by inspection only.

## U5e — Strictness and small fixes

- An unknown `oam_weight` currently falls back silently to `site`. Make it refuse, like an unknown `oam_gradient`.
- `OAM_IO.md`, at the image-search refusal: a refusal means the in-plane cell vectors are not reduced. The user should re-express `C1, C2` as a reduced basis (for example, replace `C2 = (1.5, 0.3)` by `C2 − C1 = (0.5, 0.3)`).
- Measure and report in `OAM_IO.md` (a short table; no gate):
  - FFTW planning time, and the per-sample cost of spectral against FEM;
  - sizes 128², 256², 512² at NA = 2, `oam_step 10`;
  - OAM time as a fraction of total run time.
- If planning takes more than 1 s at 512², switch `FFTW_MEASURE` to `FFTW_ESTIMATE`. The harness must still pass.

## U5f — Physical-run smoke tests (report only)

1. **Damped packet:** the clean l=1 packet from U3, with damping 0.01, T = 0, 2000 steps, `oam_step 10`, `auto`. Report `N_m` and `lambda_L_centroid` against step (a table every 200 steps).
2. **Thermal run:** FM ground state at T = 10 K, 2000 steps, `Mensemble = 4`. Expected behaviour, to state in `OAM_IO.md`:
   - it doesn't crash;
   - origin columns are finite;
   - `lambda_L_centroid` is mostly NaN, because the spread guard correctly rejects delocalised thermal magnons.
   - λ is a coherent-packet observable.

Put both tables in `docs/dev/oam/STATUS.md` (U5g).

## U5g — Docs and cleanup

1. **`BRIDGE.md`:** line 4 (status), the §3 closing paragraph and the end of §6 must now say that the production kernel is validated against LSWT with spectral gradients (B5.5), and is biased with FEM.
2. **`OAM_QUESTIONS.md`:** add a one-line pointer to `BRIDGE.md` §6 under the dated B5.4 entry that quotes "targets F_n and +1".
3. **New `docs/dev/oam/STATUS.md`**, one page:
   - what is validated, and by which checks;
   - the documented scope: collinear FM reference; periodic in-plane cells for spectral; content inside the first Brillouin zone; coherent packets; MKL-FFT builds resolve `auto` to `fem`;
   - known open items: the LSWT magnitude gap O₁,av(K) = 0.228 against 0.236 published (−3.4%, sign confirmed); the carrier-referenced gradient proposal;
   - the U5f tables.
4. **Move the prompts:** archive the four rounds under `docs/dev/oam/prompts/rounds/{r1,r2,r3,r4}/`, and fix any links.
5. **Remove stray outputs:** `traj_checks/`, `*.megaTest.*`, and build or work directories. Add ignore rules for the harness work directories.
