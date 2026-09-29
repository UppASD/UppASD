# U1 — Fold spectral k to the Brillouin zone · Opus · `[OAM-U1]`

**Why.** C17 folds each reduced index into `[-N/2, N/2)` separately. That region is the unit-cell parallelogram, not the Brillouin zone.
- On square and rectangular cells the two coincide, so those cells are unaffected.
- On the hexagonal cell (`C1=(1,0)`, `C2=(1/2,√3/2)`), the parallelogram edge along x is at k = π. The zone reaches K = 4π/3 ≈ 4.19.
- Modes between the two edges get the wrong derivative. The audit measured, on a 48×48 honeycomb with a Gaussian packet at k = 2.88 along x:

| Fold | λ_centroid |
|---|---|
| Current per-axis fold | −1.50 |
| Brillouin-zone fold | 0.99999997 |
| Exact | 1.000 |

The oracle copied the same rule, so C12a could not catch it. The rest of the spectral code is verified: it matches analytic gradients to 1e-12 with NA = 2, two layers, an oblique cell and 12×9 grids.

**1. Contract.** Replace the "For `spectral`…" paragraph of C17 with this approved text:

> For `spectral`, each sublattice `it` and layer `z` is an `N1×N2` grid. FFT index `(j1, j2)` represents the wave vectors `k(a,b) = ((j1 + a*N1)/N1)*b1 + ((j2 + b*N2)/N2)*b2` for integers `a, b`, where `b1, b2` are `2*pi` times the reciprocal basis of the unit-cell vectors `C1, C2`.
>
> The derivative multiplier is the shortest member of that set: its first-Brillouin-zone (Wigner–Seitz) image.
> - When several members tie for shortest (`|k|^2` equal within a relative `1e-10`), the multiplier is their arithmetic mean. This zeroes the ambiguous component and keeps the unambiguous one.
> - On rectangular cells this reproduces the per-axis fold and Nyquist rule exactly.
> - The search covers `a, b ∈ {-2,…,2}`. `oam_init` refuses if any shortest image lies on the search boundary (`|a| = 2` or `|b| = 2`).
>
> The gradient is exact for fields whose content lies strictly inside the sublattice's first Brillouin zone. Weights, centroid, guards and every output column are unchanged from C3–C6. The `oam_traj` header records `oam_gradient = spectral (Brillouin-zone fold)`.

**2. Fortran** (`oam_setup_spectral`):
- Replace `oam_spectral_mode` and the Nyquist zeroing with the shortest-image search and tie-averaging above.
- Remove `oam_spectral_mode` if nothing else uses it.
- Add the search-boundary refusal.

**3. Oracle (authorised).** Rewrite the k-grid in `spectral_gradient()` independently in numpy. Vectorising over `(a, b)` is fine. Add a selftest assertion that on square and rectangular grids (8×8, 8×7, 9×6) the new k-grid equals the old one exactly.

**4. Harness (authorised).** Add C12c:
- Hex lattice with `high_k_boost` at k = 3.0 along x (vortex width 5, as in C12b).
- Keys: `oam_gradient spectral` and `oam_axis 0 0 1`.
- Pass when Fortran matches the oracle to 1e-6 **and** λ_centroid is within 1e-4 of 1.
- Before changing any code, run C12c against the current binary and put the failing value in the commit message.
- Expected values from the audit oracle: about −4.35 with the current fold, and 0.99999 with the Brillouin-zone fold.

**Acceptance:**
- C0–C13 plus C12c ALL PASS on an FFTW build, and C12 is skipped cleanly without FFTW.
- In Debug builds, C12a (square case) and C12b outputs are byte-identical before and after the change.
- Regression suite 20/20.
