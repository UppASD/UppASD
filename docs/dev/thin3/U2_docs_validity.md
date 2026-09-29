# U2 — Finish the T2 docs; quote λ bias · Sonnet · `[OAM-U2]`

**No source changes except output-header text.** Run after U1.

**Still missing from T2:**
1. `BRIDGE.md` §3 (around line 165) and `OAM_QUESTIONS.md` (around line 35) still say the per-sublattice mesh "cannot recover" F. Rewrite both with the T2 explanation:
   - `Im[u†∂φu] = Σ_s Im[u_s* ∂φ u_s]` is a sum of per-sublattice terms, so per-sublattice sampling is enough.
   - The limit is gradient accuracy at short wavelength.
2. `BRIDGE.md` §5: state plainly that B5.4 validates the oracle against LSWT, not the production trajectory kernel.
3. `OAM_IO.md`:
   - Document `oam_gradient fem|spectral` (C17, as amended by U1) and `oam_axis x y z` (C8b), with defaults and refusal conditions.
   - Spectral needs a build that defines `USE_FFTW`. MKL-FFT builds don't define it, so they refuse spectral.
4. `OAM_QUESTIONS.md` (around line 83): mark the `oam_axis` proposal as implemented (C8b).

**Fix the validity numbers.** The current text quotes the bias of the phase gradient (−4.6% and −16.3%). Users read λ, whose bias is much larger. The audit's λ_centroid values for the C12b vortex (width 5) boosted along x, FEM gradient, `oam_axis 0 0 1`:

| k·a | 0 | 0.25 | 0.5 | 1.0 | 1.5 |
|---|---|---|---|---|---|
| Square | 0.987 | 0.966 | 0.906 | 0.685 | 0.376 |
| Hex | 0.990 | 0.975 | 0.929 | 0.758 | 0.506 |

Reproduce these values yourself with `oracle_traj.evaluate(..., gradient="fem")` and quote your own numbers, not these.

Then replace the validity sentence in C3, `OAM_IO.md` and the `oam_traj` header with:

> FEM λ is biased low as k·a grows (square lattice: about −X% at k·a = 0.5, −Y% at 1.0, relative to k·a = 0). Use `oam_gradient spectral` on periodic cells for k·a ≳ 0.25.

**Header.**
- Print the FEM validity line only when `oam_gradient = fem`.
- In spectral mode, print one line instead: `# spectral gradient: exact inside the first Brillouin zone; see C17`.

**Frame note (`OAM_IO.md`).** For boosted, driven or restart-loaded states, recommend setting `oam_axis` to the ground-state axis. In the C13 case, the default axis tilts by about 0.2° and shifts λ_centroid by 6% (0.903 against 0.957).

**Acceptance:**
- `git grep -n "cannot recover" -- docs tests source ':!docs/dev/thin*'` returns nothing.
- Both harnesses still ALL PASS.
