# T2 — Correct BRIDGE §3 and state the validity range · Sonnet · `[OAM-T2]`

**No source changes except output-header text.**

**The correction.** `BRIDGE.md` §3 says a per-sublattice mesh cannot recover F because no triangle spans A and B. That is wrong. The intrinsic term is `Im[u†∂φu] = Σ_s Im[u_s* ∂φ u_s]`, a sum of per-sublattice terms, so per-sublattice sampling is sufficient. The limit is the accuracy of the linear-FEM gradient at short wavelength.

The audit measured, on honeycomb ring packets with D = 1 and l = 1:

| k₀ (1/a) | Expected `l − 2⟨F⟩` | Exact gradient | Linear FEM |
|---|---|---|---|
| 1.0 | 1.0007 | 1.0006 | 0.882 |
| 1.5 | 1.0073 | 1.0071 | 0.751 |
| 2.0 | 1.0415 | 1.0408 | 0.602 |

The B5.4 production run at k₀ = 2.9 (shift −0.22) is the same effect.

**Tasks.**
1. Rewrite `BRIDGE.md` §3 and the matching sentences in `OAM_QUESTIONS.md` (the B5.4 entries) and `run_bridge_checks.py`'s docstring, using the explanation above. Keep the fixed-spinor-control statement; it is correct.
2. State plainly in `BRIDGE.md` §5 that B5.4 validates the oracle against LSWT, not the production trajectory kernel.
3. Add a validity statement to `OAM_IO.md`, `CONVENTIONS_OAM.md` (end of C3) and the `oam_traj` header: "Linear-FEM gradients are accurate for k·a ≲ 0.5; the error grows with k·a (about −12% at k·a = 1)."

   Before writing the number, measure it yourself with `oracle_traj.py` on a single-sublattice square lattice: a boosted Gaussian with k·a = 0.5 and 1.0. Quote your numbers, not these.

**Acceptance:** docs are consistent (`git grep -n "cannot recover"` is empty); both harnesses still ALL PASS.
