# U3 — Rerun and report the bridge with the correct targets · Sonnet · `[OAM-U3]`

**Needs U1.** Report only; don't gate.

**Correction to T3a.** The production packet in `run_bridge_checks.py` uses a **fixed** sublattice spinor. Nothing varies with angle, so its targets are λ(l=0) ≈ 0, λ(l=1) ≈ 1, and a shift of 1. The `l − 2F_n` relation applies only to the `oracle_bridge` packet, whose spinor follows the band around the ring.

Fix the `[INFO]` line and the docstring to say this.

**Tasks.**
1. **Production packet** (N = 24, σ = 4, k₀ = 2.902). Run with `--gradient fem` and `--gradient spectral`, and report per l:
   - λ_centroid at the first sample (step 1);
   - the median λ_centroid;
   - the phase residual.

   Expected step-1 values (audit oracle, Brillouin-zone fold):

   | Gradient | l = 0 | l = 1 |
   |---|---|---|
   | Spectral | ≈ −0.001 | ≈ +0.92 |
   | FEM | ≈ +0.0004 | ≈ −0.23 |

   The residual 0.08 from 1 comes from this packet itself:
   - the supercell wraps it (the envelope is about 1% of its peak at the cell edge);
   - the l = 1 core has non-zero amplitude.
2. **Clean control packet** (authorised addition, behind a `--clean` flag):
   - Grid: N = 48.
   - Envelope: `(r/σ)^|l| · exp(−r²/2σ²)`, with r the minimum-image distance from the cell centre.
   - Wave number: k₀ snapped to the grid, `k₀ = 2π·n1/N` with `n1 = 2·round(k₀_LSWT·N/(4π))`. That gives n1 = 22 and k₀ = 2.8798.
   - Run one step.
   - Targets, spectral: λ(l=0) = 0 within 1e-6, and λ(l=1) = 1 within 1e-4 (the audit oracle gives 0.99999997).
   - Report FEM alongside; expect a large bias at k·a ≈ 2.9.
3. Write `docs/dev/oam/BRIDGE.md` §6, "Production-kernel bridge". Include:
   - a table of both packets × both gradients;
   - the phase residuals;
   - one paragraph on what the production kernel does and does not validate.

**Acceptance:**
- The report is committed.
- `run_bridge_checks.py` still passes its existing B5.4 checks.
- The clean spectral l = 1 value meets 1e-4, or the report explains why not.
