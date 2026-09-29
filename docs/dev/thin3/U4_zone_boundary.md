# U4 — Document the zone-boundary limitation · Sonnet · `[OAM-U4]`

**Docs only.** Run after U1.

**Why.** A lattice field has no unique gradient for content at or across the Brillouin-zone boundary: momentum is only defined modulo G.
- Under the U1 rule, a mode exactly at K ties three images at 120°, and their mean is zero. Its derivative is therefore dropped.
- A packet centred on K or M straddles the boundary, and no single-valued fold is correct for it.

This matters for honeycomb Dirac-point magnons.

**Tasks.**
1. **`OAM_IO.md`, near the `oam_gradient` entry:**
   - Spectral trajectory OAM is valid for fields whose content lies strictly inside the first Brillouin zone of each sublattice's Bravais lattice.
   - K- and M-centred packets are out of scope for both gradient methods.
2. **`OAM_QUESTIONS.md`:** log a proposal, not implemented, for a carrier-referenced gradient:
   - Idea: `ψ = e^{i k₀·r} φ` with `∇ψ = e^{i k₀·r}(i k₀ φ + ∇φ)`, where `k₀` comes from an input key or the spectral peak.
   - Why it helps: the carrier term cancels in λ_centroid (it is the drift term), so only `∇φ` needs to be accurate.
   - Open questions:
     - how k₀ is chosen during the dynamics;
     - how it interacts with the C8 frame;
     - how to test it (suggested: a K-centred packet with l = 1 on the honeycomb lattice).
3. **`BRIDGE.md`:** one sentence noting that the bridge packet (k₀ ≈ 2.9 < |K| = 4.19) lies inside the zone and is therefore in scope.

**Acceptance:** docs only; both harnesses still ALL PASS; `OAM_IO.md`, the C17 text and `OAM_QUESTIONS.md` are consistent.
