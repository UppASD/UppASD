# CPU-HAM-11 — Consolidate Global Hamiltonian Energy Semantics

**Model:** Luna
**Status:** complete (2026-09-02)

## Dependency

CPU-HAM-10 complete.

## Purpose

Continue the longer-term architectural goal:

> global field and global energy semantics should come from the same canonical Hamiltonian term implementations wherever computationally sensible.

The code currently contains multiple historical/global energy routes.

This task should reduce convention drift without forcing unsuitable local-update algorithms through expensive global kernels.

## Non-goals

Do not:
- rewrite all Monte Carlo;
- remove specialized single-spin `ΔE` kernels merely for uniformity;
- force every term into `-1/2 m·B`;
- change Hamiltonian physics;
- merge unrelated legacy cleanup.

## A. Global energy call graph

Trace all current production global energy calculations.

Classify each as:
- measurement/output;
- diagnostics;
- transition/acceptance logic;
- MC global reference;
- other.

Identify which routines duplicate physics already present in `HamiltonianActions`.

Produce:

`docs/CPU_HAM_11_ENERGY_CONVENTION_MAP.md`

with a table:

| energy term | current routine | canonical field routine | same physics implementation? | can derive/reuse? | performance concern | action |
|---|---|---|---|---|---|---|

## B. Canonical term relationships

For every Hamiltonian term supported by current global measurements, document the exact field→energy relationship.

Examples:

### Reciprocal bilinear pair term

\[
E_{\rm pair}=-\frac12\sum_i \mathbf m_i\cdot\mathbf B_i^{\rm pair}.
\]

### Zeeman

\[
E_Z=-\sum_i \mathbf m_i\cdot\mathbf B_i^{\rm ext}.
\]

### Onsite terms

Use the actual canonical term expression.

Do not infer a universal homogeneity factor without proving it.

### Higher-order/multispin

Document term-specific semantics.

## C. Consolidate safe global terms

For terms where:
- the canonical field is already available at measurement time, or
- deriving energy from the canonical field is cheaper than another pair traversal,

redirect global measurement energy to canonical `HamiltonianActions` results.

Prioritize:
1. isotropic exchange;
2. DMI;
3. bilinear tensor exchange if supported;
4. simple onsite terms.

Do not duplicate pair traversals merely to preserve an old energy implementation.

## D. Measurement cadence

Preserve the execution contract:

- normal LLG steps may skip energy;
- measurement steps request canonical energy;
- backend-specific field representations remain transparent to energy semantics.

DIRECT, SPARSE and CONVOLUTION must produce the same global energy.

## E. Energy result decomposition

If users currently rely on term-resolved energies, preserve them.

Do not collapse everything into only a total `m·B`.

Where pair backends return a combined pair field but users need exchange and DMI separately, retain or introduce the minimum canonical term decomposition necessary.

Avoid independently reimplementing signs/factors.

## F. MC local `ΔE` exception matrix

Inventory specialized single-spin MC energy-difference kernels.

For each supported term, add or strengthen a parity fixture:

1. create a small state;
2. calculate canonical global total/term energy;
3. propose one spin change;
4. calculate canonical global energy after the change;
5. compare:
   `DeltaE_global`
   against:
   `DeltaE_MC`.

Use random non-collinear configurations where possible.

Include:
- exchange;
- DMI;
- anisotropy;
- other important MC-supported terms.

## G. Negative controls

At least for exchange and DMI:
- deliberately introduce a sign/factor error in the specialized MC `ΔE` path in a test mutation;
- prove the parity fixture fails.

Revert all mutations.

## H. Remove duplicated global formulas selectively

Only delete an old global energy implementation when:
- its replacement is validated;
- term-resolved outputs remain available;
- performance is equal or better for intended measurement cadence;
- no hidden caller depends on the old semantics.

Git history is sufficient archival storage.

Do not keep dead duplicate physics "for safety."

## I. Performance

For representative:
- Nd;
- Fe;
- J+D skyrmion;

measure energy-measurement cost before/after consolidation.

Energy output is not the hot LLG path, but avoid introducing pathological measurement overhead.

## J. Backend parity

For eligible systems compare term-resolved and total energies under:
- DIRECT;
- SPARSE if retained;
- CONVOLUTION.

Use the same canonical requested-energy interface.

## K. Documentation

Update developer documentation to state:

- `HamiltonianActions` is the global Hamiltonian convention authority;
- global measurement energies should reuse canonical term implementations;
- specialized local `ΔE` kernels are allowed for algorithms such as single-flip MC;
- such local kernels require parity tests against canonical global energy.

## L. Checklist

- [x] Global energy call graph complete.
- [x] Term-by-term field/energy relation documented.
- [x] Bilinear pair global energy consolidated.
- [x] DMI global energy consolidated where safe.
- [x] Onsite terms reviewed.
- [x] Term-resolved outputs preserved.
- [x] Normal LLG still skips energy when not requested.
- [x] DIRECT measurement energy passes.
- [x] SPARSE energy passes if backend retained.
- [x] CONVOLUTION energy passes.
- [x] MC specialized `ΔE` inventory complete.
- [x] Exchange MC `ΔE` parity fixture passes.
- [x] DMI MC `ΔE` parity fixture passes.
- [x] Other important MC term parity added where practical.
- [x] Exchange negative control is discriminating.
- [x] DMI negative control is discriminating.
- [x] Safe duplicate global energy formulas removed.
- [x] Measurement performance rechecked.
- [x] Developer architecture documentation updated.

## Commit

`CPU-HAM-11: consolidate Hamiltonian energy conventions`
