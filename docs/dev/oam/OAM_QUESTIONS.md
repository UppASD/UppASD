# OAM questions and decisions

## 2026-09-27 — A1/A2

- Question: Should the static `hfield` be included in the LSWT Hamiltonian?
- Options: add its Zeeman contribution now, or retain the existing Hamiltonian and document the omission.
- Recommendation: retain the omission for A1/A2. The current LSWT Hamiltonian does not include `hfield`; UppASD now prints this explicitly whenever a nonzero field is present. A separate physics decision is required before adding the term.

- Question: What nonzero SA input represents a traceless tensor symmetric about the +z axis for the requested FM-chain Goldstone regression?
- Observation: `sa2tens` currently maps the three SA components to off-diagonal entries only, so it cannot express the axial tensor `diag(A,A,-2A)`. A nonzero off-diagonal tensor generally breaks the continuous spin-rotation symmetry and can gap the FM mode.
- Escalation: the added SA regression therefore checks the independent-neighbour-list Fourier phase, while the requested Goldstone assertion is held for maintainer guidance rather than guessing a physics convention.

## 2026-09-28 — B2 mesh regression

The baseline is `b60b762`; the after build is the current `[OAM-B1]` worktree. Every run used `OMP_NUM_THREADS=1`, `mseed 1`, `tseed 1`, and a two-step measurement phase. The selected skyrmion example is `examples/SpecialFeatures/SkyrmionLattice`. For `tests/triang2D`, the stock SLD initial phase crashes on the pinned baseline, so the committed overlay disables the SLD/lattice phase and retains the 20×20 periodic spin data for this mesh regression.

The skyrmion value is the iteration-zero `Skx num`; chirality is the scalar `<C>` magnitude. Relative change is `(after-before)/|before|`; when both values are zero it is reported as zero. The parenthesized chirality annotation is the estimated fraction of triangles in cells touching a periodic wrap (or, for HeisStripe, the fraction removed by the open-direction mesh rule).

| Case | Observable | Before | After | Relative change (wrap fraction for chirality) |
|---|---|---:|---:|---:|
| `tests/kagome` | skyrmion number | 0.00000000 | 0.00000000 | 0 |
| `tests/kagome` | chirality `<C>` | 0.00000000E+00 | 0.00000000E+00 | 0 (wrap 0.159722) |
| `tests/HeisStripe` | skyrmion number | 0.00000000 | 0.00000000 | 0 |
| `tests/HeisStripe` | chirality `<C>` | 3.25840689E-18 | NaN | NaN (open-direction removal 1.000000; `ntri` 2000 → 0) |
| `tests/triang2D` | skyrmion number | 0.00000000 | 0.00000000 | 0 |
| `tests/triang2D` | chirality `<C>` | 0.00000000E+00 | 0.00000000E+00 | 0 (wrap 0.097500) |
| `SkyrmionLattice` | skyrmion number | 112.00000000 | 112.00000000 | 0 |
| `SkyrmionLattice` | chirality `<C>` | 2.83021835E-19 | 2.83021835E-19 | 0 (wrap 0.015564) |

Skyrmion number is unchanged in all four cases. The B1 mesh diagnostic exposes a regression for `tests/HeisStripe`: its open y direction leaves no xy triangles (`Mesh2D: ntri= 0`), and chirality becomes NaN. This requires maintainer direction at gate G1 before B3; no source change is made in B2.
