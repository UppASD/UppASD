# R4 — A5 rerun with corrected D mapping · Sonnet · `[OAM-A5b]`

Requires R1. **No source changes.**

A5 mapped d = 0.1 to D/J = 0.15. Check the normalisation from the K-point gap instead:
- UppASD gives `E1(K) - E2(K) = 6√3·D·S`: 84.84 meV at D = 0.15, i.e. 1.5588 in units of `4·ry_ev`;
- the paper gives `12√3|D|S`;
- so `D_UppASD = 2·D_paper`, and d = 0.1 means **D/J = 0.30**.

**Tasks.**
1. Confirm the gap relation from your own run output first, and report the numbers.
2. Rerun `docs/dev/oam/sign_pin/plus_D` and `minus_D` with every `dmfile` entry doubled (±0.15 → ±0.30).
3. Report F₁(K) and O₁,av(K) against the published 0.236ħ. The independent oracle expects about 0.228, i.e. −3%.
4. Note, without resolving it: the paper describes F₁ oscillating in sign over [0, K]; check whether yours does.

Append the result to `sign_pin/A5_report.md`. Stop for the maintainer's GA decision.
