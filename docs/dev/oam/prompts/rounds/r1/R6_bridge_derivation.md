# R6 — Bridge derivation, no code · Opus · `[OAM-R6]`

Requires R1 and R4. Write `docs/dev/oam/BRIDGE.md` before anyone implements B5.4.

The trajectory mesh triangulates **each sublattice separately** (triangles never mix A and B), whereas pyswatter triangulates all sites together. For NA > 1, the trajectory λ therefore sees intracell structure only through the lever arm between sublattice centroids, `Σ_s |u_s|² (τ_s × k)_z`.

Derive, for a narrow packet in band n at k₀:
1. What `lambda_L_centroid` converges to on the per-sublattice mesh, in terms of `T_n(k₀)`, the sublattice positions `τ_s` and the envelope winding l.
2. Whether that equals `F_n(k₀)` (C10: the Berry phase / 4π, in the site-position Fourier convention of `setup_ektij`), and if not, what the correct comparison quantity is.
3. The HP mapping from `T_n(k₀)` to site `m_x + i m_y` in the C1 frame, including hole components and `sqrt(2S)`.

Check each result numerically with `oracle_honey.py` before claiming it. End with a recommendation: keep the per-sublattice mesh, switch to an all-site mesh, or offer both. **Stop for G3.**
