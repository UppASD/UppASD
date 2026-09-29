# OAM questions and decisions

## 2026-09-28 — GA/C13 accepted

The maintainer accepts GA. The corrected A5 rerun establishes that positive
UppASD `D` maps to positive Fishman OAM and negative `D` reverses the sign,
with `D_UppASD/J = 2 D_paper/J` for the tested convention. The positive
`O_1,av(K)/hbar ≈ 0.228` agrees with the independent oracle. The approximately
3.4% difference from the published `0.236 hbar` is retained as a magnitude
discrepancy; it does not reopen the absolute-sign decision. B5.4 may proceed
once its remaining G3 prerequisites are satisfied.

## 2026-09-28 — B5.4 diagnostic run

Added `tests/SpinWaves/oam_lswt/run_bridge_checks.py`. It generates the
Haldane honeycomb `do_oam_lswt` reference, constructs small-amplitude
Holstein–Primakoff packets for `l=0` and `l=1`, runs zero-damping trajectories,
and extracts the precession frequency from the total moment output.

Command used:

```text
python3 tests/SpinWaves/oam_lswt/run_bridge_checks.py \
  --binary "$PWD/build-oam/bin/uppasd" \
  --workdir /private/tmp/oam_bridge_harness --nstep 1000 --force
```

The LSWT reference selected `k0=2.90207898`, `E=251.61986468 meV`, and
`F_n=0.09578900`. The measured packet frequency was `381.0094 rad/ps` versus
`382.2779 rad/ps` (`0.33%`, PASS). The initial packet projection was not
stable enough for the strict phase criterion, and the current per-sublattice
trajectory mesh returned `lambda_centroid(l=0)=0.00239` and an `l=1` shift of
`-0.22025`, rather than the B5.4 targets `F_n` and `+1`. The `F_n` bridge is
therefore not accepted; this is consistent with R6's warning that the
current per-sublattice mesh cannot recover the band-spinor derivative. The
generated run is retained outside the repository at the path above, and no
source or oracle change was made.

## 2026-09-28 — B5.4 packet/convention escalation

The oracle was extended with reciprocal-supercell packets, explicit particle
and hole terms, analytic spatial derivatives, and a WLS-versus-analytic
convergence check. Command:

```text
PYTHONPATH=tests/SpinWaves/oam_lswt python3 \
  tests/SpinWaves/oam_lswt/oracle_bridge.py \
  --convergence --k0 2.9 --sigma-k 0.18
```

For `n=12,16,20,24`, the analytic packet gives a stable `l=1` minus `l=0`
shift of approximately `0.99`, so the reciprocal packet and its angular
winding are working. However, the `l=0` value converges to approximately
`-0.206`, while the exact Wilson-loop result is `F_n=+0.09671`; the result is
consistent with the particle-only Fourier identity

```text
lambda = l + Im[u^dagger d u/dphi] = l - 2 F_n
```

for the documented `exp(+i k.r)` reconstruction and the Fishman definition
`F_n = -1/2 Im[u^dagger d u/dphi]`. The WLS gradient does not converge at this
large lattice wave vector, but the analytic-gradient result already exposes a
factor/sign conflict with C14's `l + F_n` target. Per C16, B5.4 stops here:
the maintainer must decide whether C14 should use the trajectory field's
`l - 2F_n` result, a conjugate/Fourier convention, or a separately defined
bridge observable. No Fortran bridge diagnostic was added and no oracle was
adjusted to force the C14 target.

The standalone all-site development oracle is now
`tests/SpinWaves/oam_lswt/oracle_bridge.py`. Its exact gauge-invariant mode is
available with `--reference-only`; for the confirmed positive-D mapping at
`k0=2.9`, it gives `F_n≈+0.09671`. The fixed-spinor control
packet evaluates to approximately zero for `l=0`, demonstrating that merely
changing the spatial mesh cannot generate `F_n`; the packet must include the
angularly varying band spinor. The oracle therefore reports both the fixed
control and the annular k-space packet, and leaves the contract-versus-Fourier
comparison (`l+F_n` versus `l-2F_n`) visible until the bridge normalization is
resolved.

## 2026-09-28 — R3

- Proposal (not implemented): add an explicit `oam_axis` contract key for restart-loaded or driven states, with the default frame axis equal to the normalized average moment at initialization. The C8 frame remains implicit in the current contract and oracle.

## 2026-09-28 — B5.4 held at GA

B5.1–B5.3 are prepared, but B5.4 is not started because its prompt requires
the maintainer's C13/GA sign decision before deriving the LSWT-to-trajectory
bridge mapping. A5 found that the positive-D Fishman sign is internally
consistent in UppASD, while neither tested sign reproduces the published
magnitude (the positive-D result is about 19.1% below it). Starting B5.4 now
would therefore bake an unresolved sign or normalization choice into the
bridge. The required HP derivation and bridge run remain a G3 follow-up after
GA.

## 2026-09-28 — B3

- Restatement: C1 defines the transverse complex field as `psi = m_x + i*m_y` in the single frame fixed at `oam_init`. In the +z Holstein–Primakoff convention this is proportional to the magnon annihilation field, so positive `lambda_L` denotes magnon OAM pointing along +z. It is the conjugate of the `m_x - i*m_y` convention used in parts of the literature.
- Restatement: C2 does not subtract a reference spin from the moment. The transverse field is made directly from the moment components, so its intensity is `|psi|^2 = 1 - m_z^2`; the reference state only fixes the frame.
- Restatement: C3 reports an origin-referenced value and a centroid-referenced value. The latter removes the packet-drift contribution `(R x P)_z`, but it does not remove the extrinsic envelope winding `l`.
- Origin dependence: with `P = sum Im[conj(psi)*grad(psi)] w`, changing the fixed origin from `x0` to `x0'` changes the numerator by `((x0 - x0') x P)_z`, hence `lambda_L(x0') - lambda_L(x0) = ((x0 - x0') x P)_z / N`. A stationary vortex has `P = 0`, so this difference vanishes for every origin.
- Gate G2 question: the C3 harness texture has an initial average transverse moment of approximately `(0.001048, 0.003601)` even though its intended frame is Cartesian +z. C8 requires `e_z` to follow the normalized average moment; applying that frame gives `lambda_L_origin = -5.059702` and `-5.123670` for the two shifts, while the independent oracle (which evaluates `m_x + i*m_y` in the unrotated Cartesian frame) expects `-4.835643` and `-4.896778`. Options are (a) update the fixture/oracle to rotate the prescribed texture into the C8 frame, or (b) treat this residual transverse mean as numerical and use the Cartesian frame, which would violate C8. Recommendation: maintainer decision is required; no tolerance or oracle change was made.

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

## 2026-09-28 — G1 follow-up on the empty open-direction mesh

The current OAM binary was rerun with the `tests/HeisStripe` input reduced to
two measurement steps and `do_chiral Y`. It reports
`Mesh2D: ntri= 0` and the explicit empty-mesh warning, while
`chirality.HeisStri.out` contains a finite zero row for all chirality
components and the magnitude. Thus the open one-dimensional case has no
defined 2D triangle chirality; the implementation's policy is to return zero,
not NaN. The earlier NaN entry above is retained as the historical B2 result,
but is not reproduced by the current binary. G1 can therefore be accepted by
classifying this case as a zero/N/A 2D-chirality regression, with no mesh or
OAM source change required.

## 2026-09-28 — G3 readiness summary

B5.1 is green in CTest (`oam-traj` and `oam-lswt`: 6/6). B5.2 is green for
`ell+1` and `w_area`; `shift_b` is N/A because pyswatter lacks the required
periodic minimum-image triangulation. B5.3 is green: the counter-clockwise
`m_x+i m_y` vortex gives positive trajectory `lambda`. B5.4 is not signable
until C14 selects the bridge normalization described above.

## 2026-09-28 — C14 bridge convention resolved

The maintainer resolves the former `l+F_n` versus `l-2F_n` conflict in favour
of the particle-field Fourier identity. With C1's `psi=m_x+i m_y`, UppASD's
`exp(+i k.r)` real-space reconstruction, and the particle-only HP packet,

```text
lambda_L_centroid(l,n,k0) = l - 2 F_n(k0)/hbar.
```

The factor two is the C9/Fishman normalization; the sign is fixed by the
documented C1 and Fourier conventions. The fixed-spinor packet remains a
control and is not expected to contain the Berry term. The reciprocal-
supercell packet with analytic spatial derivatives now validates the target
and the `l=1` minus `l=0` shift of `+1` through
`tests/SpinWaves/oam_lswt/oracle_bridge.py --validate`.

The current production trajectory mesh is per-sublattice, so its direct
`lambda_L_centroid` is retained as a discretisation diagnostic and is not
claimed to be `F_n`. B5.4/G3 can therefore sign the bridge convention using
the independent all-site oracle, while a future all-site production
observable remains a separate implementation item.

## 2026-09-29 — G1 and G3 accepted

The maintainer accepts both remaining gates. G1 closes with the documented
zero/N/A policy for the empty `HeisStripe` 2D mesh. G3 closes with B5.1–B5.4
accepted: the CTest harness is green, the two supported pyswatter comparisons
pass, `shift_b` remains N/A for the current pyswatter implementation, the
positive vortex sign passes, and the C14 particle-field bridge oracle passes
with `lambda = l - 2 F_n`.
