# RCG-04D E2E-MOVING-OFF-FINE: feature-off half

**Purpose:** the ordinary (AdaptiveCG-disabled) half of the RCG-04D moving
feature-off/all-fine parity pair. Paired fixture:
`../moving_all_fine`. Driven together by
`tests/coarse_graining/run_moving_off_fine.py`. Full analysis, defect
finding, and evidence: `docs/RCG-04_MOVING_E2E_EVIDENCE.md` (RCG-04D
section).

**Why a new fixture instead of reusing `../feature_off`:** the existing
`feature_off`/`static_all_fine` pair (RCG-04A evidence, section 3.1/3.5) uses
a uniform, exactly-aligned `momfile` -- an exact stationary fixed point of
isotropic exchange, with zero initial LLG torque. That pair is retained
unmodified as a zero-torque smoke test. It cannot support a moving-dynamics
parity claim: a frozen all-fine path and a genuinely correct one are
indistinguishable when the state never moves either way, which is exactly
what let a real defect (see below) go undetected until this slice.

**Initial state:** a conical spin spiral (RCG-04B
`moving_state_generator.conical_spiral_state`), not a planar spiral (the
degenerate construction `initmag_spin_spiral` uses, flagged in the RCG-04A
evidence and independently confirmed zero-torque by
`test_torque_oracle.py::ConicalSpiralDegeneracyTests.test_planar_spiral_zero_torque`).
Parameters (identical to `../moving_all_fine`; also recorded in
`GENERATOR_MANIFEST.json`):

```python
from moving_state_generator import Geometry, conical_spiral_state
geometry = Geometry(na=2, n1=6, n2=2, n3=2,
                     basis=((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)))
state = conical_spiral_state(
    geometry, cone_angle_deg=40.0, turns=1, axis=(0.0, 0.0, 1.0),
    phase_deg=0.0, modulation_cell_axis=0,
    moment_magnitude=2.23, landeg=2.0,
)
```

`cone_angle_deg=40` is comfortably inside the generator's non-degenerate
range (>5 degrees from 0/90/180); `run_moving_off_fine.py` independently
recomputes (not reads back from production) the initial exchange torque via
`torque_oracle.py` before accepting either run and asserts it exceeds a
documented floor -- observed `max_torque=rms_torque~0.6725` (uniform across
atoms, an exact-conical-spiral property), `max_field_misalignment_deg~3.17`.

`momfile` here and in `../moving_all_fine` are byte-identical (verified by
`run_moving_off_fine.py`, not merely by convention): both consume the same
generator output, and this file was produced by, and only by, the Python
snippet above. `posfile`/`jfile` are the shared
`tests/coarse_graining/e2e/{posfile,jfile}` every other fixture in this
directory uses (exchange-only Hamiltonian, no DMI/anisotropy/dipole/external
field -- the only intentional physical difference from `../moving_all_fine`
is AdaptiveCG enablement).

**Integration parameters** (identical in both fixtures): `damping 0.05`,
`timestep 1.0e-16`, `Nstep 50`, `do_tottraj Y`/`tottraj_step 5` (11 sampled
steps, `0,5,...,50`). Chosen from an out-of-band convergence check (not
included in this fixture, recorded in the evidence document): halving the
timestep (`Nstep 100`, `timestep 0.5e-16`, same total simulated time) changes
the final feature-off state by only ~0.003 degrees against a total motion of
~7.1 degrees over the accepted run -- comfortably inside a stable, converged
integration regime, not a large-angle/near-instability one.

**A real production defect was found and fixed while building this
fixture:** before the fix, `../moving_all_fine` (AdaptiveCG all-fine) left
every spin essentially frozen (displacement ~1e-13 rad after 50 steps,
scaling linearly with step count at a rate ~12 orders of magnitude below
this fixture's own dynamics) for this exact same initial state, while this
`feature_off` fixture evolved normally. Root cause:
`source/CoarseGraining/adaptivecgproduction.f90:adaptive_cg_cpu_step` called
`llg_rhs(direction, field, Landeg(1), lambda1_array(1), rhs)` for every
atomistic-fine atom -- passing only the dimensionless Landé g-factor as the
LLG "gamma" argument, omitting the physical gyromagnetic-ratio constant
`gama = 1.760859644e11` (`source/Parameters/constants.f90`) entirely, and
hardcoding atom index 1 instead of the actual atom's own `Landeg`/damping.
The coarse-block path (`coarse_llg_rhs`, using
`operator%channel_gamma_per_t_s`) and the GPU path (`gpuAdaptiveRuntime.cpp`,
using `kernels.gammaPerTs`) were not affected -- both already use a properly
physically-scaled rate. No existing fixture (all uniform/zero-torque) could
have caught this: with zero torque, `rhs=0` regardless of the missing
factor, so feature-off and all-fine trivially "agreed" by both doing
nothing. Fixed by passing `gama*Landeg(atom)`/`lambda1_array(atom)` at all
three CPU call sites. Full before/after evidence, and confirmation that the
fix introduces no regression in the existing `cg13-cpu`/`asd-tests` suites,
is in the RCG-04D evidence document.
