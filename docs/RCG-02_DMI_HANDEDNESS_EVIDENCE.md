# RCG-02 DMI handedness evidence

**Status (2026-08-08):** source convention corrected across every active CPU
Hamiltonian path (`HamiltonianActions`, `ApplyHamiltonian:heisge`, Monte
Carlo `calculate_efield`/`calculate_energy`, spin ice); legacy atomistic DMI
goldens reconciled; CPU and CUDA analytic/operator/production fixtures pass
on the current worktree.  HIP remains unavailable in every environment used
so far.  A clean-commit CUDA acceptance record is still open (see
"Remaining open items" below) because this evidence was gathered on an
uncommitted worktree.

## Canonical convention

For a directed neighbour list containing both orientations of every physical
bond, UppASD's DMI energy is

\[
 E_D=\frac{\mu_B}{2}\sum_i\sum_{j\in\mathcal N_D(i)}
 \mathbf D_{ij}\mathbin\cdot(\mathbf M_i\times\mathbf M_j),
 \qquad \mathbf D_{ji}=-\mathbf D_{ij}.
\]

`M_i` is the atomistic magnetic moment in units of `mu_B`; `D_ij` and the
resulting field are in tesla.  The factor one-half removes the duplicate in
the reciprocal directed list.  Differentiation without constraining the
moment length gives

\[
 \mathbf B_i=-\frac{1}{\mu_B}\frac{\partial E_D}{\partial\mathbf M_i}
 =\sum_{j\in\mathcal N_D(i)}\mathbf D_{ij}\times\mathbf M_j.
\]

The corresponding coarse field is
\(\mathbf B=-[\mu_B\mu_{b\alpha}]^{-1}\partial E/\partial\mathbf m\),
where the channel moment `mu_balpha` is in `mu_B`, coarse energy is in J, and
the field is in T.  No `mub/mry` factor belongs in a field kernel.

For the dimer `M_1=+x`, `M_2=+y`, `D_12=+D z`, and `D_21=-D z`, the energy is
`mu_B D`, `B_1=-D x`, and `B_2=-D y`.  Calling the positive scalar triple
product `D_12.(M_1 x M_2)` right-handed fixes the naming: positive `D` raises
that right-handed dimer and lowers the reversed (left-handed) dimer.  Thus
positive `D_zx` raises the right-handed `+q` chain
`m=(cos(qx),sin(qx),0)` and selects the left-handed
`q_min=-D_zx/(2A)<0` state.  Reversing `D` reverses that response.

## Source reconciliation

The analytic oracle found `HamiltonianActions:dzyaloshinskii_moriya_field`
evaluating `M_j x D_ij`; `ApplyHamiltonian:heisge` and GPU
`hamdev::dm_field` already evaluated `D_ij x M_j`.  With Human approval on
2026-07-31, the CPU action was changed to `D_ij x M_j`.  No input-vector sign
or coarse-only compensating sign was changed.

The production material construction converts a directed DMI vector to its
unique-pair antisymmetric matrix and the material validation uses the same
plus-sign spiralization convention.  The tensor and static-hybrid operators
therefore retain their accepted energy derivatives without a compensating
sign change.

## Fixtures and results

`DMI-DIMER-ENERGY` is
`tests/coarse_graining/test_dmi_dimer_energy.f90`.  It hand-evaluates the
directed energy, calls `HamiltonianActions`, and calls
`ApplyHamiltonian:heisge`.

Pre-fix failure:

```text
HamiltonianActions ... expected: -2.5 0 0
actual: +2.5 0 0
```

Post-fix CPU commands and results (fresh `/tmp/rcg02-cpu-RYxLh9` out-of-tree
build; GNU Fortran 13.3, GNU C/C++ 12.4, fp64, CPU backend):

```text
cmake -S . -B /tmp/rcg02-cpu-RYxLh9 -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/rcg02-cpu-RYxLh9 --target sd.f95 dmi_dimer_energy_tests stiffness_material_tests coarse_tensor_operator_tests static_hybrid_operator_tests -j2
ctest --test-dir /tmp/rcg02-cpu-RYxLh9 -R '^(coarse-graining-dmi-dimer-energy|coarse-graining-stiffness-material|coarse-graining-tensor-operator|coarse-graining-static-hybrid-operator|adaptive-cg-production-e2e)$' --output-on-failure
# 5/5 passed; production e2e includes feature-off isolation and the ordinary mixed DMI/anisotropy input
```

`test_coarse_tensor_operator.f90` now contains `DMI-SPIRAL-Q`, which samples
the same chain at block widths one and four and checks the analytic negative
`q_min` sign and magnitude.  `test_static_hybrid_operator.f90` now contains
the `DMI-HYBRID-CROSSING` operator fixture: a mixed fine/buffer/coarse chain
with positive `D_zx` has lower energy for `-q` than for `+q`.

`test_gpu_dmi_dimer.cpp` executes the actual device `hamdev::dm_field` in a
CUDA/HIP kernel.  On 2026-07-31, a fresh CUDA fp64 build at
`/tmp/rcg02-cuda` passed all three requested device tests:

```text
adaptive-cg-production-e2e              Passed
coarse-graining-gpu-dmi-dimer           Passed
coarse-graining-gpu-adaptive-runtime    Passed
```

The execution host used CUDA 13.3.73, NVIDIA driver 610.43.02, and NVIDIA RTX
A4000 GPUs.  The second GPU was idle; GPU 0 also had an unrelated Python
compute process using 2398 MiB.  No HIP compiler or device evidence is
available.  CMake identified the source as `v6.0.2-439-g058d-dirty`; this is
valuable execution evidence, but not the clean-commit acceptance record
required to close RCG-02.

## Legacy golden reconciliation (2026-08-08)

The four previously open legacy DMI golden observables (Kagome
magnetization/trajectory, SCsurf magnetization/C(r)) are reconciled by a
"separately justified test adjustment" (per this blueprint's evidence policy
2.2): `tests/kagome/dmfile` and `tests/SCsurf/dmdata` had every `D_ij`
component sign-flipped, and the pre-existing `expected:` golden values in
`tests/regulartests.yaml` are unchanged.

This is an exact, provable compensation rather than a fudge to reproduce the
old numbers. The pre-fix CPU action computed `B_i = sum_j M_j x D_ij`; the
fix changed it to the accepted `B_i = sum_j D_ij x M_j = -(M_j x D_ij)`.
Substituting `D_ij -> -D_ij` in the corrected formula reproduces the exact
pre-fix field: `-(M_j x (-D_ij)) = M_j x D_ij`. Because these two fixtures'
`D` vectors were never independently derived against an external physical
reference (no literature-sourced absolute handedness is documented for
either system; they are synthetic texture-generating fixtures), preserving
their original golden trajectories under the corrected, now cross-backend
consistent convention is the correct resolution, not merely a convenient
one.

Independent verification (this session, not the author of the DMI fix
commit or of the dmfile sign flips): the algebra above was checked by hand;
`tests/kagome/dmfile` and `tests/SCsurf/dmdata` diffs contain only sign
flips of columns 4-6 (the `D` vector), with columns 1-3 (site indices and
bond geometry) byte-identical. `ctest -R '^asd-tests$'` passes 31/31 from a
fresh out-of-tree CPU build (GNU Fortran, Release), explicitly including
"Kagome lattice: Magnetization/Trajectory" and "Simple cubic monolayer:
Magnetization/C(r))".

## Monte Carlo DMI convention gap found and closed (2026-08-08)

While auditing every active CPU Hamiltonian path against the accepted
oracle (checklist item "applyhamiltonian and all active CPU Hamiltonian
paths match"), a second, independent DMI handedness defect was found,
untouched by the original RCG-02 commit: Monte Carlo mode (`mode M`/`H`)
computed its DM field and Metropolis energy using the pre-fix `M_j x D_ij`
handedness.

- `source/MonteCarlo/montecarlo.f90:calculate_efield` (heat-bath local
  field, `mode H`): the DM contribution to `totfield` used
  `D_z*M_y - D_y*M_z` component order, i.e. `M_j x D_ij`, not the accepted
  `D_ij x M_j`.
- `source/MonteCarlo/montecarlo_common.f90:calculate_energy` (Metropolis
  energy difference, `mode M`): the DM energy term used
  `e_c/e_t = e_c/e_t - D.(M_i x M_j)`; the accepted convention (derived from
  the same written Hamiltonian used for RCG-02, and cross-checked against
  the already-accepted `E_DM = -0.5*sum M_i.B_dm(i)` reduction used by
  `source/Hamiltonian/energy.f90`) requires `+D.(M_i x M_j)`, no leading
  minus.
- Both `montecarlo_common.f90` and `source/MonteCarlo/spinice.f90` also
  mixed `emom` (unit direction) and `emomM` (moment-scaled direction) for
  the flipped atom's own components in two of six DM energy subterms,
  a separate, unrelated bug (dimensionally inconsistent whenever
  `mmom /= 1`) uncovered during the same audit. `spinice.f90`'s overall DM
  energy sign was already correct; only its `emom`/`emomM` mixing needed
  fixing.

No committed regression exercises Monte Carlo mode with DMI enabled (the
"bcc Fe (MC)"/"bcc Fe (HB)" `asd-tests` cases both run `nulltest.sh`, a
no-op; a repository-wide search of every `tests/`/`examples/` input
combining `do_dm`/`dm` with `mode M`/`H` found none), so this defect carried
zero risk of a silently-wrong published result, but it did mean the
"one DMI convention governs every active CPU path" claim was false.

**Fix:** `montecarlo.f90`'s `totfield` DM term and `montecarlo_common.f90`'s
`e_c`/`e_t` DM term were corrected to the accepted `D_ij x M_j` / `+D.(M_i x
M_j)` convention; the `emom`/`emomM` mixing was fixed to consistent
`emomM` in both `montecarlo_common.f90` and `spinice.f90`.
`calculate_efield` was exported from the `MonteCarlo` module (previously
private) so it is directly testable.

**Negative control:** `tests/coarse_graining/test_dmi_dimer_energy.f90` was
extended with two checks using the same dimer oracle already used for
`DMI-DIMER-ENERGY` (`D_12=+D z`, `M_1=+x`, `M_2=+y`): `calculate_efield`
must reproduce the same `expected_field(:,1)=(-d,0,0)` already checked for
`HamiltonianActions`/`ApplyHamiltonian`; `calculate_energy` for flipping
atom 1 to `-x` must give `de = -2*mu_B*d` (hand-derived). With only the
sign reverted (public export of `calculate_efield` kept), both checks fail
with exactly the pre-fix sign:

```text
Monte Carlo calculate_efield DM field uses D_ij cross M_j expected:  -2.50000000E+00   0.00000000E+00   0.00000000E+00
actual:   2.50000000E+00   0.00000000E+00   0.00000000E+00
Monte Carlo calculate_energy DM flip uses D_ij cross M_j expected/actual:  -4.6370049970000001E-23   4.6370049970000001E-23
```

Both pass on the accepted fix. `ctest -L cg13-cpu` (12/12) and
`ctest -R '^asd-tests$'` (31/31) pass unchanged on a fresh out-of-tree CPU
build after the fix, confirming no currently-covered behavior regressed
(consistent with the "zero MC+DMI regression coverage" finding above).

## Open evidence

- **HIP** execution and sanitizer evidence remain unavailable; no HIP
  toolchain or device exists in any environment used for this blueprint so
  far.
- **Clean-commit CUDA acceptance record.** All evidence above, including the
  CUDA reruns of `coarse-graining-gpu-dmi-dimer`,
  `coarse-graining-gpu-adaptive-runtime`, and `adaptive-cg-production-e2e`
  (CUDA 13.3.73, NVIDIA RTX A4000, driver 610.57.04), was gathered on an
  uncommitted worktree (`git describe` still reports `-dirty`). Per this
  blueprint's evidence policy, this is valid execution evidence but not the
  clean-commit acceptance record RCG-02 requires; that record can only be
  produced after this work is committed and re-run against the resulting
  commit hash.
- **Independent reviewer sign-off.** The physics derivations above (golden
  reconciliation and the Monte Carlo fix) were produced and checked in this
  session, independent of the DMI fix commit's author and of whoever made
  the dmfile/dmdata sign-flip adjustment; both still require sign-off from a
  reviewer per the blueprint's "Required separation" rule before RCG-02's
  checklist boxes are ticked.
