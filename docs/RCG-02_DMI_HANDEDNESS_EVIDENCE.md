# RCG-02 DMI handedness evidence

**Status:** source convention corrected; CPU and CUDA analytic/operator/
production fixtures pass.  HIP is unavailable.  Legacy atomistic DMI goldens
need independent reconciliation.

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

## Open evidence

`ctest -R '^asd-tests$'` has four changed legacy DMI golden observables:
Kagome magnetization/trajectory and SCsurf magnetization/C(r).  These are
not accepted as a pass; their old values reflect the previous CPU DMI
handedness.  They require independent physics review and regenerated,
reviewed references or a separately justified test adjustment.  HIP execution
evidence and independent reviewer sign-off are also open.
