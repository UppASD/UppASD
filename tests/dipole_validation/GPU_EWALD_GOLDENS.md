# WP4b GPU Ewald oracle fixtures

`gpu_ewald_goldens_v1.hpp` is generated only from
`periodic_ewald_reference.py` using its converged direct-particle Ewald
evaluator. It contains full periodic geometry, moments, dimensionless fields,
and total dimensionless energy. It is never generated during CTest and is not
derived from Builder A, an uploaded GPU tensor, `directConvolution`, or legacy
`do_dip`.

After a deliberate independent-oracle review, regenerate and inspect it:

```bash
python3 tests/dipole_validation/generate_gpu_ewald_goldens.py \
  > tests/dipole_validation/gpu_ewald_goldens_v1.hpp
python3 tests/dipole_validation/test_periodic_ewald_reference.py
```

The v1 generator fixes tin-foil 3D periodic Ewald, NA=1, fp64 oracle
convergence (`tol=1e-13`, `max_shell=16`), and C++/GPU n1-fast cell order. An
approved convention change requires a new `v2` fixture; do not silently
rewrite v1.
