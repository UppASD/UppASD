# Coarse-graining validation suites

The release suites are CTest label sets. After configuring a build with
`BUILD_TESTING=ON`, run:

```sh
ctest --test-dir <build> -L cg13-cpu --output-on-failure
ctest --test-dir <build> -L cg13-cuda --output-on-failure
ctest --test-dir <build> -L cg13-hip --output-on-failure
```

The backend-specific command is meaningful only for the matching build. The
CPU label contains the analytic/reference operators, production executable
cases, and setup-rejection matrix. CUDA and HIP add the same production
matrix, backend runtime parity tests, and backend-specific dipole checks.

The packaging dependency check is independent of a built executable:

```sh
python3 tests/coarse_graining/audit_fixture_dependencies.py
ctest --test-dir <build> -R adaptive-cg-fixture-dependencies --output-on-failure
```

It resolves each harness case and every input file named by its `inpsd.dat`
against `git ls-files`; generated restart/output files are not source
dependencies. A failure means the executable suites are not reproducible from
a clean clone.

Run `ctest --test-dir <build> -N` to inspect what the current configuration
provides. Hardware timing is intentionally not a pass/fail CTest: use the
`gpu_adaptive_runtime_benchmark` binary and archive its complete output as
described in `docs/CG-13_RELEASE_VALIDATION.md`.
