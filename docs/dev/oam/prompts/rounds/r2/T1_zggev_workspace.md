# T1 — Fix the `zggev` workspace · Sonnet · `[OAM-T1]`

`fallback_bosonic_diag` (`diamag.f90:714`) allocates `rwork(hdim)`. LAPACK `zggev` requires `RWORK(8*N)`. The forced-fallback test `lswt-diagonalizer-fallback` aborts on Linux/glibc with `malloc(): corrupted top size`; it only passed on macOS by luck.

**Tasks.**
1. Allocate `rwork(8*hdim)`, keeping the `memocc` accounting.
2. Check every other LAPACK call in `diamag.f90` and `chern_number.f90` against its documented minimum workspace (`zheev`, `zpotrf`, `ztrtri` and others). Fix any that fall short, and list each call with required and allocated sizes in the report.

**Acceptance:** on Linux, `ctest -L "oam"` gives 6/6. `valgrind` or `-fsanitize=address` is clean for `lswt-diagonalizer-fallback`.
