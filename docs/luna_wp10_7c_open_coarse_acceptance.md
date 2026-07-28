# WP10.7c Luna acceptance — uniform divisible `OPEN_FFT` coarse blocks

Date: 2026-07-28
Evaluated revision: `8d101f0` plus the acceptance-harness changes in the
working tree.

Decision: **GO for the CUDA fp64 production scope.**

The independent host physics gate and CUDA production coarse-block matrix
pass.  HIP numerical execution remains untested because no HIP toolchain or
device was available.

## Authority and checks

Expected coarse values were formed as an explicit finite point-pair
projection:

```text
Kcoarse(A,B) = 1/(nA nB) sum(i in A,j in B) Kpoint(i,j)
M_A          = sum(i in A) m_i
B_A          = sum_B Kcoarse(A,B) M_B
E            = -1/2 sum_A M_A dot B_A
```

The point self pair is excluded; distinct intra-block pairs remain on the
coarse diagonal.  No CPU `do_dip=2`, periodic Ewald result, or production
builder output is used as expected data.

| Gate | Result | Evidence |
|---|---:|---|
| Independent Python open oracle | PASS | 9/9 tests |
| Host builder vs direct finite point sums | PASS | max `7.815970093361102e-14` |
| Explicit projection, non-cubic geometry | PASS | max `4.441e-16` |
| Explicit projection, skew two-basis geometry | PASS | max `9.770e-15` |
| Uniform fine moments: projected field | PASS | max `4.441e-16` |
| Uniform fine moments: projected energy | PASS | max `1.110e-16` |
| Projected diagonal, block `(2,1,1)` | PASS | `Kxx=1`, `Kyy=Kzz=-0.5` exactly |
| Block-one recovery | PASS | field `<3e-14`, energy `<3e-14`; observed `6.661e-16`, `1.388e-17` |
| Nonuniform convergence | PASS | see table below |
| Host rejection: partial x/y/z edges | PASS | deterministic `invalid_argument` |
| CUDA production coarse matrix | PASS | CTest `2/2` |
| CUDA Compute Sanitizer | PASS | `ERROR SUMMARY: 0 errors` |
| HIP production matrix | UNTESTED | no HIP toolchain/device on runner |

### Nonuniform convergence table

The fine reference is the direct finite point Hamiltonian.  Coarse moments are
block sums, and the reported field error is the maximum difference at member
atoms; energy is absolute error per fine atom.

| Block width | max field error | energy error / atom |
|---:|---:|---:|
| 4 | `5.144e-01` | `4.058e-04` |
| 2 | `4.377e-01` | `2.833e-04` |
| 1 | `6.661e-16` | `1.388e-17` |

The required monotonic improvement and exact block-one limit pass.

### CUDA production results

The CUDA build exercised the production seam for NA=1/2, M=1/4, non-cubic
and skew geometries, block-one recovery, impulses, scatter, energy, and
partial/population rejection.  The reported maxima were:

| CUDA check | Maximum |
|---|---:|
| Block-one/basis production field | `2.274e-13` |
| Block-one/basis production energy | `2.757e-10` |
| Coarse projected production field | `9.948e-14` |
| Coarse projected production energy | `3.411e-13` |
| Energy derivative residual | `1.168e-9` absolute |

The energy derivative was `-2.0000000011677344` against negative field `-2`.
Both registered CUDA tests passed in `2.91 s`.

## Rejection evidence

The production Fortran preflight was exercised with `OPEN_FFT` for eight
invalid inputs; all rejected before GPU initialization:

| Case | Result |
|---|---:|
| Partial x edge (`3x1x1`, block `2x1x1`) | PASS |
| Mixed boundary conditions | PASS |
| Legacy `do_dip=1` | PASS |
| Monte Carlo mode | PASS |
| Non-default surface | PASS |
| Nonzero mesh | PASS |
| Nonzero alpha/rcut | PASS |
| Non-default tolerance | PASS |

The C++ pre-allocation path additionally checks regular coarse quotient,
basis-resolved macro count, every declared population, map range, and map
histogram agreement.  The test seam now exercises population mismatch and
partial edges independently on all three axes.

## Remaining limitation

Run the equivalent suite on an AMD/HIP host:

```text
ctest --test-dir <build> --output-on-failure \
  -R 'dipole-open-fft-layout|dipole-open-fft-coarse-e2e'
```

CUDA production enablement is accepted; HIP parity remains a WP10.8 follow-up.
