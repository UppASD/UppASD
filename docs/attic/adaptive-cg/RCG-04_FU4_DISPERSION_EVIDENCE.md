# RCG-04-FU4 coarse-block dispersion evidence

**Status:** complete for the CPU reference coarse tensor operator

**Date:** 2026-08-15

## Conclusion

The measured coarse-block dependence collapses onto the lattice symbols of the
implemented operator. In the controlled periodic sweep, exchange energy and
linearized field curvature agree with the exact discrete symbol to a maximum
relative residual of `8.881E-16`; the corresponding DMI residual is
`9.992E-16`. The fitted multiplicative energy factors are 1.000000 for both
exchange and DMI. Within double-precision and finite-period uncertainty, no
additional coarse scale factor is required.

This is an operator-level rate reconciliation: with zero damping, the
linearized LLG frequency contribution of a Fourier eigenmode is the production
field curvature multiplied by the channel gyromagnetic ratio. Therefore the
field-curvature ratios below are also frequency ratios for this benchmark.

## A. Exchange symbol from the production stencil

The source of truth is `source/CoarseGraining/coarsetensoroperator.f90`, where
`physical_forward_gradient` computes, for block coordinate n,

```text
G_p(n) = sum_d B[p,d] * (S(n+e_d) - S(n))
B = transpose(inverse(block_vectors))
```

`add_physical_gradient_transpose` is the periodic discrete transpose of this
forward difference. For a Fourier mode, define

```text
delta_d = exp(i*theta_d) - 1
theta_d = q . a_d
v_p(q) = sum_d B[p,d] * delta_d
```

The implemented exchange energy has mode symbol `lambda_A(q) = v(q)^H A v(q)`
and the exchange field derivative contains the adjoint factor, namely the
second-derivative symbol `D^H A D`. For the orthogonal benchmark with
`B[x,x]=1/h`, `q=(q,0,0)`, and `A[x,x]=A`,

```text
lambda_A(q) = A * |exp(i*q*h)-1|^2 / h^2
             = A * 4*sin(q*h/2)^2 / h^2
```

Relative to the continuum `A*q^2`, both the energy and exchange frequency
curvature have the exact factor

```text
R_A(qh) = [2*sin(qh/2)/(qh)]^2
```

## B. DMI symbol from the production stencil

The production loop evaluates

```text
E_D = V * sum_(k,p) D[k,p] * e_k . (S x G_p)
```

For the benchmark `D[3,1]=D`, a transverse helix
`S=(epsilon*cos(q*x), epsilon*sin(q*x), sqrt(1-epsilon^2))` gives

```text
E_D / V = D * epsilon^2 * sin(qh) / h
```

Thus the first-derivative lattice symbol and its continuum-normalized factor
are

```text
L_D(q) = sin(qh) / h
R_D(qh) = sin(qh) / (qh)
```

The field curvature obtained by differentiating the energy has the expected
factor of two, `2*V*D*sin(qh)/(h*mu_B*m_block)`; that factor cancels in the
measured/discrete comparison.

## C-E. Controlled periodic benchmark

`tests/coarse_graining/test_coarse_dispersion.f90` calls the production
`setup_coarse_tensor_operator` and `evaluate_coarse_tensor_operator` routines.
It uses a 64-block periodic chain, four orthogonal block sizes
`h = {0.1, 0.2, 0.4, 0.8} nm`, and six periodic modes
`m = {1, 2, 4, 8, 16, 24}`, with `qh = 2*pi*m/64`.

It runs independent exchange-only and DMI-only operators, measures the named
energy term, and extracts the transverse field curvature from the production
field. The same `qh` values occur at every `h`, so the block-size dependence
collapses when represented against `qh`.

| qh | exchange R_A | exchange measured/discrete | DMI R_D | DMI measured/discrete |
| ---: | ---: | ---: | ---: | ---: |
| 0.0982 | 0.999197 | 1.000000 | 0.998394 | 1.000000 |
| 0.1963 | 0.996791 | 1.000000 | 0.993587 | 1.000000 |
| 0.3927 | 0.987215 | 1.000000 | 0.974495 | 1.000000 |
| 0.7854 | 0.949641 | 1.000000 | 0.900316 | 1.000000 |
| 1.5708 | 0.810569 | 1.000000 | 0.636620 | 1.000000 |
| 2.3562 | 0.614991 | 1.000000 | 0.300105 | 1.000000 |

The four block-size rows at each `qh` are identical to the displayed
precision. The executable reports all 24 individual combinations and checks
both energy and field curvature against the discrete formulas at `2E-12`
relative tolerance.

Run used for the recorded result:

```text
cmake -S . -B build_rcg05g_cpu -DBUILD_TESTING=ON -DUSE_CUDA=OFF -DUSE_HIP=OFF
cmake --build build_rcg05g_cpu --target coarse_dispersion_tests -j2
OMP_NUM_THREADS=1 ./build_rcg05g_cpu/bin/coarse_dispersion_tests
```

Final summary:

```text
max exchange relative residual:   8.881784E-16
max DMI relative residual:        9.992007E-16
exchange energy multiplicative fit: 1.000000E+00
DMI energy multiplicative fit:       1.000000E+00
coarse dispersion tests passed
```

## F-G. Residual scale-factor test and decision

For every row, the measured energy was divided by the exact discrete-symbol
prediction; the same comparison was made for field curvature. A common
multiplicative fit was computed over all 24 block-size/wavelength cases. The
fits are unity to the printed precision and the residuals are at the
round-off floor. Consequently no empirical block-size renormalization,
missing normalization, or material-mapping factor is introduced by RCG-04-FU4.
If a future production frequency comparison leaves a residual after applying
`R_A(qh)` and `R_D(qh)`, that residual is evidence for a separate
implementation or material-mapping issue rather than a justification for
fitting another coarse scale factor.
