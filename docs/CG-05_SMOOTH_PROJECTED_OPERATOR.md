# CG-05: smooth projected CPU reference operator

## Scope and data convention

`SmoothProjectedOperator` is a double-precision, periodic, CPU reference
operator. It does not own or rewrite atomistic interaction data. Coarse
variables have shape `(spin, dynamic_channel, spatial_block)`, and atom arrays
retain the caller's atom index.

The setup input `atom_fractional_block_coordinate(:,a)` uses the regular block
grid as its coordinate system. Integer coordinate `(i,j,k)` is the coarse
variable at that block centre. Coordinates may cross a periodic face; stencil
indices wrap independently in each grid direction. Open and other boundary
rules are rejected rather than inferred.

For

\[
 x_a=b_a+\xi_a,\qquad 0\leq \xi_{ap}<1 ,
\]

the eight shape weights are

\[
 w_{a,\delta}=\prod_{p=1}^3
 \begin{cases}
 1-\xi_{ap},&\delta_p=0,\\
 \xi_{ap},&\delta_p=1,
 \end{cases}
\]

with the node \(b_a+\delta\) reduced modulo the block-grid extent. Thus
\(w\geq0\) and \(\sum_\delta w_{a,\delta}=1\), including collapsed dimensions
and wrapped faces.

## Normalization and adjoint field

An atom uses the dynamical channel recorded by `BlockTopology`:

\[
 v_a=\sum_I w_{aI}M_{I,c(a)},\qquad
 s_a=\frac{v_a}{r_a},\qquad r_a=|v_a|.
\]

The setup has a positive normalization floor. An interpolant at or below that
floor is rejected because its direction and derivative are undefined. No
tangent-plane approximation is used. The implemented Jacobian is

\[
 \frac{\partial s_a}{\partial M_{I,c(a)}}=
 \frac{w_{aI}}{r_a}\left(\mathsf I-s_as_a^\mathsf T\right).
\]

Let \(\mu_a\) be the atomic moment and let \(\mu_{Ic}\) be the sum of moments
assigned to block/channel \((I,c)\). For atomistic fields in tesla, the coarse
field is

\[
 B_{Ic}=\sum_{a:c(a)=c}
 \frac{\mu_a}{\mu_{Ic}}\frac{w_{aI}}{r_a}
 \left(\mathsf I-s_as_a^\mathsf T\right)B_a .
\]

This satisfies the weighted virtual-work identity

\[
 \sum_a\mu_a B_a\mathbin{\cdot}\delta s_a
 =\sum_{Ic}\mu_{Ic}B_{Ic}\mathbin{\cdot}\delta M_{Ic}.
\]

Because the prolongated spin is normalized, radial atomistic field components
are a gauge and disappear from this adjoint. Comparisons with the tensor
operator therefore compare physical tangent fields (or torque), not an
arbitrary radial component.

## Supported atomistic bilinear energy

Each physical bond is supplied once as atom indices `(i,j)` and a read-only
matrix \(K_{ij}\), in joules:

\[
 E_{ij}=-s_i^\mathsf T K_{ij}s_j .
\]

This single representation covers scalar Heisenberg exchange
\(K=J\mathsf I\), antisymmetric DMI matrices, symmetric anisotropic exchange,
and their sum. The atomistic fields are evaluated from both matrix actions,
\(K_{ij}s_j\) and \(K_{ij}^\mathsf T s_i\), then passed through the adjoint
above. Self-bonds, nonmagnetic endpoints, malformed arrays, and non-finite
couplings reject cleanly. Source bond indices and matrices are `intent(in)`;
the regression test also snapshots them to catch mutation.

## Numerical evidence

`test_smooth_projected_operator.f90` checks:

- partition of unity, periodic wrapping, collapsed grid dimensions, preserved
  atom order, and distinct block/channel identities;
- the block-size-one identity;
- projected energy derivatives and the standalone weighted-adjoint identity
  by centred finite differences;
- scalar, accepted-handedness DMI, and generic nonsymmetric bilinear matrices
  without source-array mutation;
- a 256-site nearest-neighbour spiral at block widths 2, 4, and 8.

For the spiral, \(A=J/(2a)\) is used by the CG-04 tensor reference. The smooth
projected excess energies are mesh-independent to better than 0.01% in the
fixture. Their differences from the tensor stencil are approximately 0.015%,
0.075%, and 0.317% as \(qh\) grows. On a weakly phase-modulated spiral, tangent
field differences are approximately 0.04%, 0.24%, and 1.54%. These are the
expected finite-\(qh\) interpolation/stencil errors; both descriptions share
the same long-wavelength limit.

The test-only piecewise-constant projection gives excess energies in the ratio
1:2:4 for those increasing block widths (up to the finite-\(qh\) correction).
It concentrates exchange rotation at block faces and therefore has an
incorrect leading dependence on block size. It is intentionally absent from
the production operator API.
