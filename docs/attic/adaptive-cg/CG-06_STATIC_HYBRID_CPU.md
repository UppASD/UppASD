# CG-06: static hybrid CPU coupling

## Scope

`StaticHybridOperator` is the double-precision CPU reference dispatch for a
fixed fine/coarse block mask.  It adds no selector, transition, reconstruction,
restart, GPU, finite-temperature, or μASD behavior.

The initial short-range input is a read-only unique-bond list with one
`3 x 3` energy matrix per bond.  A scalar diagonal represents Heisenberg
exchange and its antisymmetric part represents DMI.  Fine and buffer on-site
energies/fields are accepted as a paired input from the normal atomistic
Hamiltonian evaluation; coarse anisotropy and static external-field terms use
the accepted tensor operator.  The regular-grid FFT dipole remains a separate
all-grid owner and is rejected if embedded in this short-range dispatch.

## Static state and buffer construction

The caller supplies a Boolean fine seed for every canonical spatial block.
The operator computes the largest Euclidean radius among active atomistic
bond displacements,

\[
 r_{\max}=\max_b |\Delta\mathbf r_b|.
\]

For block-vector matrix \(H_b\), the conservative periodic dilation in
fractional block direction \(p\) is

\[
 w_p=\left\lceil r_{\max}
       \left\|(H_b^{-1})_{p,:}\right\|_2\right\rceil+w_{\rm safety}.
\]

This bounds a real-space ball even for skew block vectors.  Dilation uses
minimum-image block-coordinate distance and therefore wraps periodic
boundaries.  Fine seeds retain `FINE` state; other dilated blocks are
`BUFFER`; all remaining blocks are `COARSE`.  The safety dilation is a
nonnegative, static setup parameter.

## Non-overlapping ownership

The short-range ownership rules are executable data in the operator:

1. A unique atomistic bond is owned exactly when it is active and at least one
   endpoint is a real fine/buffer atom.
2. An exchange tensor density \(A_{pq}G_p\mathbf m\cdot G_q\mathbf m\) is
   owned only if the base block and every block used by both physical-gradient
   stencils are coarse.
3. A DMI tensor density is owned only if the complete corresponding physical
   gradient stencil is coarse.
4. Fine/buffer on-site terms are taken only from atomistic inputs.  Coarse
   on-site terms are taken only from the tensor operator.
5. FFT dipole energy and field are evaluated once by their existing uniform
   grid owner, outside this short-range sum.

Consequently no atomistic bond or tensor density is evaluated by both owners.
Tensor fields retain all derivative reactions on the blocks in each owned
stencil.

The dispatcher first prolongates the full coarse state with CG-05's normalized
trilinear operator.  Real fine/buffer directions replace the prolongated values
in their blocks.  Owned bonds are then evaluated once.  Fields on real atoms
are returned directly; fields on coarse ghost atoms are accumulated with the
exact moment-weighted derivative adjoint of normalized prolongation.  Only
active coarse-block reactions are returned.

Both all-fine and all-coarse configurations execute this same dispatch.  The
former owns every active atomistic bond and no tensor terms.  The latter owns
no atomistic bonds and reduces exactly to the accepted CG-04 tensor reference.

## Validation

`tests/coarse_graining/test_static_hybrid_operator.f90` covers:

- all-fine energy/field identity with a direct atomistic exchange, DMI, and
  static on-site reference;
- all-coarse identity with the accepted tensor energy and field;
- exact complementary bond ownership and interaction-radius plus safety
  dilation, including periodic wrapping;
- uniform and constant long-wave spiral interface patch tests;
- tangent finite differences on both a real fine spin and a coarse ghost
  reaction degree of freedom;
- translated periodic domain-wall-pair and skyrmion profiles across fixed
  interfaces at block size one;
- decreasing long-wave interface energy error when the block width is refined
  from four cells to two.

The CG-02 through CG-06 CPU tests pass together through CTest.  No file under
`source/Multiscale` (the protected μASD implementation) is used or modified by
CG-06.
