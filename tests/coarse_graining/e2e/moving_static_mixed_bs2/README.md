# RCG-04F: E2E-MOVING-STATIC, block_size_x=2 (spatial-refinement baseline)

Same construction as `../moving_static_mixed_bs1` (identical `momfile`,
`ncell 24 2 2`, physical/integration parameters), but `block_size_x=2` with
FINE seed blocks 1-3 (of 12 total blocks). `buffer_width_blocks_x =
ceil(1.6583/2 - 64*eps) + 0 = 1` block, giving interface/buffer blocks
{4,12} (2 blocks) and COARSE blocks {5..11} (7 blocks) -- see
`static_topology_oracle.py`. This covers the *same physical* 48-fine/
32-interface/112-coarse-atom partition as `../moving_static_mixed_bs1`,
discretised at half the block resolution along x, giving a genuine
fixed-physical-evolution spatial refinement pair (`bs2` -> `bs1`) rather
than a re-run of a different experiment.
