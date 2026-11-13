<!--
SPDX-FileCopyrightText: 2026 CSC – IT Center for Science

SPDX-License-Identifier: Apache-2.0
-->

# TODO

## Correctness

- [ ] Build and run `test-host` and `test-device` with the new test data, which includes interpolated neighbours of higher and lower index.

## Upstream UppASD issues

Report only. Submodules are not modified.

- [ ] `multiscaleInterpolateInterfaces()` updates in place while reading interpolated neighbours. The result depends on the atom order and races with OpenMP.
- [ ] `GpuInterpolation::interpolateFirst()` and `interpolateSecond()` compute the block count with `threads.y` instead of `threads.x`, launching 256 times too many blocks.
- [ ] `createPaddingInterpolationWeights()` adds a matrix entry for a corner with index 0, with a stale `tmpWeights(j)`.

## Benchmarking

- [ ] Add `--repeat N` to launch the kernel several times. One launch includes warm-up cost.
- [ ] Add options for the fraction of interpolated atoms and the neighbours per row. Typical UppASD values are unknown.

## Implementation

- [ ] Autotune the thread block dimensions.
- [ ] Review the memory ordering of the moment arrays for coalesced access.
- [ ] Profile the kernel to select the next optimizations.
