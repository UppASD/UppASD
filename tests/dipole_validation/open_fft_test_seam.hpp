#pragma once

// Test-only contract requested from Terra by WP10.1.  This header is not a
// production GPU API and intentionally has no implementation in this change.
// Terra may adapt the spelling, but the observable values and storage order
// below must remain available to the red test.

#include <array>
#include <cstddef>
#include <vector>

namespace luna_open_fft_test {

struct Input {
   std::array<std::size_t, 3> active_grid{};
   std::array<std::size_t, 3> fft_grid{};
   unsigned int basis = 0;
   unsigned int ensembles = 0;
   // Column-major [C1 C2 C3], followed by basis offsets in Cartesian order.
   std::array<double, 9> primitive_vectors{};
   std::vector<std::array<double, 3>> basis_offsets;
   // alpha + 3*(macro + active_macros*ensemble), with macro=a+NA*cell.
   std::vector<double> active_moments;
   // Explicit test tensor in fft_cell + fft_cells*kernel_batch order.
   // The seam deliberately has no open physical-kernel builder.
   std::vector<double> padded_real_kernel;
};

struct Result {
   std::size_t active_cells = 0;
   std::size_t active_macros = 0;
   std::size_t fft_cells = 0;
   std::size_t field_batches = 0;
   // These are host copies of the exact test seam buffers.
   std::vector<double> packed_moments;       // fft_cell + fft_cells*field_batch
   std::vector<double> padded_fields;        // same layout as packed_moments
   // Active fields: alpha + 3*(macro + active_macros*ensemble).
   std::vector<double> active_fields;
   // One dimensionless -1/2 sum(M dot B) per ensemble.
   std::vector<double> dimensionless_energy;
   // Independent allocation-inventory checks for the test descriptor.
   std::size_t persistent_bytes = 0;
   std::size_t persistent_inventory_bytes = 0;
   std::size_t construction_bytes = 0;
   std::size_t construction_inventory_bytes = 0;
   std::size_t workspace_bytes = 0;
   std::size_t total_bytes = 0;
};

// Terra must provide this test-only adapter around the open layout/runtime.
// It must not call periodic Ewald construction and must not apply the Tesla
// prefactor to active_fields or dimensionless_energy.
Result run(const Input& input);

}  // namespace luna_open_fft_test
