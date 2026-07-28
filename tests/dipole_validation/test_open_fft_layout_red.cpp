// WP10.1 red tests.  This file is intentionally not part of the current CMake
// targets: it becomes buildable/linkable when Terra supplies
// open_fft_test_seam.hpp's luna_open_fft_test::run implementation.

#include "open_fft_test_seam.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
   if(!condition) throw std::runtime_error(message);
}

std::size_t cell(std::size_t x, std::size_t y, std::size_t z,
                 const std::array<std::size_t, 3>& grid) {
   return x + grid[0] * (y + grid[1] * z);
}

void requireZeroPadding(const luna_open_fft_test::Result& result,
                        const luna_open_fft_test::Input& input) {
   for(std::size_t z = 0; z < input.fft_grid[2]; ++z)
      for(std::size_t y = 0; y < input.fft_grid[1]; ++y)
         for(std::size_t x = 0; x < input.fft_grid[0]; ++x) {
            const bool active = x < input.active_grid[0] && y < input.active_grid[1] && z < input.active_grid[2];
            if(active) continue;
            const std::size_t q = cell(x, y, z, input.fft_grid);
            for(std::size_t field_batch = 0; field_batch < result.field_batches; ++field_batch)
               require(result.packed_moments[q + result.fft_cells * field_batch] == 0.0,
                       "OPEN_FFT packed moment padding is nonzero");
         }
}

luna_open_fft_test::Input simpleInput() {
   luna_open_fft_test::Input input{};
   input.active_grid = {2, 1, 1};
   input.fft_grid = {4, 1, 1};
   input.basis = 2;
   input.ensembles = 4;
   input.primitive_vectors = {1.0, 0.0, 0.0,
                              0.0, 1.0, 0.0,
                              0.0, 0.0, 1.0};
   input.basis_offsets = {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}};
   input.active_moments.resize(3 * 4 * 4);
   for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
      for(std::size_t macro = 0; macro < 4; ++macro)
         for(unsigned int component = 0; component < 3; ++component)
            input.active_moments[component + 3 * (macro + 4 * ensemble)] =
               0.11 + 0.17 * ensemble + 0.07 * macro + 0.03 * component;
   return input;
}

void runLayoutAndBatchTests() {
   const auto input = simpleInput();
   const auto result = luna_open_fft_test::run(input);
   require(result.active_cells == 2 && result.active_macros == 4, "active storage is not Nactive times basis");
   require(result.fft_cells == 4 && result.field_batches == 24, "unexpected padded/batch counts");
   require(input.active_moments.size() == 3 * result.active_macros * input.ensembles,
           "active macro input does not have Nactive entries");
   require(result.packed_moments.size() == result.fft_cells * result.field_batches,
           "packed moment buffer has the wrong padded extent");
   requireZeroPadding(result, input);
   for(unsigned int ensemble = 1; ensemble < input.ensembles; ++ensemble)
      require(result.active_fields[3 * (result.active_macros * ensemble)] != result.active_fields[0],
              "ensemble batch was not kept independent");

   // A one-hot probe must land in exactly its documented component/basis/
   // ensemble batch. This catches a basis-slow or ensemble-interleaved pack.
   for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
      for(unsigned int basis = 0; basis < input.basis; ++basis)
         for(unsigned int component = 0; component < 3; ++component) {
            auto probe = input;
            probe.active_moments.assign(input.active_moments.size(), 0.0);
            const std::size_t active_index = component + 3 * (basis + result.active_macros * ensemble);
            probe.active_moments[active_index] = 1.0;
            const auto probe_result = luna_open_fft_test::run(probe);
            const std::size_t expected_batch = component + 3 * (basis + input.basis * ensemble);
            for(std::size_t field_batch = 0; field_batch < result.field_batches; ++field_batch)
               for(std::size_t q = 0; q < result.fft_cells; ++q) {
                  const double expected = q == 0 && field_batch == expected_batch ? 1.0 : 0.0;
                  require(probe_result.packed_moments[q + result.fft_cells * field_batch] == expected,
                          "component/basis/ensemble batch is not documented order");
               }
         }

   // Physical axial pair: a unit x impulse at cell 0 produces K_xx=2 at cell
   // 1 before the physical prefactor.  This catches source/target reversal and
   // verifies that the same kernel is used for field and energy.
   luna_open_fft_test::Input impulse = input;
   impulse.ensembles = 1;
   impulse.basis = 1;
   impulse.basis_offsets = {{0.0, 0.0, 0.0}};
   impulse.active_moments.assign(6, 0.0);
   impulse.active_moments[0] = 1.0;
   const auto impulse_result = luna_open_fft_test::run(impulse);
   require(std::abs(impulse_result.active_fields[3] - 2.0) < 1.0e-12,
           "same-kernel axial convolution differs before physical prefactor");
   require(std::abs(impulse_result.active_fields[0]) < 1.0e-12,
           "delta source wrapped to the opposite face");
   require(std::abs(impulse_result.dimensionless_energy[0] + 1.0) < 1.0e-12,
           "dimensionless energy is not -1/2 sum(M dot B)");

   // With both moments present, the centered finite difference of energy with
   // respect to m(cell 0,x) is -B(cell 0,x), before any physical prefactor.
   luna_open_fft_test::Input pair = impulse;
   pair.active_moments[3] = 1.0;
   const double step = 1.0e-7;
   const auto pair_base = luna_open_fft_test::run(pair);
   auto pair_plus = pair;
   auto pair_minus = pair;
   pair_plus.active_moments[0] += step;
   pair_minus.active_moments[0] -= step;
   const double derivative = (luna_open_fft_test::run(pair_plus).dimensionless_energy[0] -
                              luna_open_fft_test::run(pair_minus).dimensionless_energy[0]) / (2.0 * step);
   require(std::abs(derivative + pair_base.active_fields[0]) < 1.0e-9,
           "energy finite-difference derivative disagrees with field");
}

}  // namespace

int main() {
   runLayoutAndBatchTests();
   return 0;
}
