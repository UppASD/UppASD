#include "dipoleOpenKernel.hpp"

#include <array>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

std::size_t cell(std::size_t x, std::size_t y, std::size_t z,
                 const std::array<std::size_t, 3>& grid) {
   return x + grid[0] * (y + grid[1] * z);
}

std::size_t batch(unsigned int row, unsigned int column,
                  unsigned int target, unsigned int source,
                  unsigned int basis) {
   return row + 3 * (column + 3 * (target + basis * source));
}

std::size_t activeIndex(unsigned int component, unsigned int basis,
                        std::size_t active_cell, unsigned int channels,
                        std::size_t active_cells, unsigned int ensemble) {
   return component + 3 * (basis + channels * (active_cell + active_cells * ensemble));
}

std::array<std::size_t, 3> coordinates(std::size_t active_cell,
                                      const std::array<std::size_t, 3>& grid) {
   const std::size_t x = active_cell % grid[0];
   active_cell /= grid[0];
   const std::size_t y = active_cell % grid[1];
   return {x, y, active_cell / grid[1]};
}

template <typename Value>
void read(Value& value, const char* field) {
   if(!(std::cin >> value)) throw std::invalid_argument(std::string("failed to read ") + field);
}

} // namespace

int main() {
   try {
      DipoleOpenGeometry geometry{};
      unsigned int ensembles = 0;
      for(auto& value : geometry.active_grid) read(value, "active_grid");
      geometry.atomistic_grid = geometry.active_grid;
      for(auto& value : geometry.fft_grid) read(value, "fft_grid");
      for(auto& value : geometry.block_shape) read(value, "block_shape");
      read(geometry.basis, "basis");
      read(ensembles, "ensembles");
      for(auto& value : geometry.primitive_vectors) read(value, "primitive_vectors");
      geometry.basis_offsets.resize(geometry.basis);
      for(auto& offset : geometry.basis_offsets)
         for(auto& value : offset) read(value, "basis_offsets");

      const auto built = buildOpenDipoleDisplacementKernel(geometry);
      const std::size_t active_cells = built.diagnostics.active_cells;
      const std::size_t value_count =
         3 * active_cells * static_cast<std::size_t>(geometry.basis) * ensembles;
      std::vector<double> moments(value_count);
      for(auto& value : moments) read(value, "moments");
      std::vector<double> fields(value_count, 0.0);

      for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
         for(std::size_t target_cell = 0; target_cell < active_cells; ++target_cell) {
            const auto target_xyz = coordinates(target_cell, geometry.active_grid);
            for(std::size_t source_cell = 0; source_cell < active_cells; ++source_cell) {
               const auto source_xyz = coordinates(source_cell, geometry.active_grid);
               const std::array<std::size_t, 3> q{{
                  (target_xyz[0] + geometry.fft_grid[0] - source_xyz[0]) % geometry.fft_grid[0],
                  (target_xyz[1] + geometry.fft_grid[1] - source_xyz[1]) % geometry.fft_grid[1],
                  (target_xyz[2] + geometry.fft_grid[2] - source_xyz[2]) % geometry.fft_grid[2]}};
               const std::size_t qcell = cell(q[0], q[1], q[2], geometry.fft_grid);
               for(unsigned int target = 0; target < geometry.basis; ++target)
                  for(unsigned int source = 0; source < geometry.basis; ++source)
                     for(unsigned int row = 0; row < 3; ++row)
                        for(unsigned int column = 0; column < 3; ++column)
                           fields[activeIndex(row, target, target_cell, geometry.basis,
                                              active_cells, ensemble)] +=
                              built.kernel[qcell + built.diagnostics.fft_cells *
                                 batch(row, column, target, source, geometry.basis)] *
                              moments[activeIndex(column, source, source_cell, geometry.basis,
                                                  active_cells, ensemble)];
            }
         }

      const auto& diagnostics = built.diagnostics;
      std::cout << std::setprecision(17)
                << fields.size() << ' '
                << diagnostics.active_cells << ' '
                << diagnostics.fft_cells << ' '
                << diagnostics.kernel_batches << ' '
                << diagnostics.nonfinite_values << ' '
                << (diagnostics.all_finite ? 1 : 0) << ' '
                << diagnostics.minimum_nonself_r2 << ' '
                << diagnostics.max_reciprocity_error << ' '
                << diagnostics.max_point_self_abs << ' '
                << diagnostics.max_padding_gap_abs << '\n';
      for(const double value : fields) std::cout << value << '\n';
      return 0;
   } catch(const std::exception& error) {
      std::cerr << "open_kernel_host_driver: " << error.what() << '\n';
      return 1;
   }
}
