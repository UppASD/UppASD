#include "dipoleOpenKernel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
   if(!condition) throw std::runtime_error(message);
}

std::size_t cell(std::size_t x, std::size_t y, std::size_t z,
                 const std::array<std::size_t, 3>& grid) {
   return x + grid[0] * (y + grid[1] * z);
}

std::size_t kernelBatch(unsigned int row, unsigned int column,
                        unsigned int target, unsigned int source,
                        unsigned int basis) {
   return row + 3 * (column + 3 * (target + basis * source));
}

std::size_t activeIndex(unsigned int component, unsigned int basis,
                        std::size_t active_cell, unsigned int channels,
                        std::size_t active_cells, unsigned int ensemble = 0) {
   return component + 3 * (basis + channels * (active_cell + active_cells * ensemble));
}

std::array<std::size_t, 3> coordinates(std::size_t active_cell,
                                      const std::array<std::size_t, 3>& grid) {
   const std::size_t x = active_cell % grid[0];
   active_cell /= grid[0];
   const std::size_t y = active_cell % grid[1];
   return {x, y, active_cell / grid[1]};
}

std::vector<double> convolveBuiltKernel(const DipoleOpenGeometry& geometry,
                                        const DipoleOpenKernelResult& built,
                                        const std::vector<double>& moments,
                                        unsigned int ensembles) {
   const std::size_t active_cells = built.diagnostics.active_cells;
   const std::size_t fft_cells = built.diagnostics.fft_cells;
   std::vector<double> fields(moments.size(), 0.0);
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
                           built.kernel[qcell + fft_cells *
                              kernelBatch(row, column, target, source, geometry.basis)] *
                           moments[activeIndex(column, source, source_cell, geometry.basis,
                                               active_cells, ensemble)];
         }
      }
   return fields;
}

std::array<double, 3> position(const DipoleOpenGeometry& geometry,
                               std::size_t active_cell, unsigned int basis) {
   const auto xyz = coordinates(active_cell, geometry.active_grid);
   std::array<double, 3> result = geometry.basis_offsets[basis];
   for(unsigned int component = 0; component < 3; ++component)
      for(unsigned int axis = 0; axis < 3; ++axis)
         result[component] += static_cast<double>(xyz[axis]) *
                              geometry.primitive_vectors[3 * axis + component];
   return result;
}

std::vector<double> directFiniteFields(const DipoleOpenGeometry& geometry,
                                       const std::vector<double>& moments,
                                       unsigned int ensembles) {
   const std::size_t active_cells =
      geometry.active_grid[0] * geometry.active_grid[1] * geometry.active_grid[2];
   std::vector<double> fields(moments.size(), 0.0);
   for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
      for(std::size_t target_cell = 0; target_cell < active_cells; ++target_cell)
         for(unsigned int target = 0; target < geometry.basis; ++target) {
            const auto target_position = position(geometry, target_cell, target);
            for(std::size_t source_cell = 0; source_cell < active_cells; ++source_cell)
               for(unsigned int source = 0; source < geometry.basis; ++source) {
                  if(target_cell == source_cell && target == source) continue;
                  const auto source_position = position(geometry, source_cell, source);
                  const std::array<double, 3> r{{
                     target_position[0] - source_position[0],
                     target_position[1] - source_position[1],
                     target_position[2] - source_position[2]}};
                  const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                  require(r2 > 0.0, "direct test geometry contains a nonself overlap");
                  const double inverse_r3 = 1.0 / (r2 * std::sqrt(r2));
                  const double inverse_r5 = inverse_r3 / r2;
                  for(unsigned int row = 0; row < 3; ++row)
                     for(unsigned int column = 0; column < 3; ++column)
                        fields[activeIndex(row, target, target_cell, geometry.basis,
                                           active_cells, ensemble)] +=
                           (3.0 * r[row] * r[column] * inverse_r5 -
                            (row == column ? inverse_r3 : 0.0)) *
                           moments[activeIndex(column, source, source_cell, geometry.basis,
                                               active_cells, ensemble)];
               }
         }
   return fields;
}

DipoleOpenGeometry skewTwoBasisGeometry() {
   DipoleOpenGeometry geometry{};
   geometry.active_grid = {2, 2, 2};
   // Deliberately larger than the minimum (3,3,3) in every axis.
   geometry.fft_grid = {6, 7, 5};
   geometry.primitive_vectors = {1.0, 0.15, 0.05,
                                 0.10, 1.25, 0.20,
                                 0.05, 0.10, 0.85};
   geometry.basis = 2;
   geometry.basis_offsets = {{0.0, 0.0, 0.0}, {0.23, 0.17, 0.11}};
   return geometry;
}

struct ImpulseSummary {
   double maximum_error = 0.0;
   double reciprocity_error = 0.0;
   double minimum_nonself_r2 = 0.0;
};

ImpulseSummary runCompleteImpulseMatrix() {
   const DipoleOpenGeometry geometry = skewTwoBasisGeometry();
   const DipoleOpenKernelResult built = buildOpenDipoleDisplacementKernel(geometry);
   const std::size_t active_cells = built.diagnostics.active_cells;
   const std::size_t active_values = 3 * active_cells * geometry.basis;
   double maximum_error = 0.0;

   // Every source cell, source basis, and source component is made one-hot.
   // Comparing every target basis and target component exercises all nine
   // tensor components and every ordered basis pair.
   for(std::size_t source_cell = 0; source_cell < active_cells; ++source_cell)
      for(unsigned int source_basis = 0; source_basis < geometry.basis; ++source_basis)
         for(unsigned int source_component = 0; source_component < 3; ++source_component) {
            std::vector<double> moments(active_values, 0.0);
            moments[activeIndex(source_component, source_basis, source_cell,
                                geometry.basis, active_cells)] = 1.0;
            const auto actual = convolveBuiltKernel(geometry, built, moments, 1);
            const auto expected = directFiniteFields(geometry, moments, 1);
            for(std::size_t index = 0; index < actual.size(); ++index)
               maximum_error = std::max(maximum_error, std::abs(actual[index] - expected[index]));
         }

   require(built.diagnostics.active_cells == 8, "active_cells diagnostic is wrong");
   require(built.diagnostics.fft_cells == 210, "fft_cells diagnostic is wrong");
   require(built.diagnostics.kernel_batches == 36, "kernel_batches diagnostic is wrong");
   require(built.kernel.size() == 210 * 36, "padded kernel element count is wrong");
   require(built.diagnostics.all_finite && built.diagnostics.nonfinite_values == 0,
           "finite-value diagnostics failed");
   require(built.diagnostics.minimum_nonself_r2 > 0.0,
           "minimum nonself distance must be positive");
   require(built.diagnostics.max_reciprocity_error < 2.0e-12,
           "reciprocity diagnostic exceeds fp64 construction headroom");
   require(built.diagnostics.max_point_self_abs == 0.0,
           "point-self diagnostic must be exact");
   require(built.diagnostics.max_padding_gap_abs == 0.0,
           "oversized padding gap must remain exactly zero");
   require(maximum_error < 2.0e-12,
           "built-tensor direct convolution differs from the independent point sum");
   return {maximum_error, built.diagnostics.max_reciprocity_error,
           built.diagnostics.minimum_nonself_r2};
}

double runHandWorkedBasisPhase() {
   DipoleOpenGeometry geometry{};
   geometry.active_grid = {2, 1, 1};
   geometry.fft_grid = {4, 1, 1};
   geometry.primitive_vectors = {1.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0,
                                 0.0, 0.0, 1.0};
   geometry.basis = 2;
   geometry.basis_offsets = {{0.0, 0.0, 0.0}, {0.25, 0.0, 0.0}};
   const auto built = buildOpenDipoleDisplacementKernel(geometry);
   const std::size_t forward_batch = kernelBatch(0, 0, 0, 1, 2);
   const std::size_t reciprocal_batch = kernelBatch(0, 0, 1, 0, 2);
   const double wanted = 128.0 / 27.0;
   const double reversed_basis_phase = 128.0 / 125.0;
   const double forward = built.kernel[1 + 4 * forward_batch];
   const double reciprocal = built.kernel[3 + 4 * reciprocal_batch];
   require(std::abs(forward - wanted) < 2.0e-15,
           "target-minus-source basis phase is not 128/27");
   require(std::abs(reciprocal - wanted) < 2.0e-15,
           "reciprocal basis phase is not 128/27");
   require(std::abs(forward - reversed_basis_phase) > 1.0,
           "basis phase test does not distinguish reversed offsets");
   return std::max(std::abs(forward - wanted), std::abs(reciprocal - wanted));
}

template <typename Exception, typename Function>
void requireThrows(Function&& function, const std::string& field) {
   try {
      function();
   } catch(const Exception& error) {
      require(std::string(error.what()).find(field) != std::string::npos,
              "exception does not name field " + field + ": " + error.what());
      return;
   }
   throw std::runtime_error("expected exception was not thrown for " + field);
}

void runValidationTests() {
   const DipoleOpenGeometry valid = skewTwoBasisGeometry();

   DipoleOpenGeometry too_small = valid;
   too_small.fft_grid[0] = 2;
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(too_small); }, "fft_grid[0]");

   DipoleOpenGeometry coarse = valid;
   coarse.block_shape = {2, 1, 1};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(coarse); }, "block_shape[0]");

   DipoleOpenGeometry zero = valid;
   zero.active_grid[2] = 0;
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(zero); }, "active_grid[2]");

   DipoleOpenGeometry mismatch = valid;
   mismatch.basis_offsets.pop_back();
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(mismatch); }, "basis_offsets");

   DipoleOpenGeometry left_handed = valid;
   left_handed.primitive_vectors[0] = -1.0;
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(left_handed); }, "primitive_vectors determinant");

   DipoleOpenGeometry nonfinite = valid;
   nonfinite.basis_offsets[1][2] = std::numeric_limits<double>::infinity();
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(nonfinite); }, "basis_offsets[1][2]");

   DipoleOpenGeometry overlap{};
   overlap.active_grid = {2, 1, 1};
   overlap.fft_grid = {3, 1, 1};
   overlap.primitive_vectors = {1.0, 0.0, 0.0,
                                0.0, 1.0, 0.0,
                                0.0, 0.0, 1.0};
   overlap.basis = 2;
   overlap.basis_offsets = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(overlap); }, "overlap");

   DipoleOpenGeometry extent_overflow = valid;
   extent_overflow.active_grid[0] = std::numeric_limits<std::size_t>::max() / 2 + 2;
   requireThrows<std::overflow_error>(
      [&] { (void)buildOpenDipoleDisplacementKernel(extent_overflow); }, "active_grid");

   DipoleOpenGeometry singleton{};
   singleton.active_grid = {1, 1, 1};
   singleton.fft_grid = {1, 1, 1};
   singleton.primitive_vectors = {1.0, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0};
   singleton.basis = 1;
   singleton.basis_offsets = {{0.0, 0.0, 0.0}};
   const auto single = buildOpenDipoleDisplacementKernel(singleton);
   require(single.kernel == std::vector<double>(9, 0.0),
           "singleton kernel must be exact zero self");
   require(single.diagnostics.minimum_nonself_r2 == 0.0,
           "singleton has no nonself distance");
}

} // namespace

int main() {
   const ImpulseSummary summary = runCompleteImpulseMatrix();
   const double phase_error = runHandWorkedBasisPhase();
   runValidationTests();
   std::printf("open-host-builder direct_convolution_max=%.17g reciprocity_max=%.17g "
               "minimum_nonself_r2=%.17g phase_error=%.17g\n",
               summary.maximum_error, summary.reciprocity_error,
               summary.minimum_nonself_r2, phase_error);
   return 0;
}
