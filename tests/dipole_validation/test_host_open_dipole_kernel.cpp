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

double energy(const std::vector<double>& moments, const std::vector<double>& fields) {
   require(moments.size() == fields.size(), "energy field/moment size mismatch");
   double result = 0.0;
   for(std::size_t index = 0; index < moments.size(); ++index)
      result -= 0.5 * moments[index] * fields[index];
   return result;
}

std::array<double, 3> finePosition(
   const DipoleOpenGeometry& geometry,
   const std::array<std::size_t, 3>& fine_cell, unsigned int basis) {
   std::array<double, 3> result = geometry.basis_offsets[basis];
   for(unsigned int component = 0; component < 3; ++component)
      for(unsigned int axis = 0; axis < 3; ++axis)
         result[component] += static_cast<double>(fine_cell[axis]) *
                              geometry.primitive_vectors[3 * axis + component];
   return result;
}

// Independent test projection of Luna's finite block-one operator.  It forms
// the complete point-pair restriction P^T K_point P for one representative
// target/source block at every signed coarse displacement.  It deliberately
// does not call or inspect the production builder while producing expected
// values.
std::vector<double> explicitPointOperatorProjection(const DipoleOpenGeometry& geometry) {
   const std::size_t fft_cells =
      geometry.fft_grid[0] * geometry.fft_grid[1] * geometry.fft_grid[2];
   const std::size_t batches =
      9 * static_cast<std::size_t>(geometry.basis) * geometry.basis;
   const std::size_t population =
      geometry.block_shape[0] * geometry.block_shape[1] * geometry.block_shape[2];
   const double normalization =
      1.0 / (static_cast<double>(population) * static_cast<double>(population));
   std::vector<double> result(fft_cells * batches, 0.0);

   for(long d3 = -static_cast<long>(geometry.active_grid[2] - 1);
       d3 <= static_cast<long>(geometry.active_grid[2] - 1); ++d3)
      for(long d2 = -static_cast<long>(geometry.active_grid[1] - 1);
          d2 <= static_cast<long>(geometry.active_grid[1] - 1); ++d2)
         for(long d1 = -static_cast<long>(geometry.active_grid[0] - 1);
             d1 <= static_cast<long>(geometry.active_grid[0] - 1); ++d1) {
            const std::array<long, 3> displacement{{d1, d2, d3}};
            std::array<std::size_t, 3> target_block{}, source_block{}, q{};
            for(unsigned int axis = 0; axis < 3; ++axis) {
               if(displacement[axis] >= 0)
                  target_block[axis] = static_cast<std::size_t>(displacement[axis]);
               else
                  source_block[axis] = static_cast<std::size_t>(-displacement[axis]);
               q[axis] = displacement[axis] >= 0
                  ? static_cast<std::size_t>(displacement[axis])
                  : geometry.fft_grid[axis] -
                    static_cast<std::size_t>(-displacement[axis]);
            }
            const std::size_t qcell = cell(q[0], q[1], q[2], geometry.fft_grid);
            for(unsigned int target = 0; target < geometry.basis; ++target)
               for(unsigned int source = 0; source < geometry.basis; ++source)
                  for(std::size_t u3 = 0; u3 < geometry.block_shape[2]; ++u3)
                     for(std::size_t u2 = 0; u2 < geometry.block_shape[1]; ++u2)
                        for(std::size_t u1 = 0; u1 < geometry.block_shape[0]; ++u1)
                           for(std::size_t v3 = 0; v3 < geometry.block_shape[2]; ++v3)
                              for(std::size_t v2 = 0; v2 < geometry.block_shape[1]; ++v2)
                                 for(std::size_t v1 = 0; v1 < geometry.block_shape[0]; ++v1) {
                                    const std::array<std::size_t, 3> u{{u1, u2, u3}};
                                    const std::array<std::size_t, 3> v{{v1, v2, v3}};
                                    if(displacement == std::array<long, 3>{{0, 0, 0}} &&
                                       target == source && u == v) continue;
                                    std::array<std::size_t, 3> target_fine{}, source_fine{};
                                    for(unsigned int axis = 0; axis < 3; ++axis) {
                                       target_fine[axis] =
                                          target_block[axis] * geometry.block_shape[axis] + u[axis];
                                       source_fine[axis] =
                                          source_block[axis] * geometry.block_shape[axis] + v[axis];
                                    }
                                    const auto target_position =
                                       finePosition(geometry, target_fine, target);
                                    const auto source_position =
                                       finePosition(geometry, source_fine, source);
                                    const std::array<double, 3> r{{
                                       target_position[0] - source_position[0],
                                       target_position[1] - source_position[1],
                                       target_position[2] - source_position[2]}};
                                    const double r2 =
                                       r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                                    require(r2 > 0.0,
                                            "explicit projected point operator contains overlap");
                                    const double inverse_r3 = 1.0 / (r2 * std::sqrt(r2));
                                    const double inverse_r5 = inverse_r3 / r2;
                                    for(unsigned int column = 0; column < 3; ++column)
                                       for(unsigned int row = 0; row < 3; ++row)
                                          result[qcell + fft_cells *
                                             kernelBatch(row, column, target, source,
                                                         geometry.basis)] += normalization *
                                             (3.0 * r[row] * r[column] * inverse_r5 -
                                              (row == column ? inverse_r3 : 0.0));
                                 }
         }
   return result;
}

double maximumDifference(const std::vector<double>& left,
                         const std::vector<double>& right) {
   require(left.size() == right.size(), "kernel comparison shape mismatch");
   double maximum = 0.0;
   for(std::size_t index = 0; index < left.size(); ++index)
      maximum = std::max(maximum, std::abs(left[index] - right[index]));
   return maximum;
}

DipoleOpenGeometry skewTwoBasisGeometry() {
   DipoleOpenGeometry geometry{};
   geometry.atomistic_grid = {2, 2, 2};
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
   geometry.atomistic_grid = {2, 1, 1};
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

struct ProjectionSummary {
   double noncubic_error = 0.0;
   double skew_error = 0.0;
   double uniform_field_error = 0.0;
   double uniform_energy_error = 0.0;
};

ProjectionSummary runCoarseProjectionSuite() {
   DipoleOpenGeometry noncubic{};
   noncubic.atomistic_grid = {4, 2, 1};
   noncubic.active_grid = {2, 2, 1};
   noncubic.fft_grid = {5, 4, 1};
   noncubic.block_shape = {2, 1, 1};
   noncubic.primitive_vectors = {1.1, 0.0, 0.0,
                                 0.0, 0.8, 0.0,
                                 0.0, 0.0, 1.4};
   noncubic.basis = 1;
   noncubic.basis_offsets = {{0.0, 0.0, 0.0}};
   const auto noncubic_built = buildOpenDipoleDisplacementKernel(noncubic);
   const double noncubic_error = maximumDifference(
      noncubic_built.kernel, explicitPointOperatorProjection(noncubic));
   require(noncubic_error < 3.0e-14,
           "non-cubic coarse tensor differs from explicit block-one projection");

   // One isolated two-point block has two ordered distinct pairs.  With
   // M=m0+m1, Kxx=(2+2)/2^2=1 and Kyy=Kzz=(-1-1)/2^2=-1/2.
   DipoleOpenGeometry diagonal = noncubic;
   diagonal.atomistic_grid = {2, 1, 1};
   diagonal.active_grid = {1, 1, 1};
   diagonal.fft_grid = {1, 1, 1};
   diagonal.primitive_vectors = {1.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0,
                                 0.0, 0.0, 1.0};
   const auto diagonal_built = buildOpenDipoleDisplacementKernel(diagonal);
   require(diagonal_built.kernel[kernelBatch(0, 0, 0, 0, 1)] == 1.0 &&
           diagonal_built.kernel[kernelBatch(1, 1, 0, 0, 1)] == -0.5 &&
           diagonal_built.kernel[kernelBatch(2, 2, 0, 0, 1)] == -0.5,
           "finite intra-block point pairs did not produce the projected diagonal");
   require(diagonal_built.diagnostics.max_point_self_abs == 0.0 &&
           diagonal_built.diagnostics.max_projected_diagonal_abs == 1.0,
           "coarse diagnostics conflate point self with the finite projected diagonal");

   DipoleOpenGeometry skew{};
   skew.atomistic_grid = {2, 4, 2};
   skew.active_grid = {2, 2, 1};
   skew.fft_grid = {4, 5, 1};
   skew.block_shape = {1, 2, 2};
   skew.primitive_vectors = {1.0, 0.17, 0.04,
                             0.13, 1.2, 0.21,
                             0.08, 0.11, 0.9};
   skew.basis = 2;
   skew.basis_offsets = {{0.0, 0.0, 0.0}, {0.24, 0.16, 0.12}};
   const auto skew_built = buildOpenDipoleDisplacementKernel(skew);
   const double skew_error =
      maximumDifference(skew_built.kernel, explicitPointOperatorProjection(skew));
   require(skew_error < 2.0e-13,
           "skew two-basis coarse tensor differs from explicit block-one projection");
   require(skew_built.diagnostics.block_population == 4 &&
           skew_built.diagnostics.atomistic_cells == 16 &&
           skew_built.diagnostics.max_projected_diagonal_abs > 0.0,
           "coarse shape/population/diagonal diagnostics are incomplete");

   // Verify P^T K P at operator level: lift arbitrary block moments uniformly
   // to fine atoms, average the fine conjugate fields back with P^T, and
   // compare both fields and energies with the coarse convolution.
   const std::size_t coarse_cells = noncubic_built.diagnostics.active_cells;
   std::vector<double> macro_moments(3 * coarse_cells, 0.0);
   for(std::size_t index = 0; index < macro_moments.size(); ++index)
      macro_moments[index] = -0.37 + 0.11 * static_cast<double>(index);
   const std::size_t fine_cells =
      noncubic.atomistic_grid[0] * noncubic.atomistic_grid[1] *
      noncubic.atomistic_grid[2];
   std::vector<double> fine_moments(3 * fine_cells, 0.0);
   for(std::size_t fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
      const auto xyz = coordinates(fine_cell, noncubic.atomistic_grid);
      const std::size_t macro = cell(
         xyz[0] / noncubic.block_shape[0],
         xyz[1] / noncubic.block_shape[1],
         xyz[2] / noncubic.block_shape[2], noncubic.active_grid);
      for(unsigned int component = 0; component < 3; ++component)
         fine_moments[activeIndex(component, 0, fine_cell, 1, fine_cells)] =
            macro_moments[activeIndex(component, 0, macro, 1, coarse_cells)] / 2.0;
   }
   DipoleOpenGeometry fine_geometry = noncubic;
   fine_geometry.active_grid = fine_geometry.atomistic_grid;
   fine_geometry.fft_grid = {7, 3, 1};
   fine_geometry.block_shape = {1, 1, 1};
   const auto fine_fields = directFiniteFields(fine_geometry, fine_moments, 1);
   std::vector<double> averaged_fine_fields(macro_moments.size(), 0.0);
   for(std::size_t fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
      const auto xyz = coordinates(fine_cell, noncubic.atomistic_grid);
      const std::size_t macro = cell(
         xyz[0] / noncubic.block_shape[0],
         xyz[1] / noncubic.block_shape[1],
         xyz[2] / noncubic.block_shape[2], noncubic.active_grid);
      for(unsigned int component = 0; component < 3; ++component)
         averaged_fine_fields[activeIndex(component, 0, macro, 1, coarse_cells)] +=
            fine_fields[activeIndex(component, 0, fine_cell, 1, fine_cells)] / 2.0;
   }
   const auto coarse_fields =
      convolveBuiltKernel(noncubic, noncubic_built, macro_moments, 1);
   const double uniform_field_error =
      maximumDifference(coarse_fields, averaged_fine_fields);
   const double uniform_energy_error =
      std::abs(energy(fine_moments, fine_fields) -
               energy(macro_moments, coarse_fields));
   require(uniform_field_error < 3.0e-14 && uniform_energy_error < 3.0e-14,
           "uniform block moments do not preserve projected finite field/energy");

   // At block one the projection is the original basis-resolved point
   // operator: its diagonal is exact zero and the complete impulse matrix is
   // already checked above against Luna's direct finite sum.
   DipoleOpenGeometry block_one = skew;
   block_one.atomistic_grid = block_one.active_grid;
   block_one.block_shape = {1, 1, 1};
   const auto block_one_built = buildOpenDipoleDisplacementKernel(block_one);
   require(block_one_built.diagnostics.block_population == 1 &&
           block_one_built.diagnostics.max_projected_diagonal_abs == 0.0 &&
           maximumDifference(block_one_built.kernel,
                             explicitPointOperatorProjection(block_one)) < 2.0e-13,
           "block-one projection did not recover the finite point operator");

   return {noncubic_error, skew_error, uniform_field_error, uniform_energy_error};
}

struct ConvergenceSummary {
   double block_four_field = 0.0;
   double block_two_field = 0.0;
   double block_one_field = 0.0;
   double block_four_energy = 0.0;
   double block_two_energy = 0.0;
   double block_one_energy = 0.0;
};

ConvergenceSummary runNonuniformConvergenceSuite() {
   DipoleOpenGeometry fine{};
   fine.atomistic_grid = {4, 2, 1};
   fine.active_grid = fine.atomistic_grid;
   fine.fft_grid = {7, 3, 1};
   fine.primitive_vectors = {1.0, 0.11, 0.03,
                             0.08, 1.15, 0.07,
                             0.02, 0.13, 0.9};
   fine.basis = 1;
   fine.basis_offsets = {{0.0, 0.0, 0.0}};
   const std::size_t fine_cells = 8;
   std::vector<double> moments(3 * fine_cells, 0.0);
   for(std::size_t fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
      const auto xyz = coordinates(fine_cell, fine.atomistic_grid);
      moments[activeIndex(0, 0, fine_cell, 1, fine_cells)] =
         0.23 + 0.09 * xyz[0] - 0.04 * xyz[1];
      moments[activeIndex(1, 0, fine_cell, 1, fine_cells)] =
         -0.31 + 0.06 * xyz[0] + 0.08 * xyz[1];
      moments[activeIndex(2, 0, fine_cell, 1, fine_cells)] =
         0.19 - 0.05 * xyz[0] + 0.03 * xyz[1];
   }
   const auto fine_fields = directFiniteFields(fine, moments, 1);
   const double fine_energy = energy(moments, fine_fields);

   struct Error { double field = 0.0; double energy = 0.0; };
   const auto approximate = [&](std::size_t width) {
      DipoleOpenGeometry coarse = fine;
      coarse.block_shape = {width, 1, 1};
      coarse.active_grid = {fine.atomistic_grid[0] / width,
                            fine.atomistic_grid[1], fine.atomistic_grid[2]};
      coarse.fft_grid = {2 * coarse.active_grid[0] - 1,
                         2 * coarse.active_grid[1] - 1,
                         2 * coarse.active_grid[2] - 1};
      const auto built = buildOpenDipoleDisplacementKernel(coarse);
      const std::size_t macro_cells = built.diagnostics.active_cells;
      std::vector<double> macro_moments(3 * macro_cells, 0.0);
      for(std::size_t fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
         const auto xyz = coordinates(fine_cell, fine.atomistic_grid);
         const std::size_t macro =
            cell(xyz[0] / width, xyz[1], xyz[2], coarse.active_grid);
         for(unsigned int component = 0; component < 3; ++component)
            macro_moments[activeIndex(component, 0, macro, 1, macro_cells)] +=
               moments[activeIndex(component, 0, fine_cell, 1, fine_cells)];
      }
      const auto macro_fields =
         convolveBuiltKernel(coarse, built, macro_moments, 1);
      Error error{};
      for(std::size_t fine_cell = 0; fine_cell < fine_cells; ++fine_cell) {
         const auto xyz = coordinates(fine_cell, fine.atomistic_grid);
         const std::size_t macro =
            cell(xyz[0] / width, xyz[1], xyz[2], coarse.active_grid);
         for(unsigned int component = 0; component < 3; ++component)
            error.field = std::max(
               error.field,
               std::abs(fine_fields[activeIndex(component, 0, fine_cell, 1, fine_cells)] -
                        macro_fields[activeIndex(component, 0, macro, 1, macro_cells)]));
      }
      error.energy =
         std::abs(fine_energy - energy(macro_moments, macro_fields)) / fine_cells;
      return error;
   };

   const Error block_four = approximate(4);
   const Error block_two = approximate(2);
   const Error block_one = approximate(1);
   require(block_four.field > block_two.field && block_two.field > 1.0e-8,
           "nonuniform field fixture does not improve as block width decreases");
   require(block_four.energy > block_two.energy && block_two.energy > 1.0e-10,
           "nonuniform energy fixture does not improve as block width decreases");
   require(block_one.field < 3.0e-14 && block_one.energy < 3.0e-14,
           "nonuniform block-one limit did not recover the finite atomistic operator");
   return {block_four.field, block_two.field, block_one.field,
           block_four.energy, block_two.energy, block_one.energy};
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

   DipoleOpenGeometry zero_block = valid;
   zero_block.block_shape = {0, 1, 1};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(zero_block); }, "block_shape[0]");

   DipoleOpenGeometry partial = valid;
   partial.atomistic_grid = {5, 2, 2};
   partial.block_shape = {2, 1, 1};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(partial); }, "partial edge blocks");

   DipoleOpenGeometry partial_y = valid;
   partial_y.atomistic_grid = {2, 5, 2};
   partial_y.block_shape = {1, 2, 1};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(partial_y); }, "partial edge blocks");

   DipoleOpenGeometry partial_z = valid;
   partial_z.atomistic_grid = {2, 2, 5};
   partial_z.block_shape = {1, 1, 2};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(partial_z); }, "partial edge blocks");

   DipoleOpenGeometry unequal_shape = valid;
   unequal_shape.atomistic_grid = {4, 2, 2};
   unequal_shape.block_shape = {2, 1, 1};
   unequal_shape.active_grid = {3, 2, 2};
   requireThrows<std::invalid_argument>(
      [&] { (void)buildOpenDipoleDisplacementKernel(unequal_shape); },
      "every active macro channel");

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
   overlap.atomistic_grid = {2, 1, 1};
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
   singleton.atomistic_grid = {1, 1, 1};
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
   const ProjectionSummary projection = runCoarseProjectionSuite();
   const ConvergenceSummary convergence = runNonuniformConvergenceSuite();
   runValidationTests();
   std::printf("open-host-builder direct_convolution_max=%.17g reciprocity_max=%.17g "
               "minimum_nonself_r2=%.17g phase_error=%.17g "
               "projection_noncubic=%.3e projection_skew=%.3e "
               "uniform_field=%.3e uniform_energy=%.3e\n",
               summary.maximum_error, summary.reciprocity_error,
               summary.minimum_nonself_r2, phase_error,
               projection.noncubic_error, projection.skew_error,
               projection.uniform_field_error, projection.uniform_energy_error);
   std::printf("open-nonuniform-convergence block4(field=%.3e energy/atom=%.3e) "
               "block2(field=%.3e energy/atom=%.3e) "
               "block1(field=%.3e energy/atom=%.3e)\n",
               convergence.block_four_field, convergence.block_four_energy,
               convergence.block_two_field, convergence.block_two_energy,
               convergence.block_one_field, convergence.block_one_energy);
   return 0;
}
