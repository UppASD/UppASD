// WP10.1 red tests.  This file is intentionally not part of the current CMake
// targets: it becomes buildable/linkable when Terra supplies
// open_fft_test_seam.hpp's luna_open_fft_test::run implementation.

#include "open_fft_test_seam.hpp"
#include "dipoleOpenKernel.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
   if(!condition) throw std::runtime_error(message);
}

std::size_t cell(std::size_t x, std::size_t y, std::size_t z,
                 const std::array<std::size_t, 3>& grid) {
   return x + grid[0] * (y + grid[1] * z);
}

std::size_t kernelBatch(unsigned int row, unsigned int column, unsigned int target,
                        unsigned int source, unsigned int basis) {
   return row + 3 * (column + 3 * (target + basis * source));
}

void uploadPointTestKernel(luna_open_fft_test::Input& input) {
   const std::size_t fft_cells = input.fft_grid[0] * input.fft_grid[1] * input.fft_grid[2];
   input.padded_real_kernel.assign(fft_cells * 9 * input.basis * input.basis, 0.0);
   for(std::size_t gz = 0; gz < input.active_grid[2]; ++gz)
      for(std::size_t gy = 0; gy < input.active_grid[1]; ++gy)
         for(std::size_t gx = 0; gx < input.active_grid[0]; ++gx)
            for(std::size_t hz = 0; hz < input.active_grid[2]; ++hz)
               for(std::size_t hy = 0; hy < input.active_grid[1]; ++hy)
                  for(std::size_t hx = 0; hx < input.active_grid[0]; ++hx)
                     for(unsigned int target = 0; target < input.basis; ++target)
                        for(unsigned int source = 0; source < input.basis; ++source) {
                           const long long dx = static_cast<long long>(gx) - static_cast<long long>(hx);
                           const long long dy = static_cast<long long>(gy) - static_cast<long long>(hy);
                           const long long dz = static_cast<long long>(gz) - static_cast<long long>(hz);
                           const std::size_t qx = dx >= 0 ? static_cast<std::size_t>(dx) : input.fft_grid[0] - static_cast<std::size_t>(-dx);
                           const std::size_t qy = dy >= 0 ? static_cast<std::size_t>(dy) : input.fft_grid[1] - static_cast<std::size_t>(-dy);
                           const std::size_t qz = dz >= 0 ? static_cast<std::size_t>(dz) : input.fft_grid[2] - static_cast<std::size_t>(-dz);
                           const std::size_t q = cell(qx, qy, qz, input.fft_grid);
                           double r[3]{};
                           for(unsigned int component = 0; component < 3; ++component)
                              r[component] = static_cast<double>(dx) * input.primitive_vectors[component] +
                                 static_cast<double>(dy) * input.primitive_vectors[3 + component] +
                                 static_cast<double>(dz) * input.primitive_vectors[6 + component] +
                                 input.basis_offsets[target][component] - input.basis_offsets[source][component];
                           const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                           if(r2 == 0.0) continue;
                           const double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
                           const double inv_r5 = inv_r3 / r2;
                           for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column)
                              input.padded_real_kernel[q + fft_cells * kernelBatch(row, column, target, source, input.basis)] =
                                 3.0 * r[row] * r[column] * inv_r5 - (row == column ? inv_r3 : 0.0);
                        }
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
   uploadPointTestKernel(input);
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
   require(result.persistent_bytes == result.persistent_inventory_bytes &&
           result.construction_bytes == result.construction_inventory_bytes &&
           result.total_bytes == result.persistent_bytes + result.workspace_bytes + result.construction_bytes,
           "OPEN_FFT memory estimate does not equal its allocation inventory");
   requireZeroPadding(result, input);
   for(unsigned int ensemble = 1; ensemble < input.ensembles; ++ensemble)
      require(result.active_fields[3 * (result.active_macros * ensemble)] != result.active_fields[0],
              "ensemble batch was not kept independent");

   // The seam also accepts an arbitrary padded real tensor; this local delta
   // is deliberately not a physical dipole tensor and proves that OPEN_FFT
   // remains a test upload API rather than a production kernel builder.
   auto uploaded = input;
   uploaded.padded_real_kernel.assign(result.fft_cells * 9 * input.basis * input.basis, 0.0);
   for(unsigned int channel = 0; channel < input.basis; ++channel)
      for(unsigned int component = 0; component < 3; ++component) {
         const std::size_t batch = component + 3 * (component + 3 * (channel + input.basis * channel));
         uploaded.padded_real_kernel[result.fft_cells * batch] = 1.0;
      }
   const auto uploaded_result = luna_open_fft_test::run(uploaded);
   for(std::size_t index = 0; index < input.active_moments.size(); ++index)
      require(std::abs(uploaded_result.active_fields[index] - input.active_moments[index]) < 1.0e-12,
              "arbitrary padded real-kernel upload did not preserve active coordinates");

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
   uploadPointTestKernel(impulse);
   const auto impulse_result = luna_open_fft_test::run(impulse);
   require(std::abs(impulse_result.active_fields[3] - 2.0) < 1.0e-12,
           "same-kernel axial convolution differs before physical prefactor");
   require(std::abs(impulse_result.active_fields[0]) < 1.0e-12,
           "delta source wrapped to the opposite face");
   // The sole source has no point self term, so its active energy is zero
   // even though its field at the other active cell is K_xx=2.
   require(std::abs(impulse_result.dimensionless_energy[0]) < 1.0e-12,
           "point-self exclusion leaked into dimensionless energy");

   // With both moments present, the centered finite difference of energy with
   // respect to m(cell 0,x) is -B(cell 0,x), before any physical prefactor.
   luna_open_fft_test::Input pair = impulse;
   pair.active_moments[3] = 1.0;
   const double step = 1.0e-7;
   const auto pair_base = luna_open_fft_test::run(pair);
   require(std::abs(pair_base.dimensionless_energy[0] + 2.0) < 1.0e-12,
           "dimensionless pair energy is not -1/2 sum(M dot B)");
   auto pair_plus = pair;
   auto pair_minus = pair;
   pair_plus.active_moments[0] += step;
   pair_minus.active_moments[0] -= step;
   const double derivative = (luna_open_fft_test::run(pair_plus).dimensionless_energy[0] -
                              luna_open_fft_test::run(pair_minus).dimensionless_energy[0]) / (2.0 * step);
   std::printf("open-layout finite-difference derivative=%.17g negative_field=%.17g\n",
               derivative, -pair_base.active_fields[0]);
   std::fflush(stdout);
   require(std::abs(derivative + pair_base.active_fields[0]) < 5.0e-9,
           "energy finite-difference derivative disagrees with field");
}

std::size_t macroIndex(unsigned int component, unsigned int basis, std::size_t cell,
                       unsigned int channels, std::size_t active_cells, unsigned int ensemble) {
   return component + 3 * (basis + channels * (cell + active_cells * ensemble));
}

std::vector<double> directFiniteFields(const luna_open_fft_test::Input& input) {
   const std::size_t active_cells = input.active_grid[0] * input.active_grid[1] * input.active_grid[2];
   std::vector<double> fields(input.active_moments.size(), 0.0);
   for(std::size_t tz = 0; tz < input.active_grid[2]; ++tz)
      for(std::size_t ty = 0; ty < input.active_grid[1]; ++ty)
         for(std::size_t tx = 0; tx < input.active_grid[0]; ++tx) {
            const std::size_t target_cell = cell(tx, ty, tz, input.active_grid);
            for(std::size_t sz = 0; sz < input.active_grid[2]; ++sz)
               for(std::size_t sy = 0; sy < input.active_grid[1]; ++sy)
                  for(std::size_t sx = 0; sx < input.active_grid[0]; ++sx) {
                     const std::size_t source_cell = cell(sx, sy, sz, input.active_grid);
                     const double d[3] = {static_cast<double>(tx) - static_cast<double>(sx),
                                          static_cast<double>(ty) - static_cast<double>(sy),
                                          static_cast<double>(tz) - static_cast<double>(sz)};
                     for(unsigned int target = 0; target < input.basis; ++target)
                        for(unsigned int source = 0; source < input.basis; ++source) {
                           double r[3]{};
                           for(unsigned int component = 0; component < 3; ++component)
                              r[component] = d[0] * input.primitive_vectors[component] +
                                 d[1] * input.primitive_vectors[3 + component] +
                                 d[2] * input.primitive_vectors[6 + component] +
                                 input.basis_offsets[target][component] - input.basis_offsets[source][component];
                           const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
                           if(r2 == 0.0) continue;
                           const double inv_r3 = 1.0 / (r2 * std::sqrt(r2));
                           const double inv_r5 = inv_r3 / r2;
                           for(unsigned int row = 0; row < 3; ++row)
                              for(unsigned int column = 0; column < 3; ++column)
                                 for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
                                    fields[macroIndex(row, target, target_cell, input.basis, active_cells, ensemble)] +=
                                       (3.0 * r[row] * r[column] * inv_r5 - (row == column ? inv_r3 : 0.0)) *
                                       input.active_moments[macroIndex(column, source, source_cell, input.basis,
                                                                      active_cells, ensemble)];
                        }
                  }
         }
   return fields;
}

std::vector<double> directFiniteEnergies(const luna_open_fft_test::Input& input,
                                         const std::vector<double>& fields) {
   const std::size_t active_cells = input.active_grid[0] * input.active_grid[1] * input.active_grid[2];
   std::vector<double> energies(input.ensembles, 0.0);
   for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
      for(std::size_t cell_index = 0; cell_index < active_cells; ++cell_index)
         for(unsigned int basis = 0; basis < input.basis; ++basis)
            for(unsigned int component = 0; component < 3; ++component) {
               const std::size_t index = macroIndex(component, basis, cell_index, input.basis,
                                                    active_cells, ensemble);
               energies[ensemble] -= 0.5 * input.active_moments[index] * fields[index];
            }
   return energies;
}

void compareProductionResult(const luna_open_fft_test::Input& input,
                             const luna_open_fft_test::Result& result,
                             double& maximum_field_error, double& maximum_energy_error) {
   const auto expected_fields = directFiniteFields(input);
   const auto expected_energies = directFiniteEnergies(input, expected_fields);
   require(result.active_fields.size() == expected_fields.size(), "production field shape is wrong");
   require(result.dimensionless_energy.size() == expected_energies.size(), "production energy shape is wrong");
   require(result.scattered_fields.size() == expected_fields.size(), "production scatter field shape is wrong");
   require(result.energy_per_atom.size() == expected_energies.size(), "production energy shape is wrong");
   for(std::size_t index = 0; index < expected_fields.size(); ++index)
      maximum_field_error = std::max(maximum_field_error,
                                     std::abs(result.active_fields[index] - expected_fields[index]));
   for(std::size_t index = 0; index < expected_fields.size(); ++index)
      maximum_field_error = std::max(maximum_field_error,
                                     std::abs(result.scattered_fields[index] - expected_fields[index]));
   for(std::size_t ensemble = 0; ensemble < expected_energies.size(); ++ensemble)
      maximum_energy_error = std::max(maximum_energy_error,
                                      std::abs(result.dimensionless_energy[ensemble] - expected_energies[ensemble]));
   const double atoms_inverse = 1.0 / static_cast<double>(input.active_grid[0] * input.active_grid[1] *
                                                           input.active_grid[2] * input.basis);
   for(std::size_t ensemble = 0; ensemble < expected_energies.size(); ++ensemble)
      maximum_energy_error = std::max(maximum_energy_error,
                                      std::abs(result.energy_per_atom[ensemble] -
                                               expected_energies[ensemble] * atoms_inverse));
}

void runOpenProductionE2E() {
   struct Case {
      std::array<std::size_t, 3> grid;
      std::array<std::size_t, 3> padding;
      std::array<double, 9> primitive_vectors;
      unsigned int ensembles;
   };
   const std::array<Case, 2> cases{{
      {{{2, 1, 1}}, {{3, 1, 1}},
       {{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}}, 1},
      {{{2, 3, 1}}, {{5, 7, 3}},
       {{1.0, 0.15, 0.05, 0.10, 1.25, 0.20, 0.05, 0.10, 0.85}}, 4}
   }};
   double maximum_field_error = 0.0;
   double maximum_energy_error = 0.0;
   for(const auto& test_case : cases) {
      luna_open_fft_test::Input input{};
      input.active_grid = test_case.grid;
      input.fft_grid = test_case.padding;
      input.basis = 2;
      input.ensembles = test_case.ensembles;
      input.primitive_vectors = test_case.primitive_vectors;
      input.basis_offsets = {{0.0, 0.0, 0.0}, {0.23, 0.17, 0.11}};
      const std::size_t active_cells = input.active_grid[0] * input.active_grid[1] * input.active_grid[2];
      input.active_moments.assign(3 * active_cells * input.basis * input.ensembles, 0.0);

      // A nonuniform basis-resolved state checks both M=1 and M=4 production
      // batches against Luna's independent finite point sum.
      for(std::size_t index = 0; index < input.active_moments.size(); ++index)
         input.active_moments[index] = -0.31 + 0.019 * static_cast<double>(index);
      const auto nonuniform = luna_open_fft_test::runProduction(input);
      compareProductionResult(input, nonuniform, maximum_field_error, maximum_energy_error);

      DipoleOpenGeometry geometry{};
      geometry.atomistic_grid = input.active_grid;
      geometry.active_grid = input.active_grid;
      geometry.fft_grid = input.fft_grid;
      geometry.primitive_vectors = input.primitive_vectors;
      geometry.basis = input.basis;
      geometry.basis_offsets = input.basis_offsets;
      const auto built = buildOpenDipoleDisplacementKernel(geometry);
      require(built.diagnostics.max_reciprocity_error < 2.0e-12,
              "basis-resolved production kernel violates reciprocity");

      // Relabeling basis channels must not change any physical field or energy.
      auto swapped = input;
      std::swap(swapped.basis_offsets[0], swapped.basis_offsets[1]);
      for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
         for(std::size_t cell_index = 0; cell_index < active_cells; ++cell_index)
            for(unsigned int component = 0; component < 3; ++component) {
               const std::size_t first = macroIndex(component, 0, cell_index, input.basis, active_cells, ensemble);
               const std::size_t second = macroIndex(component, 1, cell_index, input.basis, active_cells, ensemble);
               swapped.active_moments[first] = input.active_moments[second];
               swapped.active_moments[second] = input.active_moments[first];
            }
      const auto permuted = luna_open_fft_test::runProduction(swapped);
      for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
         for(std::size_t cell_index = 0; cell_index < active_cells; ++cell_index)
            for(unsigned int component = 0; component < 3; ++component) {
               const std::size_t first = macroIndex(component, 0, cell_index, input.basis, active_cells, ensemble);
               const std::size_t second = macroIndex(component, 1, cell_index, input.basis, active_cells, ensemble);
               maximum_field_error = std::max(maximum_field_error,
                  std::abs(permuted.active_fields[first] - nonuniform.active_fields[second]));
               maximum_field_error = std::max(maximum_field_error,
                  std::abs(permuted.active_fields[second] - nonuniform.active_fields[first]));
            }
      for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble)
         maximum_energy_error = std::max(maximum_energy_error,
            std::abs(permuted.dimensionless_energy[ensemble] - nonuniform.dimensionless_energy[ensemble]));

      // Every source basis/component is impulsed in every active cell and
      // every ensemble.  This fixes target-minus-source phase by a direct
      // finite comparison rather than relying on the tensor's same-basis symmetry.
      for(std::size_t z = 0; z < input.active_grid[2]; ++z)
         for(std::size_t y = 0; y < input.active_grid[1]; ++y)
            for(std::size_t x = 0; x < input.active_grid[0]; ++x) {
               const std::size_t source_cell = cell(x, y, z, input.active_grid);
               for(unsigned int basis = 0; basis < input.basis; ++basis)
                  for(unsigned int component = 0; component < 3; ++component)
                     for(unsigned int ensemble = 0; ensemble < input.ensembles; ++ensemble) {
                        std::fill(input.active_moments.begin(), input.active_moments.end(), 0.0);
                        input.active_moments[macroIndex(component, basis, source_cell, input.basis, active_cells, ensemble)] =
                           1.0 + 0.125 * ensemble;
                        compareProductionResult(input, luna_open_fft_test::runProduction(input),
                                                maximum_field_error, maximum_energy_error);
                     }
            }

      // On the shared production field/energy operator, dE/dM is -B.
      input = swapped;
      const auto base = luna_open_fft_test::runProduction(input);
      // This Hamiltonian is quadratic in M, so a centered derivative has no
      // truncation error.  Keep the probe well above fp64 cancellation in the
      // atomically reduced per-atom energy rather than loosening its oracle
      // gate for a numerically under-resolved 1e-7 difference.
      constexpr double step = 1.0e-5;
      const std::size_t differentiated = macroIndex(1, 1, active_cells - 1, input.basis, active_cells, 0);
      auto plus = input;
      auto minus = input;
      plus.active_moments[differentiated] += step;
      minus.active_moments[differentiated] -= step;
      const auto plus_result = luna_open_fft_test::runProduction(plus);
      const auto minus_result = luna_open_fft_test::runProduction(minus);
      const double derivative = (plus_result.dimensionless_energy[0] -
                                 minus_result.dimensionless_energy[0]) / (2.0 * step);
      maximum_energy_error = std::max(maximum_energy_error,
                                      std::abs(derivative + base.active_fields[differentiated]));
      const double per_atom_derivative = (plus_result.energy_per_atom[0] -
                                          minus_result.energy_per_atom[0]) / (2.0 * step);
      maximum_energy_error = std::max(maximum_energy_error,
         std::abs(per_atom_derivative + base.scattered_fields[differentiated] /
                  static_cast<double>(active_cells * input.basis)));
   }
   std::printf("open-production-e2e impulses field_max=%.17g energy_max=%.17g\n",
               maximum_field_error, maximum_energy_error);
   // The field/energy comparisons are fp64 FFT checks; the aggregate energy
   // figure also includes the centered finite-difference derivative.
   require(maximum_field_error < 2.0e-11 && maximum_energy_error < 5.0e-9,
           "OPEN_FFT basis-resolved production E2E differs from the finite direct convolution");
}

}  // namespace

int main() {
   runLayoutAndBatchTests();
   runOpenProductionE2E();
   return 0;
}
