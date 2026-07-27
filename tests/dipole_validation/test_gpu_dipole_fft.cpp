#include "dipoleEwaldKernel.hpp"
#include "gpuDipoleConvolution.hpp"
#include "gpu_ewald_goldens_v1.hpp"
#include "gpu_ewald_multibasis_goldens_v1.hpp"
#include "tensor.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Grid = std::array<std::size_t, 3>;

std::size_t cells(const Grid& grid) {
   return grid[0] * grid[1] * grid[2];
}

std::size_t kernelBatch(unsigned int row, unsigned int column, unsigned int target, unsigned int source,
                        unsigned int basis) {
   return row + 3 * (column + 3 * (target + basis * source));
}

std::size_t macroIndex(unsigned int component, unsigned int basis_channel, std::size_t cell,
                       unsigned int basis, std::size_t grid_cells, unsigned int ensemble) {
   return component + 3 * (basis_channel + basis * (cell + grid_cells * ensemble));
}

std::array<double, 9> reciprocal(const std::array<double, 9>& h, double volume) {
   const double scale = 2.0 * std::acos(-1.0) / volume;
   return {scale * (h[4] * h[8] - h[7] * h[5]), scale * (h[5] * h[6] - h[3] * h[8]),
           scale * (h[3] * h[7] - h[4] * h[6]), scale * (h[7] * h[2] - h[8] * h[1]),
           scale * (h[8] * h[0] - h[6] * h[2]), scale * (h[6] * h[1] - h[7] * h[0]),
           scale * (h[1] * h[5] - h[2] * h[4]), scale * (h[2] * h[3] - h[0] * h[5]),
           scale * (h[0] * h[4] - h[1] * h[3])};
}

double volume(const std::array<double, 9>& h) {
   return h[0] * (h[4] * h[8] - h[7] * h[5]) - h[3] * (h[1] * h[8] - h[7] * h[2]) +
          h[6] * (h[1] * h[5] - h[4] * h[2]);
}

GpuDipoleConvolutionDescriptor descriptor(const Grid& grid, unsigned int basis, unsigned int ensembles,
                                          const std::array<double, 9>& h,
                                          const std::vector<std::array<double, 3>>& offsets = {}) {
   static_assert(sizeof(real) == sizeof(double), "GPU FFT validation is fp64 only");
   // The descriptor stores primitive vectors and derives H by multiplying by
   // atomistic grid extents.  The test uses the same grid for both notions.
   static thread_local std::array<real, 3> c1, c2, c3;
   c1 = {h[0] / grid[0], h[1] / grid[0], h[2] / grid[0]};
   c2 = {h[3] / grid[1], h[4] / grid[1], h[5] / grid[1]};
   c3 = {h[6] / grid[2], h[7] / grid[2], h[8] / grid[2]};
   GpuDipoleConvolutionDescriptor result{};
   result.atomistic_grid = {grid[0], grid[1], grid[2]};
   result.macro_grid = result.atomistic_grid;
   result.basis = basis;
   result.ensembles = ensembles;
   result.boundary = GpuDipoleBoundaryMode::Periodic3D;
   result.discretization = GpuDipoleDiscretization::MacrospinGrid;
   result.c1 = c1.data(); result.c2 = c2.data(); result.c3 = c3.data();
   result.basis_offsets = offsets;
   if(result.basis_offsets.empty()) result.basis_offsets.assign(basis, {0.0, 0.0, 0.0});
   // No macro-centre data are needed for this pure convolution test.
   result.alat = 1.0;
   result.tolerance = 1.0e-10;
   result.field_prefactor = 1.0e-7 * 9.274009994e-24;
   if(!result.valid()) throw std::runtime_error("GPU delta-test descriptor is invalid");
   return result;
}

std::vector<double> directConvolution(const std::vector<double>& kernel, const std::vector<double>& moments,
                                      const Grid& grid, unsigned int basis, unsigned int ensembles) {
   const std::size_t ncell = cells(grid);
   std::vector<double> fields(moments.size(), 0.0);
   for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
      for(std::size_t tz = 0; tz < grid[2]; ++tz) for(std::size_t ty = 0; ty < grid[1]; ++ty)
         for(std::size_t tx = 0; tx < grid[0]; ++tx) {
            const std::size_t target_cell = tx + grid[0] * (ty + grid[1] * tz);
            for(std::size_t sz = 0; sz < grid[2]; ++sz) for(std::size_t sy = 0; sy < grid[1]; ++sy)
               for(std::size_t sx = 0; sx < grid[0]; ++sx) {
                  const std::size_t source_cell = sx + grid[0] * (sy + grid[1] * sz);
                  const std::size_t dx = (tx + grid[0] - sx) % grid[0];
                  const std::size_t dy = (ty + grid[1] - sy) % grid[1];
                  const std::size_t dz = (tz + grid[2] - sz) % grid[2];
                  const std::size_t displacement = dx + grid[0] * (dy + grid[1] * dz);
                  for(unsigned int target = 0; target < basis; ++target)
                     for(unsigned int source = 0; source < basis; ++source)
                        for(unsigned int row = 0; row < 3; ++row)
                           for(unsigned int column = 0; column < 3; ++column)
                              fields[macroIndex(row, target, target_cell, basis, ncell, ensemble)] +=
                                 kernel[displacement + ncell * kernelBatch(row, column, target, source, basis)] *
                                 moments[macroIndex(column, source, source_cell, basis, ncell, ensemble)];
               }
         }
   return fields;
}

std::vector<double> deltaKernel(const Grid& grid, unsigned int basis) {
   const std::size_t ncell = cells(grid);
   std::vector<double> kernel(ncell * 9 * basis * basis);
   for(std::size_t displacement = 0; displacement < ncell; ++displacement)
      for(unsigned int target = 0; target < basis; ++target)
         for(unsigned int source = 0; source < basis; ++source)
            for(unsigned int row = 0; row < 3; ++row)
               for(unsigned int column = 0; column < 3; ++column)
                  // Deliberately non-symmetric in displacement, Cartesian, and basis indices.
                  kernel[displacement + ncell * kernelBatch(row, column, target, source, basis)] =
                     0.03125 * (1.0 + displacement + 3.0 * row + 5.0 * column + 7.0 * target + 11.0 * source);
   return kernel;
}

struct Errors { double field = 0.0; double energy = 0.0; };

struct OperatorResult {
   std::vector<real> fields;
   std::vector<real> energies;
   GpuDipoleSpectrumDiagnostics spectrum;
};

double reciprocalIdentityError(const std::array<double, 9>& h, const std::array<double, 9>& b) {
   double maximum = 0.0;
   for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
      double value = 0.0;
      for(unsigned int component = 0; component < 3; ++component)
         value += h[component + 3 * row] * b[component + 3 * column];
      maximum = std::max(maximum, std::abs(value - (row == column ? 2.0 * std::acos(-1.0) : 0.0)));
   }
   return maximum;
}

std::array<double, 9> reciprocalWithLegacyTypo(const std::array<double, 9>& h, double cell_volume) {
   auto result = reciprocal(h, cell_volume);
   const double scale = 2.0 * std::acos(-1.0) / cell_volume;
   // Historical test-helper error: b2.x used h[5], which masks itself if
   // h[1] and h[5] are both zero.
   result[3] = scale * (h[7] * h[2] - h[8] * h[5]);
   return result;
}

void runReciprocalHelperRegression() {
   const std::array<double, 9> skew{{6.0, 0.45, 0.25, 0.8, 3.0, 0.35, 0.3, 0.2, 3.4}};
   const double cell_volume = volume(skew);
   const double corrected = reciprocalIdentityError(skew, reciprocal(skew, cell_volume));
   const double legacy = reciprocalIdentityError(skew, reciprocalWithLegacyTypo(skew, cell_volume));
   std::printf("reciprocal-helper corrected_residual=%.17g legacy_typo_residual=%.17g\n", corrected, legacy);
   if(corrected > 2.0e-14 || legacy < 1.0e-3)
      throw std::runtime_error("reciprocal helper skew regression did not distinguish h[1] from h[5]");
}

OperatorResult evaluateOperator(const Grid& grid, unsigned int basis, unsigned int ensembles,
                                const std::array<double, 9>& h, const std::vector<double>& kernel,
                                const std::vector<double>& moments, bool check_not_ready = false,
                                bool nonphysical_plumbing_kernel = false) {
   const std::size_t ncell = cells(grid);
   GpuTensor<real, 3> device_moments;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("GPU stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(ncell * basis), static_cast<long int>(ensembles));
      if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("GPU moment upload failed");
      GpuDipoleConvolution solver;
      if(!solver.initiate(descriptor(grid, basis, ensembles, h), stream)) throw std::runtime_error("GPU FFT plan initiation failed");
      if(check_not_ready) {
         bool guarded = false;
         try { solver.evaluate(device_moments); } catch(const std::runtime_error&) { guarded = true; }
         if(!guarded) throw std::runtime_error("evaluate did not reject a missing complete kernel");
      }
      solver.uploadCompleteKernelForTesting(kernel, !nonphysical_plumbing_kernel);
      if(!solver.kernelReady()) throw std::runtime_error("kernel_ready was not set after normalized kernel construction");
      if(solver.diagnosticConstructionStorageAllocatedForTesting())
         throw std::runtime_error("transient real kernel storage survived kernel construction");
      solver.evaluate(device_moments);
      OperatorResult result{solver.diagnosticFieldsForTesting(), solver.diagnosticEnergiesForTesting(),
                            solver.diagnosticSpectrumForTesting()};
      solver.release();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("GPU stream destruction failed");
      return result;
   } catch(...) {
      if(!device_moments.empty()) device_moments.Free();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

// Exercise the production Builder-A entry point as well as the standalone
// upload path.  The expected values used by its callers remain independent
// Python-oracle fixtures.
OperatorResult evaluateConstructedOperator(const Grid& grid, unsigned int basis, unsigned int ensembles,
                                           const std::array<double, 9>& h,
                                           const std::vector<std::array<double, 3>>& offsets,
                                           const std::vector<double>& moments) {
   const std::size_t ncell = cells(grid);
   GpuTensor<real, 3> device_moments;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("GPU constructed-kernel stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(ncell * basis), static_cast<long int>(ensembles));
      if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("GPU constructed-kernel moment upload failed");
      GpuDipoleConvolution solver;
      if(!solver.initiate(descriptor(grid, basis, ensembles, h, offsets), stream))
         throw std::runtime_error("GPU constructed-kernel plan initiation failed");
      solver.buildPeriodicKernel();
      if(!solver.kernelReady()) throw std::runtime_error("GPU constructed-kernel did not become ready");
      solver.evaluate(device_moments);
      OperatorResult result{solver.diagnosticFieldsForTesting(), solver.diagnosticEnergiesForTesting(),
                            solver.diagnosticSpectrumForTesting()};
      solver.release();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("GPU constructed-kernel stream destruction failed");
      return result;
   } catch(...) {
      if(!device_moments.empty()) device_moments.Free();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

std::array<OperatorResult, 2> evaluateConsecutive(const Grid& grid, unsigned int basis, unsigned int ensembles,
                                                   const std::array<double, 9>& h, const std::vector<double>& kernel,
                                                   const std::vector<double>& first, const std::vector<double>& second) {
   if(first.size() != second.size()) throw std::runtime_error("changed-moment sequence has incompatible sizes");
   const std::size_t ncell = cells(grid);
   GpuTensor<real, 3> device_moments;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("GPU stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(ncell * basis), static_cast<long int>(ensembles));
      GpuDipoleConvolution solver;
      if(!solver.initiate(descriptor(grid, basis, ensembles, h), stream)) throw std::runtime_error("GPU FFT plan initiation failed");
      solver.uploadCompleteKernelForTesting(kernel);
      std::array<OperatorResult, 2> result{};
      const std::array<const std::vector<double>*, 2> input{{&first, &second}};
      for(unsigned int evaluation = 0; evaluation < input.size(); ++evaluation) {
         if(GPU_MEMCPY(device_moments.data(), input[evaluation]->data(), input[evaluation]->size() * sizeof(real),
                       GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS) {
            throw std::runtime_error("GPU changed-moment upload failed");
         }
         solver.evaluate(device_moments);
         result[evaluation] = {solver.diagnosticFieldsForTesting(), solver.diagnosticEnergiesForTesting(),
                               solver.diagnosticSpectrumForTesting()};
      }
      solver.release();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("GPU stream destruction failed");
      return result;
   } catch(...) {
      if(!device_moments.empty()) device_moments.Free();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

Errors runOperator(const std::string& name, const Grid& grid, unsigned int basis, unsigned int ensembles,
                   const std::array<double, 9>& h, const std::vector<double>& kernel,
                   const std::vector<double>& moments, bool check_not_ready = false,
                   bool nonphysical_plumbing_kernel = false) {
   const std::size_t ncell = cells(grid);
   const auto result = evaluateOperator(grid, basis, ensembles, h, kernel, moments, check_not_ready,
                                        nonphysical_plumbing_kernel);
   const auto expected = directConvolution(kernel, moments, grid, basis, ensembles);
   Errors errors{};
   for(std::size_t i = 0; i < result.fields.size(); ++i) errors.field = std::max(errors.field, std::abs(result.fields[i] - expected[i]));
   for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble) {
      double expected_energy = 0.0;
      for(std::size_t cell = 0; cell < ncell; ++cell) for(unsigned int channel = 0; channel < basis; ++channel)
         for(unsigned int component = 0; component < 3; ++component) {
            const auto index = macroIndex(component, channel, cell, basis, ncell, ensemble);
            expected_energy -= 0.5 * moments[index] * expected[index];
         }
      errors.energy = std::max(errors.energy, std::abs(result.energies[ensemble] - expected_energy));
   }
   std::printf("plumbing %s max_field_error=%.17g max_energy_error=%.17g\n", name.c_str(), errors.field, errors.energy);
   return errors;
}

void require(const Errors& errors, const std::string& name) {
   constexpr double tolerance = 5.0e-11;
   if(errors.field > tolerance || errors.energy > tolerance)
      throw std::runtime_error(name + " exceeds fp64 GPU convolution tolerance");
}

void runDeltaSuite() {
   const std::array<Grid, 4> grids{{{1, 1, 1}, {2, 1, 1}, {2, 3, 1}, {3, 2, 2}}};
   const std::array<double, 9> h{{6.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0, 8.0}};
   for(const Grid& grid : grids) for(const unsigned int basis : {1U, 2U}) for(const unsigned int ensembles : {1U, 4U}) {
      const std::size_t ncell = cells(grid);
      const auto kernel = deltaKernel(grid, basis);
      // Every source basis and Cartesian component is an impulse.  For M=4
      // each ensemble has a different amplitude, exercising batch isolation.
      for(unsigned int source = 0; source < basis; ++source) for(unsigned int component = 0; component < 3; ++component) {
         std::vector<double> moments(3 * ncell * basis * ensembles, 0.0);
         for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
            moments[macroIndex(component, source, ncell - 1, basis, ncell, ensemble)] = 1.0 + 0.25 * ensemble;
         const std::string name = "delta grid=" + std::to_string(grid[0]) + "x" + std::to_string(grid[1]) + "x" +
                                  std::to_string(grid[2]) + " basis=" + std::to_string(basis) + " ensembles=" +
                                  std::to_string(ensembles) + " source=" + std::to_string(source) + "/" + std::to_string(component);
         require(runOperator(name, grid, basis, ensembles, h, kernel, moments, grid == grids.front() && basis == 1 &&
                             ensembles == 1 && source == 0 && component == 0, true), name);
      }
   }
}

void runPeriodicCase(const std::string& name, const Grid& grid, unsigned int basis, unsigned int ensembles,
                     const std::array<double, 9>& h, const std::vector<std::array<double, 3>>& offsets) {
   DipolePeriodicGeometry geometry{};
   geometry.H = h; geometry.volume = volume(h); geometry.Brec = reciprocal(h, geometry.volume);
   geometry.grid = grid; geometry.basis = basis; geometry.basis_offsets = offsets;
   const auto complete = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10});
   const std::size_t ncell = cells(grid);
   std::vector<double> moments(3 * ncell * basis * ensembles);
   for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
      for(std::size_t cell = 0; cell < ncell; ++cell) for(unsigned int channel = 0; channel < basis; ++channel)
         for(unsigned int component = 0; component < 3; ++component)
            moments[macroIndex(component, channel, cell, basis, ncell, ensemble)] =
               0.05 * (1.0 + component + 2.0 * channel + 3.0 * cell + 5.0 * ensemble);
   require(runOperator("periodic " + name, grid, basis, ensembles, h, complete.kernel, moments), "periodic " + name);
}

void runPeriodicSuite() {
   runPeriodicCase("cubic-1x1x1", {1, 1, 1}, 1, 1,
                   {3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0}, {{0.0, 0.0, 0.0}});
   runPeriodicCase("axial-2x1x1", {2, 1, 1}, 1, 1,
                   {6.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0}, {{0.0, 0.0, 0.0}});
   runPeriodicCase("noncubic-2x3x1", {2, 3, 1}, 1, 1,
                   {6.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 3.0}, {{0.0, 0.0, 0.0}});
   runPeriodicCase("skew-2x1x1", {2, 1, 1}, 1, 1,
                   {6.0, 0.0, 0.0, 0.8, 3.0, 0.0, 0.3, 0.2, 3.4}, {{0.0, 0.0, 0.0}});
   runPeriodicCase("two-basis-3x2x2-m4", {3, 2, 2}, 2, 4,
                   {9.0, 0.0, 0.0, 0.7, 6.0, 0.0, 0.3, 0.2, 6.4}, {{0.0, 0.0, 0.0}, {0.37, 0.23, 0.41}});
}

DipolePeriodicGeometry periodicGeometry(const Grid& grid, const std::array<double, 9>& h) {
   DipolePeriodicGeometry geometry{};
   geometry.H = h;
   geometry.volume = volume(h);
   geometry.Brec = reciprocal(h, geometry.volume);
   geometry.grid = grid;
   geometry.basis = 1;
   geometry.basis_offsets = {{0.0, 0.0, 0.0}};
   return geometry;
}

double maxDifference(const std::vector<real>& actual, const double* expected, std::size_t count) {
   double result = 0.0;
   for(std::size_t index = 0; index < count; ++index) result = std::max(result, std::abs(actual[index] - expected[index]));
   return result;
}

void requireSpectrum(const GpuDipoleSpectrumDiagnostics& spectrum, const std::string& name) {
   constexpr double tolerance = 1.0e-12;
   if(spectrum.max_reciprocity_error > tolerance || spectrum.max_conjugacy_error > tolerance ||
      spectrum.max_hermitian_error > tolerance) {
      throw std::runtime_error(name + " violates complete-kernel spectral validation");
   }
}

void runIndependentOracleSuite() {
   constexpr double tolerance = 3.0e-10;
   double maximum_field = 0.0, maximum_energy = 0.0, maximum_reciprocity = 0.0, maximum_hermitian = 0.0;
   for(const auto& oracle : gpu_ewald_golden_v1::cases) {
      const Grid grid{{oracle.grid[0], oracle.grid[1], oracle.grid[2]}};
      const auto geometry = periodicGeometry(grid, oracle.H);
      const auto complete = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10});
      const std::size_t count = 3 * cells(grid) * oracle.ensembles;
      const std::vector<double> moments(oracle.moments, oracle.moments + count);
      const auto actual = evaluateOperator(grid, 1, static_cast<unsigned int>(oracle.ensembles), oracle.H,
                                           complete.kernel, moments);
      const double field_error = maxDifference(actual.fields, oracle.fields, count);
      double energy_error = 0.0;
      for(std::size_t ensemble = 0; ensemble < oracle.ensembles; ++ensemble)
         energy_error = std::max(energy_error, std::abs(actual.energies[ensemble] - oracle.energies[ensemble]));
      requireSpectrum(actual.spectrum, oracle.name);
      if(field_error > tolerance || energy_error > tolerance)
         throw std::runtime_error(std::string("independent oracle mismatch for ") + oracle.name);
      maximum_field = std::max(maximum_field, field_error);
      maximum_energy = std::max(maximum_energy, energy_error);
      maximum_reciprocity = std::max(maximum_reciprocity, actual.spectrum.max_reciprocity_error);
      maximum_hermitian = std::max(maximum_hermitian, actual.spectrum.max_hermitian_error);
      std::printf("oracle %s max_field_error=%.17g max_energy_error=%.17g reciprocity=%.17g conjugacy=%.17g hermitian=%.17g\n",
                  oracle.name, field_error, energy_error, actual.spectrum.max_reciprocity_error,
                  actual.spectrum.max_conjugacy_error, actual.spectrum.max_hermitian_error);
   }
   std::printf("oracle-matrix max_field_error=%.17g max_energy_error=%.17g max_reciprocity=%.17g max_hermitian=%.17g\n",
               maximum_field, maximum_energy, maximum_reciprocity, maximum_hermitian);
}

void runMultiBasisOracleSuite() {
   constexpr double tolerance = 3.0e-10;
   double maximum_upload_field = 0.0, maximum_upload_energy = 0.0;
   double maximum_production_field = 0.0, maximum_production_energy = 0.0;
   for(const auto& oracle : gpu_ewald_multibasis_golden_v1::cases) {
      const Grid grid{{oracle.grid[0], oracle.grid[1], oracle.grid[2]}};
      const unsigned int basis = static_cast<unsigned int>(oracle.basis);
      const unsigned int ensembles = static_cast<unsigned int>(oracle.ensembles);
      std::vector<std::array<double, 3>> offsets(basis);
      for(unsigned int channel = 0; channel < basis; ++channel)
         for(unsigned int component = 0; component < 3; ++component)
            offsets[channel][component] = oracle.offsets[component + 3 * channel];
      DipolePeriodicGeometry geometry{};
      geometry.H = oracle.H;
      geometry.volume = volume(oracle.H);
      geometry.Brec = reciprocal(oracle.H, geometry.volume);
      geometry.grid = grid;
      geometry.basis = basis;
      geometry.basis_offsets = offsets;
      const auto complete = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10});
      const std::size_t count = 3 * cells(grid) * basis * ensembles;
      const std::vector<double> moments(oracle.moments, oracle.moments + count);
      const auto uploaded = evaluateOperator(grid, basis, ensembles, oracle.H, complete.kernel, moments);
      const auto constructed = evaluateConstructedOperator(grid, basis, ensembles, oracle.H, offsets, moments);
      const auto assess = [&](const OperatorResult& result, double& max_field, double& max_energy, const char* route) {
         const double field_error = maxDifference(result.fields, oracle.fields, count);
         double energy_error = 0.0;
         for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
            energy_error = std::max(energy_error, std::abs(result.energies[ensemble] - oracle.energies[ensemble]));
         requireSpectrum(result.spectrum, oracle.name);
         if(field_error > tolerance || energy_error > tolerance)
            throw std::runtime_error(std::string("multi-basis independent oracle mismatch for ") + oracle.name + " via " + route);
         max_field = std::max(max_field, field_error);
         max_energy = std::max(max_energy, energy_error);
         std::printf("multi-basis oracle %s route=%s max_field_error=%.17g max_energy_error=%.17g\n",
                     oracle.name, route, field_error, energy_error);
      };
      assess(uploaded, maximum_upload_field, maximum_upload_energy, "upload");
      assess(constructed, maximum_production_field, maximum_production_energy, "production-builder");

      const auto low_alpha = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10, 0.70});
      const auto high_alpha = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10, 1.05});
      const auto low_result = evaluateOperator(grid, basis, ensembles, oracle.H, low_alpha.kernel, moments);
      const auto high_result = evaluateOperator(grid, basis, ensembles, oracle.H, high_alpha.kernel, moments);
      const double alpha_field_error = maxDifference(low_result.fields, high_result.fields.data(), count);
      double alpha_energy_error = 0.0;
      for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
         alpha_energy_error = std::max(alpha_energy_error,
                                       std::abs(low_result.energies[ensemble] - high_result.energies[ensemble]));
      if(alpha_field_error > tolerance || alpha_energy_error > tolerance)
         throw std::runtime_error(std::string("multi-basis explicit-alpha GPU invariance failed for ") + oracle.name);
      std::printf("multi-basis property %s alpha_field=%.17g alpha_energy=%.17g\n",
                  oracle.name, alpha_field_error, alpha_energy_error);

      // Swapping basis labels and moment channels changes labels only.  It is
      // particularly sensitive to target/source block ordering and phase sign.
      if(basis == 2) {
         auto swapped_offsets = offsets;
         std::swap(swapped_offsets[0], swapped_offsets[1]);
         std::vector<double> swapped_moments(moments.size());
         for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
            for(std::size_t cell = 0; cell < cells(grid); ++cell)
               for(unsigned int component = 0; component < 3; ++component) {
                  swapped_moments[macroIndex(component, 0, cell, basis, cells(grid), ensemble)] =
                     moments[macroIndex(component, 1, cell, basis, cells(grid), ensemble)];
                  swapped_moments[macroIndex(component, 1, cell, basis, cells(grid), ensemble)] =
                     moments[macroIndex(component, 0, cell, basis, cells(grid), ensemble)];
               }
         const auto swapped = evaluateConstructedOperator(grid, basis, ensembles, oracle.H, swapped_offsets, swapped_moments);
         double permutation_error = 0.0;
         for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
            for(std::size_t cell = 0; cell < cells(grid); ++cell)
               for(unsigned int component = 0; component < 3; ++component) {
                  permutation_error = std::max(permutation_error, std::abs(
                     swapped.fields[macroIndex(component, 0, cell, basis, cells(grid), ensemble)] -
                     constructed.fields[macroIndex(component, 1, cell, basis, cells(grid), ensemble)]));
                  permutation_error = std::max(permutation_error, std::abs(
                     swapped.fields[macroIndex(component, 1, cell, basis, cells(grid), ensemble)] -
                     constructed.fields[macroIndex(component, 0, cell, basis, cells(grid), ensemble)]));
               }
         for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
            permutation_error = std::max(permutation_error, std::abs(swapped.energies[ensemble] - constructed.energies[ensemble]));
         if(permutation_error > tolerance)
            throw std::runtime_error(std::string("multi-basis permutation property failed for ") + oracle.name);
         std::printf("multi-basis property %s basis_permutation=%.17g\n", oracle.name, permutation_error);
      }
   }
   std::printf("multi-basis oracle matrix upload_field=%.17g upload_energy=%.17g production_field=%.17g production_energy=%.17g\n",
               maximum_upload_field, maximum_upload_energy, maximum_production_field, maximum_production_energy);
}

void runMultiBasisProductionSliceSuite() {
   // This is the WP6 production boundary: basis-fast macro channels are
   // deliberately scattered through a non-identity atom-to-macro map and the
   // energy reduction visits every (cell,basis) channel exactly once.
   constexpr double mu_b = 9.274009994e-24;
   constexpr double tolerance = 3.0e-10;
   const auto& oracle = gpu_ewald_multibasis_golden_v1::cases[1];
   const Grid grid{{oracle.grid[0], oracle.grid[1], oracle.grid[2]}};
   const unsigned int basis = static_cast<unsigned int>(oracle.basis);
   const unsigned int ensembles = static_cast<unsigned int>(oracle.ensembles);
   const std::size_t ncell = cells(grid), atoms = ncell * basis;
   const std::size_t count = 3 * atoms * ensembles;
   std::vector<std::array<double, 3>> offsets(basis);
   for(unsigned int channel = 0; channel < basis; ++channel)
      for(unsigned int component = 0; component < 3; ++component)
         offsets[channel][component] = oracle.offsets[component + 3 * channel];
   const std::vector<double> moments(oracle.moments, oracle.moments + count);
   // atom -> one-based macro, intentionally not atom order.  Macro ordering
   // is basis + NA*cell, matching create_pme_macrocell_layout().
   const std::vector<unsigned int> map{{3, 1, 4, 2}};
   GpuTensor<real, 3> device_moments, beff, eneff;
   GpuTensor<unsigned int, 1> device_map;
   GpuTensor<real, 2> energy;
   GPU_STREAM_T stream{};
   if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("multi-basis production stream creation failed");
   try {
      device_moments.Allocate(3L, static_cast<long int>(atoms), static_cast<long int>(ensembles));
      beff.Allocate(3L, static_cast<long int>(atoms), static_cast<long int>(ensembles));
      eneff.Allocate(3L, static_cast<long int>(atoms), static_cast<long int>(ensembles));
      device_map.Allocate(static_cast<long int>(atoms));
      energy.Allocate(static_cast<long int>(ensembles), 7L);
      if(GPU_MEMCPY(device_moments.data(), moments.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
         GPU_MEMCPY(device_map.data(), map.data(), atoms * sizeof(unsigned int), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
         throw std::runtime_error("multi-basis production input upload failed");
      beff.zeros_async(stream); eneff.zeros_async(stream); energy.zeros_async(stream);
      GpuDipoleConvolution solver;
      if(!solver.initiate(descriptor(grid, basis, ensembles, oracle.H, offsets), stream))
         throw std::runtime_error("multi-basis production plan initiation failed");
      solver.buildPeriodicKernel();
      solver.evaluate(device_moments);
      solver.addFieldsToAtoms(beff, eneff, device_map.data(), atoms);
      solver.accumulateEnergy(energy, atoms);
      if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("multi-basis production synchronization failed");
      std::vector<real> actual_beff(count), actual_eneff(count), actual_energy(7 * ensembles);
      if(GPU_MEMCPY(actual_beff.data(), beff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
         GPU_MEMCPY(actual_eneff.data(), eneff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
         GPU_MEMCPY(actual_energy.data(), energy.data(), actual_energy.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
         throw std::runtime_error("multi-basis production output download failed");
      const double prefactor = 1.0e-7 * mu_b;
      double field_error = 0.0, energy_error = 0.0;
      for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble)
         for(std::size_t atom = 0; atom < atoms; ++atom) {
            const std::size_t macro = map[atom] - 1;
            const std::size_t cell = macro / basis;
            const unsigned int channel = static_cast<unsigned int>(macro % basis);
            for(unsigned int component = 0; component < 3; ++component) {
               const std::size_t atom_index = component + 3 * (atom + atoms * ensemble);
               const std::size_t field_index = macroIndex(component, channel, cell, basis, ncell, ensemble);
               const double expected = prefactor * oracle.fields[field_index];
               field_error = std::max(field_error, std::abs(actual_beff[atom_index] - expected));
               field_error = std::max(field_error, std::abs(actual_eneff[atom_index] - expected));
            }
         }
      for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble) {
         const double expected = prefactor * oracle.energies[ensemble] / static_cast<double>(atoms);
         energy_error = std::max(energy_error, std::abs(actual_energy[ensemble + ensembles * 5] - expected));
         energy_error = std::max(energy_error, std::abs(actual_energy[ensemble + ensembles * 6] - expected));
      }
      if(field_error > tolerance * prefactor || energy_error > tolerance * prefactor)
         throw std::runtime_error("multi-basis production scatter/energy exceeds the WP6 fp64 budget");
      std::printf("multi-basis production L10-skew field_error=%.17g energy_error=%.17g\n", field_error, energy_error);
      solver.release();
      energy.Free(); device_map.Free(); eneff.Free(); beff.Free(); device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("multi-basis production stream destruction failed");
   } catch(...) {
      if(!energy.empty()) energy.Free();
      if(!device_map.empty()) device_map.Free();
      if(!eneff.empty()) eneff.Free();
      if(!beff.empty()) beff.Free();
      if(!device_moments.empty()) device_moments.Free();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
}

void runAlphaAndPropertySuite() {
   constexpr double tolerance = 3.0e-10;
   const auto& cubic = gpu_ewald_golden_v1::cases[0];
   const Grid cubic_grid{{cubic.grid[0], cubic.grid[1], cubic.grid[2]}};
   const auto cubic_geometry = periodicGeometry(cubic_grid, cubic.H);
   const std::vector<double> cubic_moments(cubic.moments, cubic.moments + 3 * cells(cubic_grid));
   const auto low = buildPeriodicEwaldDisplacementKernel(cubic_geometry, {1.0e-10, 0.50});
   const auto high = buildPeriodicEwaldDisplacementKernel(cubic_geometry, {1.0e-10, 1.00});
   const auto low_result = evaluateOperator(cubic_grid, 1, 1, cubic.H, low.kernel, cubic_moments);
   const auto high_result = evaluateOperator(cubic_grid, 1, 1, cubic.H, high.kernel, cubic_moments);
   const double alpha_field = maxDifference(low_result.fields, high_result.fields.data(), low_result.fields.size());
   const double alpha_energy = std::abs(low_result.energies[0] - high_result.energies[0]);
   if(alpha_field > tolerance || alpha_energy > tolerance)
      throw std::runtime_error("explicit-alpha GPU invariance exceeds the WP4b fp64 budget");
   std::printf("explicit-alpha low=%.17g high=%.17g field_error=%.17g energy_error=%.17g\n",
               low.diagnostics.selected.alpha, high.diagnostics.selected.alpha, alpha_field, alpha_energy);

   // Translation and sign-flip use an independently checked NA=1 fixture.
   const auto& axial = gpu_ewald_golden_v1::cases[1];
   const Grid grid{{axial.grid[0], axial.grid[1], axial.grid[2]}};
   const auto geometry = periodicGeometry(grid, axial.H);
   const auto kernel = buildPeriodicEwaldDisplacementKernel(geometry, {1.0e-10}).kernel;
   const std::vector<double> moments(axial.moments, axial.moments + 3 * cells(grid));
   const auto baseline = evaluateOperator(grid, 1, 1, axial.H, kernel, moments);
   std::vector<double> shifted(moments.size());
   for(std::size_t cell = 0; cell < cells(grid); ++cell)
      for(unsigned int component = 0; component < 3; ++component)
         shifted[macroIndex(component, 0, (cell + 1) % cells(grid), 1, cells(grid), 0)] =
            moments[macroIndex(component, 0, cell, 1, cells(grid), 0)];
   const auto translated = evaluateOperator(grid, 1, 1, axial.H, kernel, shifted);
   double translation_error = 0.0;
   for(std::size_t cell = 0; cell < cells(grid); ++cell) for(unsigned int component = 0; component < 3; ++component)
      translation_error = std::max(translation_error, std::abs(
         translated.fields[macroIndex(component, 0, (cell + 1) % cells(grid), 1, cells(grid), 0)] -
         baseline.fields[macroIndex(component, 0, cell, 1, cells(grid), 0)]));
   std::vector<double> negated(moments);
   for(double& value : negated) value = -value;
   const std::vector<double> zero(moments.size(), 0.0);
   const auto sequence = evaluateConsecutive(grid, 1, 1, axial.H, kernel, moments, negated);
   const auto zero_result = evaluateOperator(grid, 1, 1, axial.H, kernel, zero);
   double sign_error = 0.0, zero_error = std::abs(zero_result.energies[0]);
   for(std::size_t index = 0; index < baseline.fields.size(); ++index) {
      sign_error = std::max(sign_error, std::abs(sequence[1].fields[index] + baseline.fields[index]));
      zero_error = std::max(zero_error, std::abs(zero_result.fields[index]));
   }
   const double changed_energy_error = std::abs(sequence[1].energies[0] - sequence[0].energies[0]);
   if(translation_error > tolerance || sign_error > tolerance || zero_error > tolerance || changed_energy_error > tolerance)
      throw std::runtime_error("GPU translation/sign/zero/changed-moment property failed");

   // M=4 is compared to all four independent fields above.  Change only one
   // ensemble here and require every other ensemble to remain bitwise-close.
   const auto& m4 = gpu_ewald_golden_v1::cases[4];
   const Grid m4_grid{{m4.grid[0], m4.grid[1], m4.grid[2]}};
   const auto m4_geometry = periodicGeometry(m4_grid, m4.H);
   const auto m4_kernel = buildPeriodicEwaldDisplacementKernel(m4_geometry, {1.0e-10}).kernel;
   const std::size_t per_ensemble = 3 * cells(m4_grid);
   const std::vector<double> m4_moments(m4.moments, m4.moments + per_ensemble * m4.ensembles);
   auto changed_m4 = m4_moments;
   for(std::size_t index = 2 * per_ensemble; index < 3 * per_ensemble; ++index) changed_m4[index] *= 1.7;
   const auto m4_baseline = evaluateOperator(m4_grid, 1, static_cast<unsigned int>(m4.ensembles), m4.H, m4_kernel, m4_moments);
   const auto m4_changed = evaluateOperator(m4_grid, 1, static_cast<unsigned int>(m4.ensembles), m4.H, m4_kernel, changed_m4);
   double isolation_error = 0.0;
   for(std::size_t ensemble = 0; ensemble < m4.ensembles; ++ensemble) if(ensemble != 2) {
      isolation_error = std::max(isolation_error, std::abs(m4_baseline.energies[ensemble] - m4_changed.energies[ensemble]));
      for(std::size_t index = ensemble * per_ensemble; index < (ensemble + 1) * per_ensemble; ++index)
         isolation_error = std::max(isolation_error, std::abs(m4_baseline.fields[index] - m4_changed.fields[index]));
   }
   if(isolation_error > tolerance) throw std::runtime_error("GPU M=4 ensemble isolation failed");
   std::printf("gpu-properties translation=%.17g sign=%.17g zero=%.17g changed_energy=%.17g ensemble_isolation=%.17g\n",
               translation_error, sign_error, zero_error, changed_energy_error, isolation_error);
}

void runProductionSliceSuite() {
   // This is deliberately the narrow WP5 contract: NA=1, block one, and the
   // complete kernel constructed automatically by the production API.  The
   // expected numbers remain the independent dimensionless oracle fixture;
   // this test applies the physical Tesla/mRy conversion only at its boundary.
   constexpr double mu_b = 9.274009994e-24;
   constexpr double mry = 2.179872325e-21;
   constexpr double tolerance = 3.0e-10;
   double maximum_field_tesla = 0.0, maximum_energy_mry = 0.0;
   double changed_field_error = 0.0, changed_energy_error = 0.0, zero_field_error = 0.0, zero_energy_error = 0.0;
   for(const auto& oracle : gpu_ewald_golden_v1::cases) {
      const Grid grid{{oracle.grid[0], oracle.grid[1], oracle.grid[2]}};
      const std::size_t ncell = cells(grid);
      const std::size_t count = 3 * ncell * oracle.ensembles;
      const double prefactor = 1.0e-7 * mu_b; // alat=1 in descriptor().
      const std::vector<double> moments(oracle.moments, oracle.moments + count);
      // Seed values use the physical scale of this tiny fixture.  A literal
      // Tesla-scale offset would round away every dipole contribution in fp64
      // and would not test additive scatter at all.
      const real field_seed = static_cast<real>(0.125 * prefactor);
      const real total_energy_seed = static_cast<real>(0.375 * prefactor);
      std::vector<real> seeded_fields(count, field_seed);
      std::vector<real> seeded_energy(7 * oracle.ensembles, static_cast<real>(0));
      for(std::size_t ensemble = 0; ensemble < oracle.ensembles; ++ensemble)
         seeded_energy[ensemble + oracle.ensembles * 5] = total_energy_seed;
      std::vector<unsigned int> map(ncell);
      for(std::size_t cell = 0; cell < ncell; ++cell) map[cell] = static_cast<unsigned int>(cell + 1);
      GpuTensor<real, 3> device_moments, beff, eneff;
      GpuTensor<unsigned int, 1> device_map;
      GpuTensor<real, 2> energy;
      GPU_STREAM_T stream{};
      if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("GPU production-slice stream creation failed");
      try {
         device_moments.Allocate(3L, static_cast<long int>(ncell), static_cast<long int>(oracle.ensembles));
         beff.Allocate(3L, static_cast<long int>(ncell), static_cast<long int>(oracle.ensembles));
         eneff.Allocate(3L, static_cast<long int>(ncell), static_cast<long int>(oracle.ensembles));
         device_map.Allocate(static_cast<long int>(ncell));
         energy.Allocate(static_cast<long int>(oracle.ensembles), 7L);
         if(GPU_MEMCPY(device_moments.data(), moments.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(device_map.data(), map.data(), ncell * sizeof(unsigned int), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
            throw std::runtime_error("GPU production-slice input upload failed");
         if(GPU_MEMCPY(beff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(eneff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(energy.data(), seeded_energy.data(), seeded_energy.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
            throw std::runtime_error("GPU production-slice seed upload failed");
         GpuDipoleConvolution solver;
         if(!solver.initiate(descriptor(grid, 1, static_cast<unsigned int>(oracle.ensembles), oracle.H), stream))
            throw std::runtime_error("GPU production-slice plan initiation failed");
         solver.buildPeriodicKernel();
         if(!solver.kernelReady()) throw std::runtime_error("GPU production-slice kernel did not become ready");
         solver.evaluate(device_moments);
         solver.addFieldsToAtoms(beff, eneff, device_map.data(), ncell);
         solver.accumulateEnergy(energy, ncell);
         if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU production-slice synchronization failed");
         std::vector<real> actual_fields(count), actual_eneff(count), actual_energy(7 * oracle.ensembles);
         if(GPU_MEMCPY(actual_fields.data(), beff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
            GPU_MEMCPY(actual_eneff.data(), eneff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
            GPU_MEMCPY(actual_energy.data(), energy.data(), actual_energy.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
            throw std::runtime_error("GPU production-slice output download failed");
         for(std::size_t index = 0; index < count; ++index) {
            const double expected = field_seed + prefactor * oracle.fields[index];
            maximum_field_tesla = std::max(maximum_field_tesla, std::abs(actual_fields[index] - expected));
            maximum_field_tesla = std::max(maximum_field_tesla, std::abs(actual_eneff[index] - expected));
         }
         for(std::size_t ensemble = 0; ensemble < oracle.ensembles; ++ensemble) {
            const double expected_mry = oracle.energies[ensemble] * prefactor * mu_b / mry / static_cast<double>(ncell);
            const double dip_mry = actual_energy[ensemble + oracle.ensembles * 6] * mu_b / mry;
            const double total_mry = actual_energy[ensemble + oracle.ensembles * 5] * mu_b / mry;
            maximum_energy_mry = std::max(maximum_energy_mry, std::abs(dip_mry - expected_mry));
            maximum_energy_mry = std::max(maximum_energy_mry,
                                           std::abs(total_mry - (total_energy_seed * mu_b / mry + expected_mry)));
         }

         // Repeat the exact production scatter/reduction path with changed
         // current moments.  This is the buffer-lifetime analogue of the SD
         // predictor/corrector calls: ensemble zero is sign-flipped while all
         // other ensembles must stay isolated.
         auto changed = moments;
         const std::size_t per_ensemble = 3 * ncell;
         for(std::size_t index = 0; index < per_ensemble; ++index) changed[index] = -changed[index];
         if(GPU_MEMCPY(device_moments.data(), changed.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(beff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(eneff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(energy.data(), seeded_energy.data(), seeded_energy.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
            throw std::runtime_error("GPU production changed-moment upload failed");
         solver.evaluate(device_moments);
         solver.addFieldsToAtoms(beff, eneff, device_map.data(), ncell);
         solver.accumulateEnergy(energy, ncell);
         if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU production changed-moment synchronization failed");
         std::vector<real> changed_fields(count), changed_energy(7 * oracle.ensembles);
         if(GPU_MEMCPY(changed_fields.data(), beff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
            GPU_MEMCPY(changed_energy.data(), energy.data(), changed_energy.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
            throw std::runtime_error("GPU production changed-moment download failed");
         for(std::size_t index = 0; index < count; ++index) {
            const bool changed_ensemble = index < per_ensemble;
            const double expected = field_seed + (changed_ensemble ? -1.0 : 1.0) * prefactor * oracle.fields[index];
            changed_field_error = std::max(changed_field_error, std::abs(changed_fields[index] - expected));
         }
         for(std::size_t ensemble = 0; ensemble < oracle.ensembles; ++ensemble) {
            const double expected = oracle.energies[ensemble] * prefactor / static_cast<double>(ncell);
            changed_energy_error = std::max(changed_energy_error,
                                             std::abs(changed_energy[ensemble + oracle.ensembles * 6] - expected));
            changed_energy_error = std::max(changed_energy_error,
                                             std::abs(changed_energy[ensemble + oracle.ensembles * 5] -
                                                      (total_energy_seed + expected)));
         }

         // A second changed evaluation with a zero first ensemble must not
         // retain either its former field or energy contribution.
         std::fill(changed.begin(), changed.begin() + per_ensemble, 0.0);
         if(GPU_MEMCPY(device_moments.data(), changed.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(beff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(eneff.data(), seeded_fields.data(), count * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(energy.data(), seeded_energy.data(), seeded_energy.size() * sizeof(real), GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS)
            throw std::runtime_error("GPU production zero-moment upload failed");
         solver.evaluate(device_moments);
         solver.addFieldsToAtoms(beff, eneff, device_map.data(), ncell);
         solver.accumulateEnergy(energy, ncell);
         if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("GPU production zero-moment synchronization failed");
         std::vector<real> zero_fields(count), zero_energy(7 * oracle.ensembles);
         if(GPU_MEMCPY(zero_fields.data(), beff.data(), count * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS ||
            GPU_MEMCPY(zero_energy.data(), energy.data(), zero_energy.size() * sizeof(real), GPU_MEMCPY_DEVICE_TO_HOST) != GPU_SUCCESS)
            throw std::runtime_error("GPU production zero-moment download failed");
         for(std::size_t index = 0; index < count; ++index) {
            const bool zeroed_ensemble = index < per_ensemble;
            const double expected = field_seed + (zeroed_ensemble ? 0.0 : prefactor * oracle.fields[index]);
            zero_field_error = std::max(zero_field_error, std::abs(zero_fields[index] - expected));
         }
         for(std::size_t ensemble = 0; ensemble < oracle.ensembles; ++ensemble) {
            const double expected = ensemble == 0 ? 0.0 : oracle.energies[ensemble] * prefactor / static_cast<double>(ncell);
            zero_energy_error = std::max(zero_energy_error,
                                          std::abs(zero_energy[ensemble + oracle.ensembles * 6] - expected));
            zero_energy_error = std::max(zero_energy_error,
                                          std::abs(zero_energy[ensemble + oracle.ensembles * 5] -
                                                   (total_energy_seed + expected)));
         }
         solver.release();
         energy.Free(); device_map.Free(); eneff.Free(); beff.Free(); device_moments.Free();
         if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("GPU production-slice stream destruction failed");
      } catch(...) {
         if(!energy.empty()) energy.Free();
         if(!device_map.empty()) device_map.Free();
         if(!eneff.empty()) eneff.Free();
         if(!beff.empty()) beff.Free();
         if(!device_moments.empty()) device_moments.Free();
         GPU_STREAM_DESTROY(stream);
         throw;
      }
   }
   if(maximum_field_tesla > tolerance * (1.0e-7 * mu_b) ||
      maximum_energy_mry > tolerance * (1.0e-7 * mu_b) * mu_b / mry ||
      changed_field_error > tolerance * (1.0e-7 * mu_b) ||
      changed_energy_error > tolerance * (1.0e-7 * mu_b) ||
      zero_field_error > tolerance * (1.0e-7 * mu_b) ||
      zero_energy_error > tolerance * (1.0e-7 * mu_b)) {
      throw std::runtime_error("GPU production Tesla/mRy acceptance exceeds the WP5 fp64 budget");
   }
   std::printf("production-slice max_field_tesla_error=%.17g max_energy_mry_error=%.17g changed_field_error=%.17g changed_energy_error=%.17g zero_field_error=%.17g zero_energy_error=%.17g\n",
               maximum_field_tesla, maximum_energy_mry, changed_field_error, changed_energy_error,
               zero_field_error, zero_energy_error);
}

} // namespace

int main() {
   try {
      runReciprocalHelperRegression();
      runDeltaSuite();
      runPeriodicSuite();
      runIndependentOracleSuite();
      runMultiBasisOracleSuite();
      runMultiBasisProductionSliceSuite();
      runAlphaAndPropertySuite();
      runProductionSliceSuite();
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU dipole FFT convolution: %s\n", error.what());
      return 1;
   }
   std::puts("PASS GPU dipole FFT delta and periodic convolution suites");
   return 0;
}
