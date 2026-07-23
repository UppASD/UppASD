#include "dipoleEwaldKernel.hpp"
#include "gpuDipoleConvolution.hpp"
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
           scale * (h[3] * h[7] - h[4] * h[6]), scale * (h[7] * h[2] - h[8] * h[5]),
           scale * (h[8] * h[0] - h[6] * h[2]), scale * (h[6] * h[1] - h[7] * h[0]),
           scale * (h[1] * h[5] - h[2] * h[4]), scale * (h[2] * h[3] - h[0] * h[5]),
           scale * (h[0] * h[4] - h[1] * h[3])};
}

double volume(const std::array<double, 9>& h) {
   return h[0] * (h[4] * h[8] - h[7] * h[5]) - h[3] * (h[1] * h[8] - h[7] * h[2]) +
          h[6] * (h[1] * h[5] - h[4] * h[2]);
}

GpuDipoleConvolutionDescriptor descriptor(const Grid& grid, unsigned int basis, unsigned int ensembles,
                                          const std::array<double, 9>& h) {
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

Errors runOperator(const std::string& name, const Grid& grid, unsigned int basis, unsigned int ensembles,
                   const std::array<double, 9>& h, const std::vector<double>& kernel,
                   const std::vector<double>& moments, bool check_not_ready = false) {
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
      solver.uploadCompleteKernelForTesting(kernel);
      if(!solver.kernelReady()) throw std::runtime_error("kernel_ready was not set after normalized kernel construction");
      if(solver.diagnosticConstructionStorageAllocatedForTesting())
         throw std::runtime_error("transient real kernel storage survived kernel construction");
      solver.evaluate(device_moments);
      const auto fields = solver.diagnosticFieldsForTesting();
      const auto energies = solver.diagnosticEnergiesForTesting();
      const auto expected = directConvolution(kernel, moments, grid, basis, ensembles);
      Errors errors{};
      for(std::size_t i = 0; i < fields.size(); ++i) errors.field = std::max(errors.field, std::abs(fields[i] - expected[i]));
      for(unsigned int ensemble = 0; ensemble < ensembles; ++ensemble) {
         double expected_energy = 0.0;
         for(std::size_t cell = 0; cell < ncell; ++cell) for(unsigned int channel = 0; channel < basis; ++channel)
            for(unsigned int component = 0; component < 3; ++component) {
               const auto index = macroIndex(component, channel, cell, basis, ncell, ensemble);
               expected_energy -= 0.5 * moments[index] * expected[index];
            }
         errors.energy = std::max(errors.energy, std::abs(energies[ensemble] - expected_energy));
      }
      solver.release();
      device_moments.Free();
      if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("GPU stream destruction failed");
      std::printf("%s max_field_error=%.17g max_energy_error=%.17g\n", name.c_str(), errors.field, errors.energy);
      return errors;
   } catch(...) {
      if(!device_moments.empty()) device_moments.Free();
      GPU_STREAM_DESTROY(stream);
      throw;
   }
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
                             ensembles == 1 && source == 0 && component == 0), name);
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

} // namespace

int main() {
   try {
      runDeltaSuite();
      runPeriodicSuite();
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU dipole FFT convolution: %s\n", error.what());
      return 1;
   }
   std::puts("PASS GPU dipole FFT delta and periodic convolution suites");
   return 0;
}
