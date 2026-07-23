#include "dipoleEwaldKernel.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

namespace {

using Vector = std::array<double, 3>;

void close(double actual, double expected, double tolerance = 3.0e-9) {
   if(std::abs(actual - expected) > tolerance * std::max({1.0, std::abs(actual), std::abs(expected)})) {
      throw std::runtime_error("host Ewald kernel differs from Luna's independent oracle value");
   }
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

DipolePeriodicGeometry geometry(const std::array<double, 9>& h, std::array<std::size_t, 3> grid) {
   DipolePeriodicGeometry result{};
   result.H = h;
   result.volume = volume(h);
   result.Brec = reciprocal(h, result.volume);
   result.grid = grid;
   result.basis = 1;
   result.basis_offsets = {{0.0, 0.0, 0.0}};
   return result;
}

std::vector<Vector> apply(const DipoleKernelBuildResult& result, const DipolePeriodicGeometry& geometry,
                          const std::vector<Vector>& moments) {
   const auto [n1, n2, n3] = geometry.grid;
   const std::size_t cells = n1 * n2 * n3;
   std::vector<Vector> fields(cells, {0.0, 0.0, 0.0});
   for(std::size_t z = 0; z < n3; ++z) for(std::size_t y = 0; y < n2; ++y) for(std::size_t x = 0; x < n1; ++x) {
      const std::size_t target = x + n1 * (y + n2 * z);
      for(std::size_t sz = 0; sz < n3; ++sz) for(std::size_t sy = 0; sy < n2; ++sy) for(std::size_t sx = 0; sx < n1; ++sx) {
         const std::size_t source = sx + n1 * (sy + n2 * sz);
         const std::size_t dx = (x + n1 - sx) % n1;
         const std::size_t dy = (y + n2 - sy) % n2;
         const std::size_t dz = (z + n3 - sz) % n3;
         const std::size_t cell = dx + n1 * (dy + n2 * dz);
         for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
            const std::size_t batch = row + 3 * column;
            fields[target][row] += result.kernel[cell + cells * batch] * moments[source][column];
         }
      }
   }
   return fields;
}

void compareFields(const std::vector<Vector>& actual, const std::vector<Vector>& expected) {
   if(actual.size() != expected.size()) throw std::runtime_error("wrong number of host Ewald fields");
   for(std::size_t i = 0; i < actual.size(); ++i)
      for(unsigned int component = 0; component < 3; ++component) close(actual[i][component], expected[i][component]);
}

void checkDiagnostics(const DipoleKernelBuildResult& result, const DipolePeriodicGeometry& geom) {
   if(result.kernel.size() != geom.grid[0] * geom.grid[1] * geom.grid[2] * 9) {
      throw std::runtime_error("host Ewald kernel packing has the wrong size");
   }
   if(result.diagnostics.selected.alpha <= 0.0 || result.diagnostics.setup_work == 0 ||
      result.diagnostics.max_reciprocity_error > 2.0e-11 || result.diagnostics.max_hermitian_error > 2.0e-11 ||
      result.diagnostics.max_alpha_difference > 3.0e-10 ||
      result.diagnostics.real_tail_residual > 3.0e-10 || result.diagnostics.reciprocal_tail_residual > 3.0e-10) {
      throw std::runtime_error("host Ewald convergence, alpha-invariance, or reciprocity diagnostics failed");
   }
   close(periodicKernelReciprocityError(result.kernel, geom), result.diagnostics.max_reciprocity_error, 1.0e-13);
}

} // namespace

int main() {
   try {
      const DipoleKernelSettings settings{1.0e-10};

      // Luna's independently evaluated cubic 1x1x1 point-dipole reference.
      auto cubic = geometry({3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0}, {1, 1, 1});
      const auto cubic_result = buildPeriodicEwaldDisplacementKernel(cubic, settings);
      checkDiagnostics(cubic_result, cubic);
      compareFields(apply(cubic_result, cubic, {{1.0, -0.2, 0.4}}),
                    {{0.1551403779550518, -0.031028075591010243, 0.06205615118202068}});

      // Non-cubic 2x3x1 reference, with n1-fast cell ordering.
      auto non_cubic = geometry({6.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 3.0}, {2, 3, 1});
      const auto non_cubic_result = buildPeriodicEwaldDisplacementKernel(non_cubic, settings);
      checkDiagnostics(non_cubic_result, non_cubic);
      compareFields(apply(non_cubic_result, non_cubic,
                          {{0.20, -0.40, 0.10}, {0.23, -0.38, 0.09}, {0.26, -0.36, 0.08},
                           {0.29, -0.34, 0.07}, {0.32, -0.32, 0.06}, {0.35, -0.30, 0.05}}),
                    {{0.03764369536701564, -0.042224008061581005, 0.015999787603181142},
                     {0.026879869768327465, -0.03863606619530979, 0.0142058166700661},
                     {0.04804551673699256, -0.05609310321813055, 0.012532513813188912},
                     {0.03728169113830427, -0.05250516135185944, 0.01073854288007384},
                     {0.05844733810696956, -0.0699621983746795, 0.009065240023196796},
                     {0.047683512508281026, -0.06637425650840889, 0.0072712690900815224}});

      // A skew H must use its full metric in both direct and reciprocal sums.
      auto skew = geometry({6.0, 0.0, 0.0, 0.8, 3.0, 0.0, 0.3, 0.2, 3.4}, {2, 1, 1});
      const auto skew_result = buildPeriodicEwaldDisplacementKernel(skew, settings);
      checkDiagnostics(skew_result, skew);
      compareFields(apply(skew_result, skew, {{1.0, -0.2, 0.4}, {-0.3, 0.8, 0.1}}),
                    {{-0.14844761149081095, 0.04249102731392704, 0.05706032893486248},
                     {0.2799443305183004, 0.0559736351960156, -0.010870106834032328}});
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL host periodic Ewald kernel: %s\n", error.what());
      return 1;
   }
   std::puts("PASS host periodic Ewald Builder A");
   return 0;
}
