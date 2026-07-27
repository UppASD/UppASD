#include "dipoleEwaldKernel.hpp"
#include "gpu_ewald_multibasis_goldens_v1.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
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
   if(moments.size() != cells * geometry.basis) throw std::runtime_error("host Ewald moment layout has the wrong basis count");
   std::vector<Vector> fields(cells * geometry.basis, {0.0, 0.0, 0.0});
   for(std::size_t z = 0; z < n3; ++z) for(std::size_t y = 0; y < n2; ++y) for(std::size_t x = 0; x < n1; ++x) {
      const std::size_t target = x + n1 * (y + n2 * z);
      for(std::size_t sz = 0; sz < n3; ++sz) for(std::size_t sy = 0; sy < n2; ++sy) for(std::size_t sx = 0; sx < n1; ++sx) {
         const std::size_t source = sx + n1 * (sy + n2 * sz);
         const std::size_t dx = (x + n1 - sx) % n1;
         const std::size_t dy = (y + n2 - sy) % n2;
         const std::size_t dz = (z + n3 - sz) % n3;
         const std::size_t cell = dx + n1 * (dy + n2 * dz);
         for(unsigned int target_basis = 0; target_basis < geometry.basis; ++target_basis)
            for(unsigned int source_basis = 0; source_basis < geometry.basis; ++source_basis)
               for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
                  const std::size_t batch = row + 3 * (column + 3 * (target_basis + geometry.basis * source_basis));
                  fields[target_basis + geometry.basis * target][row] +=
                     result.kernel[cell + cells * batch] * moments[source_basis + geometry.basis * source][column];
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
   if(result.kernel.size() != geom.grid[0] * geom.grid[1] * geom.grid[2] * 9 * geom.basis * geom.basis) {
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

std::vector<std::complex<double>> completeAliasSpectrum(const DipoleAliasSpectrumBuildResult& result,
                                                        const DipolePeriodicGeometry& geom) {
   auto spectrum = referenceNormalizedKernelSpectrum(result.real_kernel, geom);
   if(spectrum.size() != result.reciprocal_alias_spectrum.size()) {
      throw std::runtime_error("Builder B returned a reciprocal alias spectrum with the wrong shape");
   }
   for(std::size_t index = 0; index < spectrum.size(); ++index) spectrum[index] += result.reciprocal_alias_spectrum[index];
   return spectrum;
}

double maxDifference(const std::vector<std::complex<double>>& left, const std::vector<std::complex<double>>& right) {
   if(left.size() != right.size()) throw std::runtime_error("cross-builder spectra have different shapes");
   double maximum = 0.0;
   for(std::size_t index = 0; index < left.size(); ++index) maximum = std::max(maximum, std::abs(left[index] - right[index]));
   return maximum;
}

void checkAliasBuilder(const DipoleKernelBuildResult& reference, const DipolePeriodicGeometry& geom,
                       const DipoleKernelSettings& settings) {
   const auto alias = buildPeriodicEwaldAliasSpectrum(geom, settings);
   if(alias.real_kernel.size() != reference.kernel.size() || alias.reciprocal_alias_spectrum.empty() ||
      alias.diagnostics.selected.alpha <= 0.0 || alias.diagnostics.setup_work == 0 ||
      alias.diagnostics.real_tail_residual > 3.0e-10 || alias.diagnostics.reciprocal_tail_residual > 3.0e-10 ||
      alias.diagnostics.max_reciprocity_error > 2.0e-11) {
      throw std::runtime_error("Builder B automatic construction diagnostics failed");
   }
   const double error = maxDifference(referenceNormalizedKernelSpectrum(reference.kernel, geom),
                                      completeAliasSpectrum(alias, geom));
   if(error > 4.0e-10) throw std::runtime_error("Builder A/B spectral equivalence exceeds the fp64 budget");

   // Explicit alphas remain test-only, but Builder B must preserve the same
   // split invariance as the accepted reference builder.
   const auto low = buildPeriodicEwaldAliasSpectrum(geom, {settings.tolerance, 0.70});
   const auto high = buildPeriodicEwaldAliasSpectrum(geom, {settings.tolerance, 1.05});
   const double alpha_error = maxDifference(completeAliasSpectrum(low, geom), completeAliasSpectrum(high, geom));
   if(alpha_error > 4.0e-10) throw std::runtime_error("Builder B explicit-alpha invariance exceeds the fp64 budget");
   std::printf("Builder A/B benchmark: spectrum=%.3e alpha=%.3e A=%.6fs/%zu B=%.6fs/%zu real_bytes=%zu alias_bytes=%zu\n",
               error, alpha_error, reference.diagnostics.setup_seconds, reference.diagnostics.setup_work,
               alias.diagnostics.setup_seconds, alias.diagnostics.setup_work, alias.real_kernel.size() * sizeof(double),
               alias.reciprocal_alias_spectrum.size() * sizeof(std::complex<double>));
}

} // namespace

int main() {
   try {
      const DipoleKernelSettings settings{1.0e-10};

      // Luna's independently evaluated cubic 1x1x1 point-dipole reference.
      auto cubic = geometry({3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0}, {1, 1, 1});
      const auto cubic_result = buildPeriodicEwaldDisplacementKernel(cubic, settings);
      checkDiagnostics(cubic_result, cubic);
      checkAliasBuilder(cubic_result, cubic, settings);
      compareFields(apply(cubic_result, cubic, {{1.0, -0.2, 0.4}}),
                    {{0.1551403779550518, -0.031028075591010243, 0.06205615118202068}});

      // Non-cubic 2x3x1 reference, with n1-fast cell ordering.
      auto non_cubic = geometry({6.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 3.0}, {2, 3, 1});
      const auto non_cubic_result = buildPeriodicEwaldDisplacementKernel(non_cubic, settings);
      checkDiagnostics(non_cubic_result, non_cubic);
      checkAliasBuilder(non_cubic_result, non_cubic, settings);
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
      checkAliasBuilder(skew_result, skew, settings);
      compareFields(apply(skew_result, skew, {{1.0, -0.2, 0.4}, {-0.3, 0.8, 0.1}}),
                    {{-0.14844761149081095, 0.04249102731392704, 0.05706032893486248},
                     {0.2799443305183004, 0.0559736351960156, -0.010870106834032328}});

      // Independent L1_0-style two-basis fixture.  The second case is a
      // deliberately distorted/skew supercell, so this covers basis phases
      // and the triclinic reciprocal metric in the CPU-only Builder-A gate.
      for(const auto& oracle : gpu_ewald_multibasis_golden_v1::cases) {
         auto multi = geometry(oracle.H, oracle.grid);
         multi.basis = static_cast<unsigned int>(oracle.basis);
         multi.basis_offsets.resize(multi.basis);
         for(unsigned int basis = 0; basis < multi.basis; ++basis)
            for(unsigned int component = 0; component < 3; ++component)
               multi.basis_offsets[basis][component] = oracle.offsets[component + 3 * basis];
         const auto multi_result = buildPeriodicEwaldDisplacementKernel(multi, settings);
         checkDiagnostics(multi_result, multi);
         checkAliasBuilder(multi_result, multi, settings);
         const std::size_t count = multi.grid[0] * multi.grid[1] * multi.grid[2] * multi.basis;
         std::vector<Vector> moments(count);
         std::vector<Vector> expected(count);
         for(std::size_t index = 0; index < count; ++index) {
            for(unsigned int component = 0; component < 3; ++component) {
               moments[index][component] = oracle.moments[3 * index + component];
               expected[index][component] = oracle.fields[3 * index + component];
            }
         }
         compareFields(apply(multi_result, multi, moments), expected);
      }
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL host periodic Ewald kernel: %s\n", error.what());
      return 1;
   }
   std::puts("PASS host periodic Ewald Builder A");
   return 0;
}
