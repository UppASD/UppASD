#include "dipoleEwaldKernel.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double minimum_tolerance = 128.0 * std::numeric_limits<double>::epsilon();

using Vector = std::array<double, 3>;
using Tensor = std::array<double, 9>;
using Complex = std::complex<double>;

struct CandidateResult {
   double alpha = 0.0;
   int real_extent = 0;
   int reciprocal_extent = 0;
   double real_residual = 0.0;
   double reciprocal_residual = 0.0;
   std::size_t real_work = 0;
   std::size_t reciprocal_work = 0;
};

std::size_t cells(const DipolePeriodicGeometry& geometry) {
   const auto [n1, n2, n3] = geometry.grid;
   if(n1 == 0 || n2 == 0 || n3 == 0 || n1 > std::numeric_limits<std::size_t>::max() / n2 ||
      n1 * n2 > std::numeric_limits<std::size_t>::max() / n3) {
      throw std::invalid_argument("periodic Ewald kernel requires a nonzero non-overflowing grid");
   }
   return n1 * n2 * n3;
}

Vector add(const Vector& left, const Vector& right) {
   return {left[0] + right[0], left[1] + right[1], left[2] + right[2]};
}

Vector subtract(const Vector& left, const Vector& right) {
   return {left[0] - right[0], left[1] - right[1], left[2] - right[2]};
}

double dot(const Vector& left, const Vector& right) {
   return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

Vector matrixVector(const std::array<double, 9>& matrix, const Vector& vector) {
   return {matrix[0] * vector[0] + matrix[3] * vector[1] + matrix[6] * vector[2],
           matrix[1] * vector[0] + matrix[4] * vector[1] + matrix[7] * vector[2],
           matrix[2] * vector[0] + matrix[5] * vector[1] + matrix[8] * vector[2]};
}

double maxAbs(const Tensor& tensor) {
   double result = 0.0;
   for(const double value : tensor) result = std::max(result, std::abs(value));
   return result;
}

void addTensor(Tensor& target, const Tensor& source) {
   for(std::size_t i = 0; i < target.size(); ++i) target[i] += source[i];
}

Tensor realTensor(const Vector& r, double alpha) {
   const double r2 = dot(r, r);
   if(r2 == 0.0) throw std::logic_error("zero displacement was not excluded from the real Ewald sum");
   const double radius = std::sqrt(r2);
   const double gaussian = std::exp(-alpha * alpha * r2);
   const double erfc = std::erfc(alpha * radius);
   const double root_pi = std::sqrt(pi);
   const double a = 3.0 * erfc / (r2 * r2 * radius) +
                    6.0 * alpha * gaussian / (root_pi * r2 * r2) +
                    4.0 * alpha * alpha * alpha * gaussian / (root_pi * r2);
   const double b = -erfc / (r2 * radius) - 2.0 * alpha * gaussian / (root_pi * r2);
   Tensor result{};
   for(unsigned int row = 0; row < 3; ++row)
      for(unsigned int column = 0; column < 3; ++column)
         result[row + 3 * column] = a * r[row] * r[column] + (row == column ? b : 0.0);
   return result;
}

Tensor reciprocalTensor(const Vector& wave, double volume, double alpha) {
   const double wave2 = dot(wave, wave);
   if(wave2 == 0.0) throw std::logic_error("physical reciprocal k=0 was not excluded");
   const double factor = -4.0 * pi * std::exp(-wave2 / (4.0 * alpha * alpha)) / (volume * wave2);
   Tensor result{};
   for(unsigned int row = 0; row < 3; ++row)
      for(unsigned int column = 0; column < 3; ++column)
         result[row + 3 * column] = factor * wave[row] * wave[column];
   return result;
}

double determinant(const std::array<double, 9>& matrix) {
   return matrix[0] * (matrix[4] * matrix[8] - matrix[7] * matrix[5]) -
          matrix[3] * (matrix[1] * matrix[8] - matrix[7] * matrix[2]) +
          matrix[6] * (matrix[1] * matrix[5] - matrix[4] * matrix[2]);
}

double reciprocalIdentityError(const DipolePeriodicGeometry& geometry) {
   double maximum = 0.0;
   for(unsigned int row = 0; row < 3; ++row) {
      for(unsigned int column = 0; column < 3; ++column) {
         double value = 0.0;
         for(unsigned int component = 0; component < 3; ++component)
            value += geometry.H[component + 3 * row] * geometry.Brec[component + 3 * column];
         maximum = std::max(maximum, std::abs(value - (row == column ? 2.0 * pi : 0.0)));
      }
   }
   return maximum;
}

void validate(const DipolePeriodicGeometry& geometry, const DipoleKernelSettings& settings) {
   const std::size_t grid_cells = cells(geometry);
   (void)grid_cells;
   if(geometry.basis == 0 || geometry.basis_offsets.size() != geometry.basis) {
      throw std::invalid_argument("periodic Ewald kernel basis offsets must contain one entry per basis channel");
   }
   if(!std::isfinite(settings.tolerance) || settings.tolerance < minimum_tolerance) {
      throw std::invalid_argument("gpu_dipole_tol must be finite and no smaller than fp64 construction headroom");
   }
   if(!std::isfinite(settings.explicit_alpha_for_testing)) {
      throw std::invalid_argument("test-only Ewald alpha must be finite");
   }
   if(!std::isfinite(geometry.volume) || geometry.volume <= 0.0 || !std::isfinite(determinant(geometry.H)) ||
      determinant(geometry.H) <= 0.0) {
      throw std::invalid_argument("periodic Ewald kernel requires a right-handed positive-volume H");
   }
   const double volume_error = std::abs(determinant(geometry.H) - geometry.volume);
   if(volume_error > 128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, geometry.volume)) {
      throw std::invalid_argument("periodic Ewald kernel volume disagrees with det(H)");
   }
   for(const double value : geometry.H) if(!std::isfinite(value)) throw std::invalid_argument("H must be finite");
   for(const double value : geometry.Brec) if(!std::isfinite(value)) throw std::invalid_argument("Brec must be finite");
   for(const auto& offset : geometry.basis_offsets)
      for(const double value : offset) if(!std::isfinite(value)) throw std::invalid_argument("basis offsets must be finite");
   const double identity_error = reciprocalIdentityError(geometry);
   if(identity_error > 4096.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, maxAbs(geometry.H))) {
      throw std::invalid_argument("Brec is inconsistent with H^(-T) * 2*pi");
   }
}

Vector displacement(const DipolePeriodicGeometry& geometry, std::size_t d1, std::size_t d2, std::size_t d3,
                    unsigned int target_basis, unsigned int source_basis) {
   const Vector fractional = {static_cast<double>(d1) / static_cast<double>(geometry.grid[0]),
                              static_cast<double>(d2) / static_cast<double>(geometry.grid[1]),
                              static_cast<double>(d3) / static_cast<double>(geometry.grid[2])};
   return add(matrixVector(geometry.H, fractional),
              subtract(geometry.basis_offsets[target_basis], geometry.basis_offsets[source_basis]));
}

Tensor sumReal(const DipolePeriodicGeometry& geometry, const Vector& displacement_vector, double alpha, int extent,
               std::size_t& work) {
   Tensor result{};
   for(int n3 = -extent; n3 <= extent; ++n3) for(int n2 = -extent; n2 <= extent; ++n2)
      for(int n1 = -extent; n1 <= extent; ++n1) {
         const Vector image = matrixVector(geometry.H, {static_cast<double>(n1), static_cast<double>(n2), static_cast<double>(n3)});
         const Vector r = subtract(displacement_vector, image);
         if(dot(r, r) == 0.0) continue;
         addTensor(result, realTensor(r, alpha));
         ++work;
      }
   return result;
}

Tensor sumReciprocal(const DipolePeriodicGeometry& geometry, const Vector& displacement_vector, double alpha,
                     int extent, std::size_t& work) {
   Tensor result{};
   for(int n3 = -extent; n3 <= extent; ++n3) for(int n2 = -extent; n2 <= extent; ++n2)
      for(int n1 = -extent; n1 <= extent; ++n1) {
         if(n1 == 0 && n2 == 0 && n3 == 0) continue; // Tin foil omits physical k=0 only.
         const Vector wave = matrixVector(geometry.Brec, {static_cast<double>(n1), static_cast<double>(n2), static_cast<double>(n3)});
         const double phase = std::cos(dot(wave, displacement_vector));
         Tensor term = reciprocalTensor(wave, geometry.volume, alpha);
         for(double& value : term) value *= phase;
         addTensor(result, term);
         ++work;
      }
   return result;
}

std::vector<double> sumSide(const DipolePeriodicGeometry& geometry, double alpha, int extent, bool real,
                            std::size_t& work) {
   const std::size_t grid_cells = cells(geometry);
   const std::size_t batches = 9 * static_cast<std::size_t>(geometry.basis) * geometry.basis;
   std::vector<double> result(grid_cells * batches, 0.0);
   for(std::size_t d3 = 0; d3 < geometry.grid[2]; ++d3) for(std::size_t d2 = 0; d2 < geometry.grid[1]; ++d2)
      for(std::size_t d1 = 0; d1 < geometry.grid[0]; ++d1) {
         const std::size_t cell = d1 + geometry.grid[0] * (d2 + geometry.grid[1] * d3);
         for(unsigned int target = 0; target < geometry.basis; ++target)
            for(unsigned int source = 0; source < geometry.basis; ++source) {
               const Vector r = displacement(geometry, d1, d2, d3, target, source);
               const Tensor tensor = real ? sumReal(geometry, r, alpha, extent, work)
                                          : sumReciprocal(geometry, r, alpha, extent, work);
               for(unsigned int column = 0; column < 3; ++column) for(unsigned int row = 0; row < 3; ++row) {
                  const std::size_t batch = row + 3 * (column + 3 * (target + geometry.basis * source));
                  result[cell + grid_cells * batch] = tensor[row + 3 * column];
               }
            }
      }
   return result;
}

double vectorMaxAbs(const std::vector<double>& values) {
   double result = 0.0;
   for(const double value : values) result = std::max(result, std::abs(value));
   return result;
}

double vectorMaxDifference(const std::vector<double>& left, const std::vector<double>& right) {
   if(left.size() != right.size()) throw std::logic_error("Ewald side sizes differ");
   double result = 0.0;
   for(std::size_t i = 0; i < left.size(); ++i) result = std::max(result, std::abs(left[i] - right[i]));
   return result;
}

// A lower bound for sigma_min(A) is 1 / ||A^-1||_F.  Unlike a Cartesian box
// length or individual column norm, it accounts for all skew/cancellation in
// the full lattice matrix.  Shell convergence still measures actual tensor
// changes, independently per side.
double fullMatrixLengthLowerBound(const std::array<double, 9>& matrix) {
   const double det = determinant(matrix);
   if(!std::isfinite(det) || det == 0.0) return 0.0;
   const std::array<double, 9> cofactors{{
      matrix[4] * matrix[8] - matrix[7] * matrix[5],
      matrix[7] * matrix[2] - matrix[1] * matrix[8],
      matrix[1] * matrix[5] - matrix[4] * matrix[2],
      matrix[6] * matrix[5] - matrix[3] * matrix[8],
      matrix[0] * matrix[8] - matrix[6] * matrix[2],
      matrix[3] * matrix[2] - matrix[0] * matrix[5],
      matrix[3] * matrix[7] - matrix[6] * matrix[4],
      matrix[6] * matrix[1] - matrix[0] * matrix[7],
      matrix[0] * matrix[4] - matrix[3] * matrix[1]}};
   double inverse_frobenius_squared = 0.0;
   for(const double cofactor : cofactors) inverse_frobenius_squared += (cofactor / det) * (cofactor / det);
   return inverse_frobenius_squared > 0.0 ? 1.0 / std::sqrt(inverse_frobenius_squared) : 0.0;
}

int analyticSeed(const std::array<double, 9>& matrix, double exponent_scale, double tolerance) {
   const double shortest = fullMatrixLengthLowerBound(matrix);
   if(!std::isfinite(shortest) || shortest <= 0.0) return 0;
   return std::max(0, static_cast<int>(std::ceil(std::sqrt(-std::log(tolerance)) / (exponent_scale * shortest))) - 1);
}

CandidateResult convergeCandidate(const DipolePeriodicGeometry& geometry, double alpha, double tolerance) {
   constexpr int max_shell = 64;
   const double side_tolerance = tolerance / 4.0;
   CandidateResult result{};
   result.alpha = alpha;
   const auto converge = [&](bool real, int seed, int& selected_extent, double& residual, std::size_t& work) {
      std::vector<double> previous;
      int stable = 0;
      for(int extent = std::max(0, seed - 1); extent <= max_shell; ++extent) {
         std::vector<double> current = sumSide(geometry, alpha, extent, real, work);
         const double change = previous.empty() ? std::numeric_limits<double>::infinity() : vectorMaxDifference(current, previous);
         const double scale = vectorMaxAbs(current);
         if(!previous.empty() && change <= side_tolerance * std::max(1.0, scale)) {
            ++stable;
            if(stable == 2) {
               selected_extent = extent;
               residual = change;
               return current;
            }
         } else {
            stable = 0;
         }
         previous = std::move(current);
      }
      throw std::runtime_error(real ? "real Ewald shells did not converge" : "reciprocal Ewald shells did not converge");
   };
   const int real_seed = analyticSeed(geometry.H, alpha, side_tolerance);
   const int reciprocal_seed = analyticSeed(geometry.Brec, 1.0 / (2.0 * alpha), side_tolerance);
   const auto real_values = converge(true, real_seed, result.real_extent, result.real_residual, result.real_work);
   const auto reciprocal_values = converge(false, reciprocal_seed, result.reciprocal_extent, result.reciprocal_residual,
                                           result.reciprocal_work);
   (void)real_values;
   (void)reciprocal_values;
   return result;
}

std::vector<double> completeKernel(const DipolePeriodicGeometry& geometry, const CandidateResult& candidate,
                                   std::size_t& real_work, std::size_t& reciprocal_work) {
   auto real = sumSide(geometry, candidate.alpha, candidate.real_extent, true, real_work);
   const auto reciprocal = sumSide(geometry, candidate.alpha, candidate.reciprocal_extent, false, reciprocal_work);
   for(std::size_t i = 0; i < real.size(); ++i) real[i] += reciprocal[i];
   const std::size_t grid_cells = cells(geometry);
   const double self = 4.0 * candidate.alpha * candidate.alpha * candidate.alpha / (3.0 * std::sqrt(pi));
   const std::size_t zero_cell = 0;
   for(unsigned int basis = 0; basis < geometry.basis; ++basis)
      for(unsigned int component = 0; component < 3; ++component) {
         const std::size_t batch = component + 3 * (component + 3 * (basis + geometry.basis * basis));
         real[zero_cell + grid_cells * batch] += self;
      }
   return real;
}

std::vector<double> realSelfKernel(const DipolePeriodicGeometry& geometry, const CandidateResult& candidate,
                                   std::size_t& real_work) {
   auto real = sumSide(geometry, candidate.alpha, candidate.real_extent, true, real_work);
   const std::size_t grid_cells = cells(geometry);
   const double self = 4.0 * candidate.alpha * candidate.alpha * candidate.alpha / (3.0 * std::sqrt(pi));
   const std::size_t zero_cell = 0;
   for(unsigned int basis = 0; basis < geometry.basis; ++basis)
      for(unsigned int component = 0; component < 3; ++component) {
         const std::size_t batch = component + 3 * (component + 3 * (basis + geometry.basis * basis));
         real[zero_cell + grid_cells * batch] += self;
      }
   return real;
}

std::array<double, 9> aliasReciprocalMatrix(const DipolePeriodicGeometry& geometry) {
   std::array<double, 9> result = geometry.Brec;
   for(unsigned int column = 0; column < 3; ++column)
      for(unsigned int row = 0; row < 3; ++row)
         result[row + 3 * column] *= static_cast<double>(geometry.grid[column]);
   return result;
}

int signedFrequency(std::size_t q, std::size_t extent) {
   return static_cast<int>(q <= extent / 2 ? q : q - extent);
}

std::size_t spectralCells(const DipolePeriodicGeometry& geometry) {
   const std::size_t n1 = geometry.grid[0] / 2 + 1;
   if(n1 > std::numeric_limits<std::size_t>::max() / geometry.grid[1] ||
      n1 * geometry.grid[1] > std::numeric_limits<std::size_t>::max() / geometry.grid[2]) {
      throw std::invalid_argument("periodic Ewald alias spectrum dimensions overflow");
   }
   return n1 * geometry.grid[1] * geometry.grid[2];
}

std::vector<Complex> sumReciprocalAliases(const DipolePeriodicGeometry& geometry, double alpha, int extent,
                                          std::size_t& work) {
   const std::size_t cells_spectral = spectralCells(geometry);
   const std::size_t batches = 9 * static_cast<std::size_t>(geometry.basis) * geometry.basis;
   std::vector<Complex> result(cells_spectral * batches, Complex{});
   const std::size_t stored_n1 = geometry.grid[0] / 2 + 1;
   for(std::size_t q3 = 0; q3 < geometry.grid[2]; ++q3)
      for(std::size_t q2 = 0; q2 < geometry.grid[1]; ++q2)
         for(std::size_t q1 = 0; q1 < stored_n1; ++q1) {
            const std::size_t spectral = q1 + stored_n1 * (q2 + geometry.grid[1] * q3);
            const std::array<int, 3> q{{static_cast<int>(q1), signedFrequency(q2, geometry.grid[1]),
                                        signedFrequency(q3, geometry.grid[2])}};
            for(unsigned int target = 0; target < geometry.basis; ++target)
               for(unsigned int source = 0; source < geometry.basis; ++source) {
                  const Vector offset = subtract(geometry.basis_offsets[target], geometry.basis_offsets[source]);
                  for(int l3 = -extent; l3 <= extent; ++l3)
                     for(int l2 = -extent; l2 <= extent; ++l2)
                        for(int l1 = -extent; l1 <= extent; ++l1) {
                           const std::array<int, 3> n{{q[0] + static_cast<int>(geometry.grid[0]) * l1,
                                                       q[1] + static_cast<int>(geometry.grid[1]) * l2,
                                                       q[2] + static_cast<int>(geometry.grid[2]) * l3}};
                           if(n[0] == 0 && n[1] == 0 && n[2] == 0) continue; // Tin foil: physical k=0 only.
                           const Vector wave = matrixVector(geometry.Brec,
                              {static_cast<double>(n[0]), static_cast<double>(n[1]), static_cast<double>(n[2])});
                           const Complex phase = std::polar(1.0, dot(wave, offset));
                           const Tensor tensor = reciprocalTensor(wave, geometry.volume, alpha);
                           for(unsigned int column = 0; column < 3; ++column)
                              for(unsigned int row = 0; row < 3; ++row) {
                                 const std::size_t batch = row + 3 * (column + 3 * (target + geometry.basis * source));
                                 result[spectral + cells_spectral * batch] += tensor[row + 3 * column] * phase;
                              }
                           ++work;
                        }
               }
         }
   return result;
}

double complexVectorMaxAbs(const std::vector<Complex>& values) {
   double result = 0.0;
   for(const Complex& value : values) result = std::max(result, std::abs(value));
   return result;
}

double complexVectorMaxDifference(const std::vector<Complex>& left, const std::vector<Complex>& right) {
   if(left.size() != right.size()) throw std::logic_error("reciprocal alias spectrum sizes differ");
   double result = 0.0;
   for(std::size_t i = 0; i < left.size(); ++i) result = std::max(result, std::abs(left[i] - right[i]));
   return result;
}

CandidateResult convergeAliasCandidate(const DipolePeriodicGeometry& geometry, double alpha, double tolerance) {
   constexpr int max_shell = 64;
   const double side_tolerance = tolerance / 4.0;
   CandidateResult result{};
   result.alpha = alpha;
   const auto converge_real = [&](int seed) {
      std::vector<double> previous;
      int stable = 0;
      for(int extent = std::max(0, seed - 1); extent <= max_shell; ++extent) {
         std::vector<double> current = sumSide(geometry, alpha, extent, true, result.real_work);
         const double change = previous.empty() ? std::numeric_limits<double>::infinity() :
            vectorMaxDifference(current, previous);
         if(!previous.empty() && change <= side_tolerance * std::max(1.0, vectorMaxAbs(current))) {
            if(++stable == 2) {
               result.real_extent = extent;
               result.real_residual = change;
               return;
            }
         } else stable = 0;
         previous = std::move(current);
      }
      throw std::runtime_error("real Ewald shells did not converge for reciprocal-alias builder");
   };
   const auto converge_reciprocal = [&](int seed) {
      std::vector<Complex> previous;
      int stable = 0;
      for(int extent = std::max(0, seed - 1); extent <= max_shell; ++extent) {
         std::vector<Complex> current = sumReciprocalAliases(geometry, alpha, extent, result.reciprocal_work);
         const double change = previous.empty() ? std::numeric_limits<double>::infinity() :
            complexVectorMaxDifference(current, previous);
         if(!previous.empty() && change <= side_tolerance * std::max(1.0, complexVectorMaxAbs(current))) {
            if(++stable == 2) {
               result.reciprocal_extent = extent;
               result.reciprocal_residual = change;
               return;
            }
         } else stable = 0;
         previous = std::move(current);
      }
      throw std::runtime_error("reciprocal Ewald alias shells did not converge");
   };
   converge_real(analyticSeed(geometry.H, alpha, side_tolerance));
   converge_reciprocal(analyticSeed(aliasReciprocalMatrix(geometry), 1.0 / (2.0 * alpha), side_tolerance));
   return result;
}

} // namespace

double periodicKernelReciprocityError(const std::vector<double>& kernel, const DipolePeriodicGeometry& geometry) {
   const std::size_t grid_cells = cells(geometry);
   const std::size_t expected = grid_cells * 9 * static_cast<std::size_t>(geometry.basis) * geometry.basis;
   if(kernel.size() != expected) throw std::invalid_argument("kernel size does not match periodic geometry");
   double maximum = 0.0;
   for(std::size_t d3 = 0; d3 < geometry.grid[2]; ++d3) for(std::size_t d2 = 0; d2 < geometry.grid[1]; ++d2)
      for(std::size_t d1 = 0; d1 < geometry.grid[0]; ++d1) {
         const std::size_t cell = d1 + geometry.grid[0] * (d2 + geometry.grid[1] * d3);
         const std::size_t r1 = (geometry.grid[0] - d1) % geometry.grid[0];
         const std::size_t r2 = (geometry.grid[1] - d2) % geometry.grid[1];
         const std::size_t r3 = (geometry.grid[2] - d3) % geometry.grid[2];
         const std::size_t reverse = r1 + geometry.grid[0] * (r2 + geometry.grid[1] * r3);
         for(unsigned int target = 0; target < geometry.basis; ++target)
            for(unsigned int source = 0; source < geometry.basis; ++source)
               for(unsigned int row = 0; row < 3; ++row) for(unsigned int column = 0; column < 3; ++column) {
                  const std::size_t batch = row + 3 * (column + 3 * (target + geometry.basis * source));
                  const std::size_t reciprocal_batch = column + 3 * (row + 3 * (source + geometry.basis * target));
                  maximum = std::max(maximum, std::abs(kernel[cell + grid_cells * batch] -
                                                       kernel[reverse + grid_cells * reciprocal_batch]));
               }
      }
   return maximum;
}

DipoleKernelBuildResult buildPeriodicEwaldDisplacementKernel(const DipolePeriodicGeometry& geometry,
                                                              const DipoleKernelSettings& settings) {
   validate(geometry, settings);
   const auto started = std::chrono::steady_clock::now();
   const double length_scale = std::cbrt(geometry.volume);
   const double nominal_alpha = std::sqrt(pi) / length_scale;
   const std::array<double, 3> alpha_scale{{0.75, 1.0, 1.25}};
   std::vector<CandidateResult> candidates;
   candidates.reserve(settings.explicit_alpha_for_testing > 0.0 ? 1 : alpha_scale.size());
   if(settings.explicit_alpha_for_testing > 0.0) {
      candidates.push_back(convergeCandidate(geometry, settings.explicit_alpha_for_testing, settings.tolerance));
   } else {
      for(const double scale : alpha_scale)
         candidates.push_back(convergeCandidate(geometry, nominal_alpha * scale, settings.tolerance));
   }
   const auto chosen = std::min_element(candidates.begin(), candidates.end(), [](const CandidateResult& left,
                                                                                   const CandidateResult& right) {
      return left.real_work + left.reciprocal_work < right.real_work + right.reciprocal_work;
   });
   std::size_t complete_real_work = 0, complete_reciprocal_work = 0;
   DipoleKernelBuildResult result{};
   result.kernel = completeKernel(geometry, *chosen, complete_real_work, complete_reciprocal_work);
   std::size_t check_real_work = 0, check_reciprocal_work = 0;
   if(candidates.size() > 1) {
      const auto alternate = candidates.begin() == chosen ? candidates.begin() + 1 : candidates.begin();
      const auto second_alpha_kernel = completeKernel(geometry, *alternate, check_real_work, check_reciprocal_work);
      result.diagnostics.max_alpha_difference = vectorMaxDifference(result.kernel, second_alpha_kernel);
      const double scale = std::max(1.0, vectorMaxAbs(result.kernel));
      if(result.diagnostics.max_alpha_difference > settings.tolerance * scale / 4.0) {
         throw std::runtime_error("independent alpha Ewald construction check did not meet gpu_dipole_tol");
      }
   }
   result.diagnostics.selected.alpha = chosen->alpha;
   result.diagnostics.selected.real_images = {chosen->real_extent, chosen->real_extent, chosen->real_extent};
   result.diagnostics.selected.reciprocal_images = {chosen->reciprocal_extent, chosen->reciprocal_extent,
                                                      chosen->reciprocal_extent};
   result.diagnostics.selected.tolerance = settings.tolerance;
   result.diagnostics.real_tail_residual = chosen->real_residual;
   result.diagnostics.reciprocal_tail_residual = chosen->reciprocal_residual;
   result.diagnostics.residual = std::max(chosen->real_residual, chosen->reciprocal_residual);
   result.diagnostics.max_reciprocity_error = periodicKernelReciprocityError(result.kernel, geometry);
   // In the real displacement representation the Cartesian/basis Hermitian
   // relation reduces to the same target/source reciprocity relation.  The
   // transformed spectrum is checked again at the upload boundary.
   result.diagnostics.max_hermitian_error = result.diagnostics.max_reciprocity_error;
   result.diagnostics.reciprocal_identity_error = reciprocalIdentityError(geometry);
   result.diagnostics.real_tensor_evaluations = complete_real_work + check_real_work;
   result.diagnostics.reciprocal_tensor_evaluations = complete_reciprocal_work + check_reciprocal_work;
   for(const auto& candidate : candidates) {
      result.diagnostics.real_tensor_evaluations += candidate.real_work;
      result.diagnostics.reciprocal_tensor_evaluations += candidate.reciprocal_work;
   }
   result.diagnostics.setup_work = result.diagnostics.real_tensor_evaluations +
                                   result.diagnostics.reciprocal_tensor_evaluations;
   result.diagnostics.setup_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
   return result;
}

DipoleAliasSpectrumBuildResult buildPeriodicEwaldAliasSpectrum(const DipolePeriodicGeometry& geometry,
                                                                 const DipoleKernelSettings& settings) {
   validate(geometry, settings);
   const auto started = std::chrono::steady_clock::now();
   const double nominal_alpha = std::sqrt(pi) / std::cbrt(geometry.volume);
   const std::array<double, 3> alpha_scale{{0.75, 1.0, 1.25}};
   std::vector<CandidateResult> candidates;
   candidates.reserve(settings.explicit_alpha_for_testing > 0.0 ? 1 : alpha_scale.size());
   if(settings.explicit_alpha_for_testing > 0.0) {
      candidates.push_back(convergeAliasCandidate(geometry, settings.explicit_alpha_for_testing, settings.tolerance));
   } else {
      for(const double scale : alpha_scale)
         candidates.push_back(convergeAliasCandidate(geometry, nominal_alpha * scale, settings.tolerance));
   }
   const auto chosen = std::min_element(candidates.begin(), candidates.end(), [](const CandidateResult& left,
                                                                                   const CandidateResult& right) {
      return left.real_work + left.reciprocal_work < right.real_work + right.reciprocal_work;
   });
   std::size_t complete_real_work = 0;
   DipoleAliasSpectrumBuildResult result{};
   result.real_kernel = realSelfKernel(geometry, *chosen, complete_real_work);
   std::size_t complete_reciprocal_work = 0;
   result.reciprocal_alias_spectrum = sumReciprocalAliases(geometry, chosen->alpha, chosen->reciprocal_extent,
                                                            complete_reciprocal_work);
   result.spectral_grid = {geometry.grid[0] / 2 + 1, geometry.grid[1], geometry.grid[2]};
   result.diagnostics.selected.alpha = chosen->alpha;
   result.diagnostics.selected.real_images = {chosen->real_extent, chosen->real_extent, chosen->real_extent};
   result.diagnostics.selected.reciprocal_images = {chosen->reciprocal_extent, chosen->reciprocal_extent,
                                                      chosen->reciprocal_extent};
   result.diagnostics.selected.tolerance = settings.tolerance;
   result.diagnostics.real_tail_residual = chosen->real_residual;
   result.diagnostics.reciprocal_tail_residual = chosen->reciprocal_residual;
   result.diagnostics.residual = std::max(chosen->real_residual, chosen->reciprocal_residual);
   result.diagnostics.max_reciprocity_error = periodicKernelReciprocityError(result.real_kernel, geometry);
   result.diagnostics.max_hermitian_error = result.diagnostics.max_reciprocity_error;
   result.diagnostics.reciprocal_identity_error = reciprocalIdentityError(geometry);
   result.diagnostics.real_tensor_evaluations = complete_real_work;
   result.diagnostics.reciprocal_tensor_evaluations = complete_reciprocal_work;
   for(const auto& candidate : candidates) {
      result.diagnostics.real_tensor_evaluations += candidate.real_work;
      result.diagnostics.reciprocal_tensor_evaluations += candidate.reciprocal_work;
   }
   result.diagnostics.setup_work = result.diagnostics.real_tensor_evaluations +
                                   result.diagnostics.reciprocal_tensor_evaluations;
   result.diagnostics.setup_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
   return result;
}

std::vector<std::complex<double>> referenceNormalizedKernelSpectrum(const std::vector<double>& kernel,
                                                                      const DipolePeriodicGeometry& geometry) {
   validate(geometry, {});
   const std::size_t grid_cells = cells(geometry);
   const std::size_t batches = 9 * static_cast<std::size_t>(geometry.basis) * geometry.basis;
   if(kernel.size() != grid_cells * batches) throw std::invalid_argument("kernel size does not match periodic geometry");
   const std::size_t stored_n1 = geometry.grid[0] / 2 + 1;
   std::vector<Complex> spectrum(spectralCells(geometry) * batches, Complex{});
   const double inverse_cells = 1.0 / static_cast<double>(grid_cells);
   for(std::size_t q3 = 0; q3 < geometry.grid[2]; ++q3)
      for(std::size_t q2 = 0; q2 < geometry.grid[1]; ++q2)
         for(std::size_t q1 = 0; q1 < stored_n1; ++q1) {
            const std::size_t spectral = q1 + stored_n1 * (q2 + geometry.grid[1] * q3);
            for(std::size_t d3 = 0; d3 < geometry.grid[2]; ++d3)
               for(std::size_t d2 = 0; d2 < geometry.grid[1]; ++d2)
                  for(std::size_t d1 = 0; d1 < geometry.grid[0]; ++d1) {
                     const std::size_t cell = d1 + geometry.grid[0] * (d2 + geometry.grid[1] * d3);
                     const double angle = -2.0 * pi * (static_cast<double>(q1 * d1) / geometry.grid[0] +
                        static_cast<double>(q2 * d2) / geometry.grid[1] +
                        static_cast<double>(q3 * d3) / geometry.grid[2]);
                     const Complex phase = std::polar(inverse_cells, angle);
                     for(std::size_t batch = 0; batch < batches; ++batch)
                        spectrum[spectral + spectralCells(geometry) * batch] += kernel[cell + grid_cells * batch] * phase;
                  }
         }
   return spectrum;
}
