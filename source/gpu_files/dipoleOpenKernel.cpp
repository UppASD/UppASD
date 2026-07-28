#include "dipoleOpenKernel.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using Vector = std::array<double, 3>;
using Tensor = std::array<double, 9>;

struct SignedMagnitude {
   std::size_t magnitude = 0;
   bool negative = false;
};

std::string gridText(const std::array<std::size_t, 3>& grid) {
   std::ostringstream stream;
   stream << '(' << grid[0] << ',' << grid[1] << ',' << grid[2] << ')';
   return stream.str();
}

std::size_t checkedMultiply(std::size_t left, std::size_t right, const std::string& field) {
   if(left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
      std::ostringstream stream;
      stream << field << " overflows size_t: " << left << " * " << right;
      throw std::overflow_error(stream.str());
   }
   return left * right;
}

std::size_t checkedCells(const std::array<std::size_t, 3>& grid, const char* field) {
   for(unsigned int axis = 0; axis < 3; ++axis) {
      if(grid[axis] == 0) {
         std::ostringstream stream;
         stream << field << '[' << axis << "] must be positive, got 0 in " << gridText(grid);
         throw std::invalid_argument(stream.str());
      }
   }
   return checkedMultiply(checkedMultiply(grid[0], grid[1], std::string(field) + " cell count"),
                          grid[2], std::string(field) + " cell count");
}

std::size_t minimumFftExtent(std::size_t active, unsigned int axis) {
   constexpr std::size_t largest_safe_active = std::numeric_limits<std::size_t>::max() / 2 + 1;
   if(active > largest_safe_active) {
      std::ostringstream stream;
      stream << "active_grid[" << axis << "]=" << active
             << " makes 2*active_grid-1 unrepresentable";
      throw std::overflow_error(stream.str());
   }
   return 2 * active - 1;
}

double determinant(const std::array<double, 9>& matrix) {
   return matrix[0] * (matrix[4] * matrix[8] - matrix[7] * matrix[5]) -
          matrix[3] * (matrix[1] * matrix[8] - matrix[7] * matrix[2]) +
          matrix[6] * (matrix[1] * matrix[5] - matrix[4] * matrix[2]);
}

std::size_t fftCell(const std::array<std::size_t, 3>& grid,
                    std::size_t q1, std::size_t q2, std::size_t q3) {
   return q1 + grid[0] * (q2 + grid[1] * q3);
}

std::size_t kernelBatch(unsigned int row, unsigned int column,
                        unsigned int target, unsigned int source,
                        unsigned int basis) {
   return static_cast<std::size_t>(row) +
          3 * (static_cast<std::size_t>(column) +
          3 * (static_cast<std::size_t>(target) +
          static_cast<std::size_t>(basis) * static_cast<std::size_t>(source)));
}

std::size_t encodedDisplacements(std::size_t active, unsigned int axis) {
   return minimumFftExtent(active, axis);
}

SignedMagnitude decodeDisplacement(std::size_t encoded, std::size_t active) {
   if(encoded < active) return {encoded, false};
   return {encoded - active + 1, true};
}

std::size_t embeddedCoordinate(const SignedMagnitude& displacement, std::size_t fft_extent) {
   return displacement.negative ? fft_extent - displacement.magnitude : displacement.magnitude;
}

SignedMagnitude opposite(const SignedMagnitude& displacement) {
   if(displacement.magnitude == 0) return displacement;
   return {displacement.magnitude, !displacement.negative};
}

Vector displacementVector(const DipoleOpenGeometry& geometry,
                          const std::array<SignedMagnitude, 3>& displacement,
                          const std::array<std::size_t, 3>& target_offset,
                          const std::array<std::size_t, 3>& source_offset,
                          unsigned int target, unsigned int source) {
   Vector result{};
   for(unsigned int component = 0; component < 3; ++component) {
      double lattice = 0.0;
      for(unsigned int axis = 0; axis < 3; ++axis) {
         double coarse =
            static_cast<double>(displacement[axis].magnitude) *
            static_cast<double>(geometry.block_shape[axis]);
         if(displacement[axis].negative) coarse = -coarse;
         const double intra_block =
            static_cast<double>(target_offset[axis]) -
            static_cast<double>(source_offset[axis]);
         lattice += (coarse + intra_block) *
                    geometry.primitive_vectors[3 * axis + component];
      }
      const double offset = geometry.basis_offsets[target][component] -
                            geometry.basis_offsets[source][component];
      result[component] = lattice + offset;
      if(!std::isfinite(result[component])) {
         std::ostringstream stream;
         stream << "primitive_vectors and basis_offsets produce a non-finite displacement"
                << " for component " << component << ", target basis " << target
                << ", source basis " << source;
         throw std::invalid_argument(stream.str());
      }
   }
   return result;
}

bool exactZero(const Vector& vector) {
   return vector[0] == 0.0 && vector[1] == 0.0 && vector[2] == 0.0;
}

bool zeroDisplacement(const std::array<SignedMagnitude, 3>& displacement) {
   return displacement[0].magnitude == 0 &&
          displacement[1].magnitude == 0 &&
          displacement[2].magnitude == 0;
}

bool sameOffset(const std::array<std::size_t, 3>& left,
                const std::array<std::size_t, 3>& right) {
   return left == right;
}

void validateModel(const DipoleOpenGeometry& geometry,
                   std::size_t& atomistic_cells, std::size_t& active_cells,
                   std::size_t& fft_cells, std::size_t& block_population,
                   std::size_t& kernel_batches, std::size_t& kernel_elements) {
   atomistic_cells = checkedCells(geometry.atomistic_grid, "atomistic_grid");
   active_cells = checkedCells(geometry.active_grid, "active_grid");
   fft_cells = checkedCells(geometry.fft_grid, "fft_grid");

   for(unsigned int axis = 0; axis < 3; ++axis) {
      const std::size_t minimum = minimumFftExtent(geometry.active_grid[axis], axis);
      if(geometry.fft_grid[axis] < minimum) {
         std::ostringstream stream;
         stream << "fft_grid[" << axis << "]=" << geometry.fft_grid[axis]
                << " is smaller than 2*active_grid[" << axis << "]-1=" << minimum
                << " for active_grid=" << gridText(geometry.active_grid)
                << " and fft_grid=" << gridText(geometry.fft_grid);
         throw std::invalid_argument(stream.str());
      }
      if(geometry.block_shape[axis] == 0) {
         std::ostringstream stream;
         stream << "block_shape[" << axis << "]=" << geometry.block_shape[axis]
                << " must be positive";
         throw std::invalid_argument(stream.str());
      }
      if(geometry.atomistic_grid[axis] % geometry.block_shape[axis] != 0) {
         std::ostringstream stream;
         stream << "atomistic_grid[" << axis << "]=" << geometry.atomistic_grid[axis]
                << " is not divisible by block_shape[" << axis << "]="
                << geometry.block_shape[axis] << "; partial edge blocks are unsupported";
         throw std::invalid_argument(stream.str());
      }
      const std::size_t projected =
         geometry.atomistic_grid[axis] / geometry.block_shape[axis];
      if(projected != geometry.active_grid[axis]) {
         std::ostringstream stream;
         stream << "active_grid[" << axis << "]=" << geometry.active_grid[axis]
                << " does not equal atomistic_grid[" << axis << "]/block_shape[" << axis
                << "]=" << projected << "; every active macro channel must have the same full block shape";
         throw std::invalid_argument(stream.str());
      }
   }
   block_population = checkedMultiply(
      checkedMultiply(geometry.block_shape[0], geometry.block_shape[1],
                      "block_population=block_shape[0]*block_shape[1]*block_shape[2]"),
      geometry.block_shape[2],
      "block_population=block_shape[0]*block_shape[1]*block_shape[2]");
   (void)checkedMultiply(active_cells, static_cast<std::size_t>(geometry.basis),
                         "active macro channel count=active_cells*basis");

   if(geometry.basis == 0) {
      throw std::invalid_argument("basis must be positive, got 0");
   }
   if(geometry.basis_offsets.size() != geometry.basis) {
      std::ostringstream stream;
      stream << "basis_offsets size " << geometry.basis_offsets.size()
             << " does not match basis=" << geometry.basis;
      throw std::invalid_argument(stream.str());
   }
   for(std::size_t index = 0; index < geometry.primitive_vectors.size(); ++index) {
      if(!std::isfinite(geometry.primitive_vectors[index])) {
         std::ostringstream stream;
         stream << "primitive_vectors[" << index << "] must be finite, got "
                << geometry.primitive_vectors[index];
         throw std::invalid_argument(stream.str());
      }
   }
   for(std::size_t basis = 0; basis < geometry.basis_offsets.size(); ++basis)
      for(unsigned int component = 0; component < 3; ++component)
         if(!std::isfinite(geometry.basis_offsets[basis][component])) {
            std::ostringstream stream;
            stream << "basis_offsets[" << basis << "][" << component
                   << "] must be finite, got " << geometry.basis_offsets[basis][component];
            throw std::invalid_argument(stream.str());
         }

   const double det = determinant(geometry.primitive_vectors);
   if(!std::isfinite(det) || det <= 0.0) {
      std::ostringstream stream;
      stream << "primitive_vectors determinant must be finite and positive, got " << det;
      throw std::invalid_argument(stream.str());
   }

   const std::size_t basis_size = static_cast<std::size_t>(geometry.basis);
   kernel_batches = checkedMultiply(
      checkedMultiply(9, basis_size, "kernel_batches=9*basis*basis"),
      basis_size, "kernel_batches=9*basis*basis");
   kernel_elements = checkedMultiply(fft_cells, kernel_batches,
                                     "kernel element count=fft_cells*kernel_batches");
   const std::size_t pair_count = checkedMultiply(
      block_population, block_population,
      "block pair count=block_population*block_population");
   const std::size_t displacement_count = checkedMultiply(
      checkedMultiply(encodedDisplacements(geometry.active_grid[0], 0),
                      encodedDisplacements(geometry.active_grid[1], 1),
                      "projected displacement count"),
      encodedDisplacements(geometry.active_grid[2], 2),
      "projected displacement count");
   const std::size_t pair_component_work = checkedMultiply(
      pair_count, kernel_batches, "projected block pair/component work");
   (void)checkedMultiply(displacement_count, pair_component_work,
                         "total projected tensor construction work");
   (void)checkedMultiply(kernel_elements, sizeof(double), "kernel byte count");
   if(kernel_elements > std::vector<double>{}.max_size()) {
      std::ostringstream stream;
      stream << "kernel element count " << kernel_elements
             << " exceeds std::vector<double>::max_size()";
      throw std::overflow_error(stream.str());
   }
}

bool axisInGap(std::size_t q, std::size_t active, std::size_t fft) {
   return q >= active && q <= fft - active;
}

} // namespace

DipoleOpenKernelResult buildOpenDipoleDisplacementKernel(
   const DipoleOpenGeometry& geometry) {
   std::size_t atomistic_cells = 0;
   std::size_t active_cells = 0;
   std::size_t fft_cells = 0;
   std::size_t block_population = 0;
   std::size_t kernel_batches = 0;
   std::size_t kernel_elements = 0;
   validateModel(geometry, atomistic_cells, active_cells, fft_cells,
                 block_population, kernel_batches, kernel_elements);

   DipoleOpenKernelResult result{};
   result.kernel.assign(kernel_elements, 0.0);
   result.diagnostics.atomistic_grid = geometry.atomistic_grid;
   result.diagnostics.active_grid = geometry.active_grid;
   result.diagnostics.fft_grid = geometry.fft_grid;
   result.diagnostics.block_shape = geometry.block_shape;
   result.diagnostics.block_population = block_population;
   result.diagnostics.atomistic_cells = atomistic_cells;
   result.diagnostics.active_cells = active_cells;
   result.diagnostics.fft_cells = fft_cells;
   result.diagnostics.kernel_batches = kernel_batches;

   double minimum_nonself_r2 = std::numeric_limits<double>::infinity();
   const std::array<std::size_t, 3> displacement_counts{{
      encodedDisplacements(geometry.active_grid[0], 0),
      encodedDisplacements(geometry.active_grid[1], 1),
      encodedDisplacements(geometry.active_grid[2], 2)}};

   for(std::size_t e3 = 0; e3 < displacement_counts[2]; ++e3)
      for(std::size_t e2 = 0; e2 < displacement_counts[1]; ++e2)
         for(std::size_t e1 = 0; e1 < displacement_counts[0]; ++e1) {
            const std::array<SignedMagnitude, 3> displacement{{
               decodeDisplacement(e1, geometry.active_grid[0]),
               decodeDisplacement(e2, geometry.active_grid[1]),
               decodeDisplacement(e3, geometry.active_grid[2])}};
            const std::size_t q1 = embeddedCoordinate(displacement[0], geometry.fft_grid[0]);
            const std::size_t q2 = embeddedCoordinate(displacement[1], geometry.fft_grid[1]);
            const std::size_t q3 = embeddedCoordinate(displacement[2], geometry.fft_grid[2]);
            const std::size_t cell = fftCell(geometry.fft_grid, q1, q2, q3);

            for(unsigned int target = 0; target < geometry.basis; ++target)
               for(unsigned int source = 0; source < geometry.basis; ++source) {
                  Tensor projected{};
                  for(std::size_t u3 = 0; u3 < geometry.block_shape[2]; ++u3)
                     for(std::size_t u2 = 0; u2 < geometry.block_shape[1]; ++u2)
                        for(std::size_t u1 = 0; u1 < geometry.block_shape[0]; ++u1)
                           for(std::size_t v3 = 0; v3 < geometry.block_shape[2]; ++v3)
                              for(std::size_t v2 = 0; v2 < geometry.block_shape[1]; ++v2)
                                 for(std::size_t v1 = 0; v1 < geometry.block_shape[0]; ++v1) {
                                    const std::array<std::size_t, 3> target_offset{{u1, u2, u3}};
                                    const std::array<std::size_t, 3> source_offset{{v1, v2, v3}};
                                    if(zeroDisplacement(displacement) && target == source &&
                                       sameOffset(target_offset, source_offset)) {
                                       // K_point(i,i)=0 exactly.  Other pairs
                                       // at d=0 form the finite coarse diagonal.
                                       continue;
                                    }

                                    const Vector r = displacementVector(
                                       geometry, displacement, target_offset, source_offset,
                                       target, source);
                                    if(exactZero(r)) {
                                       std::ostringstream stream;
                                       stream << "basis_offsets, primitive_vectors, and block_shape overlap"
                                              << " at a nonself point pair for target basis " << target
                                              << ", source basis " << source
                                              << ", embedded cell " << cell;
                                       throw std::invalid_argument(stream.str());
                                    }
                                    const double radius = std::hypot(r[0], r[1], r[2]);
                                    const double r2 = radius * radius;
                                    if(!std::isfinite(radius) || !std::isfinite(r2) || r2 <= 0.0) {
                                       std::ostringstream stream;
                                       stream << "minimum_nonself_r2 is unrepresentable for target basis "
                                              << target << ", source basis " << source
                                              << ", radius=" << radius << ", r2=" << r2;
                                       throw std::invalid_argument(stream.str());
                                    }
                                    minimum_nonself_r2 = std::min(minimum_nonself_r2, r2);

                                    const double inverse_radius = 1.0 / radius;
                                    const double inverse_r3 =
                                       inverse_radius * inverse_radius * inverse_radius;
                                    if(!std::isfinite(inverse_r3)) {
                                       std::ostringstream stream;
                                       stream << "primitive_vectors and basis_offsets produce a non-finite"
                                              << " 1/|r|^3 for target basis " << target
                                              << ", source basis " << source << ", r2=" << r2;
                                       throw std::invalid_argument(stream.str());
                                    }
                                    const Vector direction{{r[0] * inverse_radius,
                                                            r[1] * inverse_radius,
                                                            r[2] * inverse_radius}};
                                    for(unsigned int column = 0; column < 3; ++column)
                                       for(unsigned int row = 0; row < 3; ++row) {
                                          const double value =
                                             (3.0 * direction[row] * direction[column] -
                                              (row == column ? 1.0 : 0.0)) * inverse_r3;
                                          if(!std::isfinite(value)) {
                                             std::ostringstream stream;
                                             stream << "open dipole tensor contains a non-finite value"
                                                    << " for row " << row << ", column " << column
                                                    << ", target basis " << target
                                                    << ", source basis " << source << ", r2=" << r2;
                                             throw std::invalid_argument(stream.str());
                                          }
                                          projected[row + 3 * column] += value;
                                       }
                                 }
                  const double inverse_population_squared =
                     1.0 / (static_cast<double>(block_population) *
                            static_cast<double>(block_population));
                  for(unsigned int column = 0; column < 3; ++column)
                     for(unsigned int row = 0; row < 3; ++row) {
                        const std::size_t batch =
                           kernelBatch(row, column, target, source, geometry.basis);
                        result.kernel[cell + fft_cells * batch] =
                           projected[row + 3 * column] * inverse_population_squared;
                     }
               }
         }

   result.diagnostics.minimum_nonself_r2 =
      std::isfinite(minimum_nonself_r2) ? minimum_nonself_r2 : 0.0;

   double maximum_kernel_abs = 0.0;
   // Independently audit finite values and the explicitly unused padding gap.
   for(std::size_t q3 = 0; q3 < geometry.fft_grid[2]; ++q3)
      for(std::size_t q2 = 0; q2 < geometry.fft_grid[1]; ++q2)
         for(std::size_t q1 = 0; q1 < geometry.fft_grid[0]; ++q1) {
            const std::size_t cell = fftCell(geometry.fft_grid, q1, q2, q3);
            const bool gap = axisInGap(q1, geometry.active_grid[0], geometry.fft_grid[0]) ||
                             axisInGap(q2, geometry.active_grid[1], geometry.fft_grid[1]) ||
                             axisInGap(q3, geometry.active_grid[2], geometry.fft_grid[2]);
            for(std::size_t batch = 0; batch < kernel_batches; ++batch) {
               const double value = result.kernel[cell + fft_cells * batch];
               if(!std::isfinite(value)) ++result.diagnostics.nonfinite_values;
               maximum_kernel_abs = std::max(maximum_kernel_abs, std::abs(value));
               if(gap) result.diagnostics.max_padding_gap_abs =
                  std::max(result.diagnostics.max_padding_gap_abs, std::abs(value));
            }
         }
   result.diagnostics.all_finite = result.diagnostics.nonfinite_values == 0;

   // A projected diagonal is not point self: it contains every distinct
   // intra-block pair.  max_point_self_abs remains exact zero because those
   // terms were skipped before accumulation; expose the finite diagonal
   // separately so a later caller cannot mistake it for an ad hoc self term.
   for(unsigned int basis = 0; basis < geometry.basis; ++basis)
      for(unsigned int column = 0; column < 3; ++column)
         for(unsigned int row = 0; row < 3; ++row) {
            const std::size_t batch =
               kernelBatch(row, column, basis, basis, geometry.basis);
            result.diagnostics.max_projected_diagonal_abs = std::max(
               result.diagnostics.max_projected_diagonal_abs,
               std::abs(result.kernel[fft_cells * batch]));
         }

   // K(d,a,b) = K(-d,b,a)^T in the exact padded storage convention.
   for(std::size_t e3 = 0; e3 < displacement_counts[2]; ++e3)
      for(std::size_t e2 = 0; e2 < displacement_counts[1]; ++e2)
         for(std::size_t e1 = 0; e1 < displacement_counts[0]; ++e1) {
            const std::array<SignedMagnitude, 3> displacement{{
               decodeDisplacement(e1, geometry.active_grid[0]),
               decodeDisplacement(e2, geometry.active_grid[1]),
               decodeDisplacement(e3, geometry.active_grid[2])}};
            const std::array<SignedMagnitude, 3> reverse{{
               opposite(displacement[0]), opposite(displacement[1]), opposite(displacement[2])}};
            const std::size_t cell = fftCell(
               geometry.fft_grid,
               embeddedCoordinate(displacement[0], geometry.fft_grid[0]),
               embeddedCoordinate(displacement[1], geometry.fft_grid[1]),
               embeddedCoordinate(displacement[2], geometry.fft_grid[2]));
            const std::size_t reverse_cell = fftCell(
               geometry.fft_grid,
               embeddedCoordinate(reverse[0], geometry.fft_grid[0]),
               embeddedCoordinate(reverse[1], geometry.fft_grid[1]),
               embeddedCoordinate(reverse[2], geometry.fft_grid[2]));
            for(unsigned int target = 0; target < geometry.basis; ++target)
               for(unsigned int source = 0; source < geometry.basis; ++source)
                  for(unsigned int column = 0; column < 3; ++column)
                     for(unsigned int row = 0; row < 3; ++row) {
                        const std::size_t batch =
                           kernelBatch(row, column, target, source, geometry.basis);
                        const std::size_t reverse_batch =
                           kernelBatch(column, row, source, target, geometry.basis);
                        result.diagnostics.max_reciprocity_error = std::max(
                           result.diagnostics.max_reciprocity_error,
                           std::abs(result.kernel[cell + fft_cells * batch] -
                                    result.kernel[reverse_cell + fft_cells * reverse_batch]));
                     }
         }

   if(!result.diagnostics.all_finite) {
      throw std::logic_error("open dipole kernel finite-value audit failed after construction");
   }
   if(result.diagnostics.max_point_self_abs != 0.0) {
      throw std::logic_error("open dipole kernel exact-point-self audit is nonzero");
   }
   if(result.diagnostics.max_padding_gap_abs != 0.0) {
      throw std::logic_error("open dipole kernel unused-gap audit is nonzero");
   }
   const double reciprocity_headroom =
      256.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, maximum_kernel_abs);
   if(result.diagnostics.max_reciprocity_error > reciprocity_headroom) {
      std::ostringstream stream;
      stream << "open dipole kernel reciprocity audit " << result.diagnostics.max_reciprocity_error
             << " exceeds fp64 construction headroom " << reciprocity_headroom;
      throw std::logic_error(stream.str());
   }
   return result;
}
