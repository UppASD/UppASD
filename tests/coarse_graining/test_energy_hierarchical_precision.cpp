// CGP-00 ENERGY-HIERARCHICAL-PRECISION fixture.
//
// docs/CGP_work.md CGP-00 asks whether the RCG-09B hierarchical energy
// reduction (per-thread local contribution -> block-local shared-memory tree
// reduction -> ordered serial sum of block partials) genuinely needs FP64
// accumulation, now that its topology is nothing like the old globally
// contended atomicAdd accumulator test_energy_fp32_accum.cpp (RCG-06B,
// F-11) isolates. That earlier fixture is deliberately left untouched: its
// own header comment protects it as historical negative-control evidence for
// a different, already-fixed defect class, and CGP-00 explicitly says not to
// compare it against the new topology and call that an FP32-vs-FP64 result.
//
// This fixture instead reproduces the *current* topology exactly, sharing
// the literal block-tree primitive production uses
// (reduceAdaptiveEnergyBlockT<Acc> in gpuAtomicDouble.hpp, generalized by
// CGP-00 from gpuAdaptiveRuntime.cpp's RCG-09B reduceAdaptiveEnergyBlock so
// this fixture cannot silently drift from what production actually runs),
// and sweeps its accumulator precision at three stages independently:
//
//   local  -- per-thread contribution precision before it enters the tree;
//   block  -- shared-memory tree-reduction accumulator precision;
//   final  -- ordered serial sum of one partial per launch block.
//
// Four modes are compared (CGP-00 Part A):
//   1. F32 local / F32 block / F32 final
//   2. F32 local / F32 block / F64 final
//   3. F32 local / F64 block / F64 final  (today's SINGLE_PREC GPU build)
//   4. F64 local / F64 block / F64 final  (today's DOUBLE_PREC GPU build,
//      and the numerical reference this fixture measures error against)
//
// against a host-computed Neumaier-compensated `long double` reference sum,
// over five term distributions (CGP-00 Part B: all-positive, alternating
// sign, severe cancellation, wide dynamic range, randomized
// physically-plausible) plus a sixth case isolating the accuracy of a small
// difference between two large sums, which is the quantity a resolution-
// transition accept/reject decision actually depends on (CGP-00 Part D).
//
// Unlike test_energy_fp32_accum.cpp's atomicAdd accumulator, this reduction
// has a fixed, hardware-scheduling-independent accumulation order (grid
// size and block assignment are deterministic; the tree reduction and the
// host-side ordered final sum contain no data race), so a single run per
// (case, N, mode) is exact and reproducible -- no repeats are needed to
// average out nondeterminism.

#include "gpuAtomicDouble.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#if defined(CUDA_V)
#include <cuda_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != cudaSuccess) { \
   std::fprintf(stderr, "%s\n", cudaGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE cudaErrorNoDevice
#define GPU_GET_COUNT cudaGetDeviceCount
#define GPU_MALLOC cudaMalloc
#define GPU_MEMCPY cudaMemcpy
#define GPU_MEMCPY_D2H cudaMemcpyDeviceToHost
#define GPU_MEMCPY_H2D cudaMemcpyHostToDevice
#define GPU_FREE cudaFree
#define GPU_DEVICE_SYNC cudaDeviceSynchronize
#define GPU_GET_LAST_ERROR cudaGetLastError
#elif defined(HIP_V)
#include <hip/hip_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != hipSuccess) { \
   std::fprintf(stderr, "%s\n", hipGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE hipErrorNoDevice
#define GPU_GET_COUNT hipGetDeviceCount
#define GPU_MALLOC hipMalloc
#define GPU_MEMCPY hipMemcpy
#define GPU_MEMCPY_D2H hipMemcpyDeviceToHost
#define GPU_MEMCPY_H2D hipMemcpyHostToDevice
#define GPU_FREE hipFree
#define GPU_DEVICE_SYNC hipDeviceSynchronize
#define GPU_GET_LAST_ERROR hipGetLastError
#else
#error "This fixture requires the configured CUDA or HIP backend"
#endif

namespace {

constexpr unsigned int threadsPerBlock = 256;
constexpr std::uint64_t maxGrid = 195312; // same cap as test_energy_fp32_accum.cpp

// One thread per term. Each thread casts its term down to LocalT (modeling
// the precision the field/operator math that produced it actually ran in),
// casts that up/down to BlockT, and enters the identical block-tree
// primitive production's evaluateAdaptive* kernels use.
template <typename LocalT, typename BlockT>
__global__ void hierarchicalLocalBlockKernel(
   const double* terms, std::uint64_t n, BlockT* partials) {
   extern __shared__ unsigned char rawShared[];
   BlockT* shared = reinterpret_cast<BlockT*>(rawShared);
   const std::uint64_t index =
      static_cast<std::uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   const LocalT localValue =
      index < n ? static_cast<LocalT>(terms[index]) : LocalT(0);
   reduceAdaptiveEnergyBlockT<BlockT>(
      static_cast<BlockT>(localValue), partials, shared);
}

unsigned int gridBlocksFor(std::uint64_t n) {
   const std::uint64_t blocks =
      (n + threadsPerBlock - 1) / threadsPerBlock;
   return static_cast<unsigned int>(std::min<std::uint64_t>(maxGrid,
      std::max<std::uint64_t>(1, blocks)));
}

// Runs one (LocalT, BlockT) instantiation on device, then combines the
// block partials in ascending block order using FinalT -- the same shape as
// production's reduceAdaptiveEnergyPartials (an ordered serial sum, not
// another reduction tree). Returns the result widened to double for
// reporting; FinalT precision is genuinely used for the running sum itself.
template <typename LocalT, typename BlockT, typename FinalT>
int runHierarchicalMode(const std::vector<double>& terms, double& outResult) {
   const std::uint64_t n = terms.size();
   const unsigned int gridBlocks = gridBlocksFor(n);

   double* deviceTerms = nullptr;
   BlockT* devicePartials = nullptr;
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceTerms),
                        n * sizeof(double)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&devicePartials),
                        static_cast<std::size_t>(gridBlocks) * sizeof(BlockT)));
   GPU_CHECK(GPU_MEMCPY(deviceTerms, terms.data(), n * sizeof(double),
                        GPU_MEMCPY_H2D));

   const std::size_t sharedBytes = threadsPerBlock * sizeof(BlockT);
   hierarchicalLocalBlockKernel<LocalT, BlockT><<<gridBlocks, threadsPerBlock,
      sharedBytes>>>(deviceTerms, n, devicePartials);
   GPU_CHECK(GPU_GET_LAST_ERROR());
   GPU_CHECK(GPU_DEVICE_SYNC());

   std::vector<BlockT> hostPartials(gridBlocks);
   GPU_CHECK(GPU_MEMCPY(hostPartials.data(), devicePartials,
                        static_cast<std::size_t>(gridBlocks) * sizeof(BlockT),
                        GPU_MEMCPY_D2H));
   GPU_CHECK(GPU_FREE(deviceTerms));
   GPU_CHECK(GPU_FREE(devicePartials));

   FinalT total = FinalT(0);
   for(unsigned int block = 0; block < gridBlocks; ++block)
      total += static_cast<FinalT>(hostPartials[block]);
   outResult = static_cast<double>(total);
   return 0;
}

// Neumaier-compensated summation in long double: the high-fidelity host
// reference every mode's error is measured against.
double neumairReference(const std::vector<double>& terms) {
   long double sum = 0.0L;
   long double compensation = 0.0L;
   for(const double term : terms) {
      const long double value = static_cast<long double>(term);
      const long double t = sum + value;
      if(std::fabs(sum) >= std::fabs(value))
         compensation += (sum - t) + value;
      else
         compensation += (value - t) + sum;
      sum = t;
   }
   return static_cast<double>(sum + compensation);
}

struct ModeResult {
   const char* label;
   double absError;
   double relError;
   double perTermError;
};

struct CaseOutcome {
   std::string name;
   std::uint64_t n;
   double reference;
   std::vector<ModeResult> modes;
};

int evaluateAllModes(const std::string& caseName, std::uint64_t n,
                     const std::vector<double>& terms, CaseOutcome& outcome) {
   outcome.name = caseName;
   outcome.n = n;
   outcome.reference = neumairReference(terms);

   struct ModeSpec {
      const char* label;
      int (*run)(const std::vector<double>&, double&);
   };
   const ModeSpec specs[4] = {
      {"F32/F32/F32", &runHierarchicalMode<float, float, float>},
      {"F32/F32/F64", &runHierarchicalMode<float, float, double>},
      {"F32/F64/F64", &runHierarchicalMode<float, double, double>},
      {"F64/F64/F64", &runHierarchicalMode<double, double, double>},
   };

   for(const auto& spec : specs) {
      double result = 0.0;
      if(spec.run(terms, result) != 0) return 1;
      ModeResult mode;
      mode.label = spec.label;
      mode.absError = std::fabs(result - outcome.reference);
      mode.relError = std::fabs(outcome.reference) > 0.0 ?
         mode.absError / std::fabs(outcome.reference) : 0.0;
      mode.perTermError = mode.absError / static_cast<double>(n);
      outcome.modes.push_back(mode);
   }
   return 0;
}

void printOutcome(const CaseOutcome& outcome) {
   std::printf("%-24s N=%9llu reference=%+.10e\n", outcome.name.c_str(),
               static_cast<unsigned long long>(outcome.n), outcome.reference);
   for(const auto& mode : outcome.modes)
      std::printf("   %-12s abs_err=%14.6e rel_err=%14.6e per_term_err=%14.6e\n",
                  mode.label, mode.absError, mode.relError, mode.perTermError);
}

// -- Term-distribution generators -------------------------------------

std::vector<double> caseAllPositive(std::uint64_t n) {
   std::vector<double> terms(n, 1.0 / static_cast<double>(n));
   return terms;
}

std::vector<double> caseAlternating(std::uint64_t n) {
   std::vector<double> terms(n);
   const double magnitude = 1.0 / static_cast<double>(n);
   for(std::uint64_t i = 0; i < n; ++i)
      terms[i] = (i % 2 == 0 ? magnitude : -magnitude);
   return terms;
}

// Large near-equal pairs (+C, -C) with a small, exactly-known per-pair
// residual riding on top -- the classic catastrophic-cancellation stress
// case: the terms are individually ~1e6x the quantity that actually matters.
std::vector<double> caseSevereCancellation(std::uint64_t n) {
   std::vector<double> terms(n);
   const double large = 1.0e6;
   const double residualStep = 1.0e-3 / static_cast<double>(n);
   for(std::uint64_t i = 0; i + 1 < n; i += 2) {
      const double residual = residualStep * static_cast<double>(i / 2);
      terms[i] = large + residual;
      terms[i + 1] = -large;
   }
   if(n % 2 == 1) terms[n - 1] = 0.0;
   return terms;
}

std::vector<double> caseWideDynamicRange(std::uint64_t n) {
   std::vector<double> terms(n);
   std::mt19937_64 rng(0xCC9000D1D2ULL);
   std::uniform_real_distribution<double> exponent(-6.0, 6.0);
   std::uniform_real_distribution<double> sign01(0.0, 1.0);
   for(std::uint64_t i = 0; i < n; ++i) {
      const double magnitude = std::pow(10.0, exponent(rng));
      terms[i] = (sign01(rng) < 0.5 ? -magnitude : magnitude);
   }
   return terms;
}

std::vector<double> caseRandomizedPhysical(std::uint64_t n) {
   std::vector<double> terms(n);
   std::mt19937_64 rng(0xCC9000A70FULL ^ n);
   // Order-1 Gaussian-ish magnitudes, representative of dimensionless
   // per-bond exchange-energy contributions in the adaptive-CG kernels
   // (kernels.exchangeStiffness * dot(gradient, gradient) etc.), not tied
   // to a specific physical unit system.
   std::normal_distribution<double> gaussian(0.0, 1.0);
   for(std::uint64_t i = 0; i < n; ++i)
      terms[i] = gaussian(rng);
   return terms;
}

// -- Small-difference-of-large-sums case (CGP-00 Part B + D relevance) --

// Builds two term sets of equal size whose exact sums differ by a small,
// known delta -- modeling energy_before/energy_after in a resolution
// transition, where the quantity that actually gates accept/reject is that
// small difference, not either total in isolation.
struct DifferencePair {
   std::vector<double> before;
   std::vector<double> after;
   double exactDelta;
};

DifferencePair buildDifferencePair(std::uint64_t n, double delta) {
   DifferencePair pair;
   pair.before.resize(n);
   pair.after.resize(n);
   std::mt19937_64 rng(0xCC9000DE17AULL ^ n);
   std::normal_distribution<double> gaussian(0.0, 10.0);
   for(std::uint64_t i = 0; i < n; ++i) {
      const double value = gaussian(rng);
      pair.before[i] = value;
      pair.after[i] = value;
   }
   // Perturb one term of "after" by exactly delta so the exact difference
   // between the two sums is delta regardless of N.
   pair.after[0] += delta;
   pair.exactDelta = delta;
   return pair;
}

int evaluateDifferenceCase(std::uint64_t n, double delta,
                           std::vector<CaseOutcome>& outcomes) {
   const DifferencePair pair = buildDifferencePair(n, delta);

   struct ModeSpec {
      const char* label;
      int (*run)(const std::vector<double>&, double&);
   };
   const ModeSpec specs[4] = {
      {"F32/F32/F32", &runHierarchicalMode<float, float, float>},
      {"F32/F32/F64", &runHierarchicalMode<float, float, double>},
      {"F32/F64/F64", &runHierarchicalMode<float, double, double>},
      {"F64/F64/F64", &runHierarchicalMode<double, double, double>},
   };

   CaseOutcome outcome;
   outcome.name = "energy_difference(delta=" + std::to_string(delta) + ")";
   outcome.n = n;
   outcome.reference = delta;

   for(const auto& spec : specs) {
      double before = 0.0, after = 0.0;
      if(spec.run(pair.before, before) != 0) return 1;
      if(spec.run(pair.after, after) != 0) return 1;
      const double computedDelta = after - before;
      ModeResult mode;
      mode.label = spec.label;
      mode.absError = std::fabs(computedDelta - delta);
      mode.relError = std::fabs(delta) > 0.0 ? mode.absError / std::fabs(delta) : 0.0;
      mode.perTermError = mode.absError / static_cast<double>(n);
      outcome.modes.push_back(mode);
   }
   outcomes.push_back(outcome);
   return 0;
}

} // namespace

int main() {
   int count = 0;
   const auto countStatus = GPU_GET_COUNT(&count);
   if(countStatus == GPU_NO_DEVICE || count == 0) {
      std::puts("ENERGY-HIERARCHICAL-PRECISION unavailable: no backend device");
      return 77;
   }
   GPU_CHECK(countStatus);

   std::puts("ENERGY-HIERARCHICAL-PRECISION: local/block/final accumulator "
             "precision swept through the literal RCG-09B block-tree "
             "reduction (reduceAdaptiveEnergyBlockT), against a Neumaier "
             "long-double host reference. Deterministic: one run per "
             "(case, N, mode).");

   const std::vector<std::uint64_t> sizes = {
      100ULL, 1000ULL, 10000ULL, 100000ULL, 1000000ULL, 3000000ULL
   };

   struct CaseGenerator {
      const char* name;
      std::vector<double> (*generate)(std::uint64_t);
   };
   const CaseGenerator generators[] = {
      {"all_positive", &caseAllPositive},
      {"alternating_sign", &caseAlternating},
      {"severe_cancellation", &caseSevereCancellation},
      {"wide_dynamic_range", &caseWideDynamicRange},
      {"randomized_physical", &caseRandomizedPhysical},
   };

   std::vector<CaseOutcome> allOutcomes;
   for(const auto& generator : generators) {
      for(const std::uint64_t n : sizes) {
         const std::vector<double> terms = generator.generate(n);
         CaseOutcome outcome;
         if(evaluateAllModes(generator.name, n, terms, outcome) != 0) return 1;
         printOutcome(outcome);
         allOutcomes.push_back(outcome);
      }
   }

   // Small-difference-of-large-sums case: N fixed at the production-scale
   // end of the sweep, delta spanning several orders down to a value
   // comparable to a tight adaptive-CG energy-jump-limit setting.
   const std::vector<double> deltas = {1.0e-1, 1.0e-3, 1.0e-6, 1.0e-9};
   for(const double delta : deltas) {
      if(evaluateDifferenceCase(1000000ULL, delta, allOutcomes) != 0) return 1;
      printOutcome(allOutcomes.back());
   }

   bool ok = true;

   // Negative control 1 (required by CGP-00): on well-conditioned sums (no
   // designed cancellation, reference not near zero), the FP64 reference
   // mode (F64/F64/F64) must stay near machine epsilon *relative to the
   // reference value* at every case and every N. This check deliberately
   // excludes severe_cancellation and energy_difference: those are designed
   // so the reference is a small residual of much larger terms, where even
   // exact arithmetic has relative error near (largest term magnitude /
   // reference magnitude) * epsilon -- large relative-to-residual error
   // there is the mathematics of cancellation, not a reduction defect. Their
   // own discriminating checks below use absolute-error ratios instead.
   double worstReferenceAbsError = 0.0;
   double worstReferenceRelError = 0.0;
   for(const auto& outcome : allOutcomes) {
      if(outcome.name == "severe_cancellation" ||
         outcome.name.rfind("energy_difference", 0) == 0)
         continue;
      for(const auto& mode : outcome.modes)
         if(std::strcmp(mode.label, "F64/F64/F64") == 0) {
            worstReferenceAbsError = std::max(worstReferenceAbsError, mode.absError);
            worstReferenceRelError = std::max(worstReferenceRelError, mode.relError);
         }
   }
   constexpr double referenceRelBudget = 1.0e-9;
   if(!(worstReferenceRelError < referenceRelBudget)) {
      std::printf("FAIL: F64/F64/F64 reference mode relative error %.6e "
                  "exceeded budget %.6e on a well-conditioned case -- "
                  "reference reduction regressed\n",
                  worstReferenceRelError, referenceRelBudget);
      ok = false;
   }

   // Negative control 2 (required by CGP-00): the cancellation-sensitive
   // test must discriminate a genuinely worse strategy (F32/F32/F32) from
   // the accepted reference (F64/F64/F64), by absolute error, at every
   // tested N -- not just the largest.
   bool cancellationDiscriminates = true;
   for(const auto& outcome : allOutcomes) {
      if(outcome.name != "severe_cancellation") continue;
      double f32Error = -1.0, f64Error = -1.0;
      for(const auto& mode : outcome.modes) {
         if(std::strcmp(mode.label, "F32/F32/F32") == 0) f32Error = mode.absError;
         if(std::strcmp(mode.label, "F64/F64/F64") == 0) f64Error = mode.absError;
      }
      if(!(f32Error > 2.0 * f64Error)) {
         std::printf("FAIL: severe_cancellation N=%llu did not discriminate "
                     "F32/F32/F32 (err=%.6e) from F64/F64/F64 (err=%.6e)\n",
                     static_cast<unsigned long long>(outcome.n), f32Error, f64Error);
         cancellationDiscriminates = false;
      }
   }
   if(!cancellationDiscriminates) ok = false;

   // Negative control 3 (required by CGP-00 Part D relevance): even at the
   // tightest, most cancellation-dominated energy difference tested (delta
   // comparable to a tight energy_jump_limit_j), F64/F64/F64 must still
   // discriminate from F32/F32/F32 by absolute error in the computed delta.
   bool deltaDiscriminates = true;
   for(const auto& outcome : allOutcomes) {
      if(outcome.name.rfind("energy_difference", 0) != 0) continue;
      double f32Error = -1.0, f64Error = -1.0;
      for(const auto& mode : outcome.modes) {
         if(std::strcmp(mode.label, "F32/F32/F32") == 0) f32Error = mode.absError;
         if(std::strcmp(mode.label, "F64/F64/F64") == 0) f64Error = mode.absError;
      }
      if(!(f32Error >= f64Error)) {
         std::printf("FAIL: %s did not discriminate F32/F32/F32 (err=%.6e) "
                     "from F64/F64/F64 (err=%.6e)\n",
                     outcome.name.c_str(), f32Error, f64Error);
         deltaDiscriminates = false;
      }
   }
   if(!deltaDiscriminates) ok = false;

   if(!ok) return 1;
   std::printf("ENERGY-HIERARCHICAL-PRECISION: PASS -- F64/F64/F64 stays "
               "within %.1e relative error on well-conditioned cases (worst "
               "abs=%.3e, worst rel=%.3e); severe_cancellation and "
               "energy_difference discriminate F32/F32/F32 from the "
               "reference by absolute error at every tested point\n",
               referenceRelBudget, worstReferenceAbsError, worstReferenceRelError);
   return 0;
}
