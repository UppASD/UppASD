// RCG-06B (F-11) ENERGY-FP32-ACCUM fixture.
//
// Isolates the exact defect class F-11 described -- a global energy
// accumulator reduced across many concurrent atomicAdd calls with no FP64
// storage anywhere in the device pipeline -- from the full adaptive-CG
// production physics, so its scaling behavior can be measured directly
// rather than assumed. Standalone and deliberately NOT linked against
// asdlib, matching RCG-06A's `test_stack_overflow_fault_injection.f90`
// precedent: this fixture must remain valid evidence for the accumulator
// defect class independent of later production source changes.
//
// It exercises the literal `atomicAddEnergyTerm` CAS-loop helper landed by
// this slice in gpuAtomicDouble.hpp (shared with gpuAdaptiveRuntime.cpp's
// production kernels, not a reimplementation that could silently drift)
// against a plain `atomicAdd(float*, float)` standing in for the pre-fix
// path -- native atomicAdd is only defined for `float` in the CUDA/HIP
// standard headers without opting into the (not guaranteed available, see
// gpuAtomicDouble.hpp) double intrinsic, which is exactly why F-11 was a
// real defect specifically in SINGLE_PREC (FP32 `real`) builds.
//
// Method: N threads each atomicAdd the same value = 1/N into a shared
// accumulator, once as `float` and once as `double` (computed independently
// in each type's own precision, not derived from one another, so the
// comparison isolates accumulator precision rather than input truncation).
// The exact mathematical sum is 1.0. atomicAdd interleaving order is
// hardware-scheduled and therefore not reproducible run to run; each size is
// repeated several times to show the resulting error is a stable band, not
// one lucky/unlucky draw.

#include "gpuAtomicDouble.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

#if defined(CUDA_V)
#include <cuda_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != cudaSuccess) { \
   std::fprintf(stderr, "%s\n", cudaGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE cudaErrorNoDevice
#define GPU_GET_COUNT cudaGetDeviceCount
#define GPU_MALLOC cudaMalloc
#define GPU_MEMCPY cudaMemcpy
#define GPU_MEMSET cudaMemset
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
#define GPU_MEMSET hipMemset
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

__global__ void accumulateFloatKernel(float* accumulator, float value,
                                      std::uint64_t n) {
   const std::uint64_t index =
      static_cast<std::uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   const std::uint64_t stride =
      static_cast<std::uint64_t>(gridDim.x) * blockDim.x;
   for(std::uint64_t i = index; i < n; i += stride)
      atomicAdd(accumulator, value);
}

__global__ void accumulateDoubleKernel(double* accumulator, double value,
                                       std::uint64_t n) {
   const std::uint64_t index =
      static_cast<std::uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   const std::uint64_t stride =
      static_cast<std::uint64_t>(gridDim.x) * blockDim.x;
   for(std::uint64_t i = index; i < n; i += stride)
      atomicAddEnergyTerm(accumulator, value);
}

struct Trial {
   double floatAbsError = 0.0;
   double doubleAbsError = 0.0;
};

int runTrial(std::uint64_t n, Trial& trial) {
   float* deviceFloatAcc = nullptr;
   double* deviceDoubleAcc = nullptr;
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceFloatAcc), sizeof(float)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceDoubleAcc), sizeof(double)));
   GPU_CHECK(GPU_MEMSET(deviceFloatAcc, 0, sizeof(float)));
   GPU_CHECK(GPU_MEMSET(deviceDoubleAcc, 0, sizeof(double)));

   const float valueF = 1.0f / static_cast<float>(n);
   const double valueD = 1.0 / static_cast<double>(n);

   const std::uint64_t maxGrid = 195312; // caps launch size well under limits
   const unsigned int gridBlocks = static_cast<unsigned int>(
      std::min<std::uint64_t>(maxGrid,
         (n + threadsPerBlock - 1) / threadsPerBlock));

   accumulateFloatKernel<<<gridBlocks, threadsPerBlock>>>(deviceFloatAcc, valueF, n);
   GPU_CHECK(GPU_GET_LAST_ERROR());
   accumulateDoubleKernel<<<gridBlocks, threadsPerBlock>>>(deviceDoubleAcc, valueD, n);
   GPU_CHECK(GPU_GET_LAST_ERROR());
   GPU_CHECK(GPU_DEVICE_SYNC());

   float floatResult = 0.0f;
   double doubleResult = 0.0;
   GPU_CHECK(GPU_MEMCPY(&floatResult, deviceFloatAcc, sizeof(float), GPU_MEMCPY_D2H));
   GPU_CHECK(GPU_MEMCPY(&doubleResult, deviceDoubleAcc, sizeof(double), GPU_MEMCPY_D2H));
   GPU_CHECK(GPU_FREE(deviceFloatAcc));
   GPU_CHECK(GPU_FREE(deviceDoubleAcc));

   trial.floatAbsError = std::fabs(static_cast<double>(floatResult) - 1.0);
   trial.doubleAbsError = std::fabs(doubleResult - 1.0);
   return 0;
}

} // namespace

int main() {
   int count = 0;
   const auto countStatus = GPU_GET_COUNT(&count);
   if(countStatus == GPU_NO_DEVICE || count == 0) {
      std::puts("ENERGY-FP32-ACCUM unavailable: no backend device");
      return 77;
   }
   GPU_CHECK(countStatus);

   // Every thread's atomicAdd to the single shared accumulator serializes
   // (that contention is the entire point -- it is what makes the
   // accumulation order, and therefore the rounding error, nondeterministic
   // and worth measuring rather than assuming). Serialized cost scales
   // ~linearly with N, so the top of this sweep is kept well short of the
   // full production active-atom/block count: it is large enough to show
   // the FP32 accumulator's error growing over 4 orders of magnitude while
   // keeping this a CTest-regression-suite-appropriate runtime (measured
   // ~30s uncontended; this repository's shared CI GPU can be several times
   // that under load from unrelated processes -- see the evidence doc).
   const std::vector<std::uint64_t> sizes = {
      1000ULL, 10000ULL, 100000ULL, 1000000ULL, 3000000ULL
   };
   constexpr int repeats = 3;

   std::printf("ENERGY-FP32-ACCUM: N threads atomicAdd 1/N into a shared "
               "accumulator (exact sum = 1.0); %d repeats per N\n", repeats);
   std::printf("%12s %14s %14s %14s %14s\n", "N", "float_mean", "float_max",
               "double_mean", "double_max");

   double smallestNFloatMeanError = 0.0;
   double largestNFloatMeanError = 0.0;
   double largestNFloatMaxError = 0.0;
   double worstDoubleErrorOverall = 0.0;

   for(std::size_t sizeIndex = 0; sizeIndex < sizes.size(); ++sizeIndex) {
      const std::uint64_t n = sizes[sizeIndex];
      double floatSum = 0.0, floatMax = 0.0;
      double doubleSum = 0.0, doubleMax = 0.0;
      for(int repeat = 0; repeat < repeats; ++repeat) {
         Trial trial;
         if(runTrial(n, trial) != 0) return 1;
         floatSum += trial.floatAbsError;
         floatMax = std::max(floatMax, trial.floatAbsError);
         doubleSum += trial.doubleAbsError;
         doubleMax = std::max(doubleMax, trial.doubleAbsError);
      }
      const double floatMean = floatSum / repeats;
      const double doubleMean = doubleSum / repeats;
      std::printf("%12llu %14.6e %14.6e %14.6e %14.6e\n",
                  static_cast<unsigned long long>(n), floatMean, floatMax,
                  doubleMean, doubleMax);
      std::fflush(stdout);

      if(sizeIndex == 0) smallestNFloatMeanError = floatMean;
      if(sizeIndex + 1 == sizes.size()) {
         largestNFloatMeanError = floatMean;
         largestNFloatMaxError = floatMax;
      }
      worstDoubleErrorOverall = std::max(worstDoubleErrorOverall, doubleMax);
   }

   bool ok = true;

   // Discriminating-test check: the FP32 accumulator's error must actually
   // grow with N (measured, not assumed -- atomicAdd interleaving order is
   // hardware-scheduled, so a flat/non-growing curve here would mean this
   // fixture is not exercising the defect at all).
   if(!(largestNFloatMeanError > 20.0 * smallestNFloatMeanError)) {
      std::printf("FAIL: float accumulator error did not grow with N "
                  "(smallest-N mean=%.6e, largest-N mean=%.6e) -- fixture is "
                  "not discriminating\n",
                  smallestNFloatMeanError, largestNFloatMeanError);
      ok = false;
   }

   // Fix negative control: the FP64 accumulator must stay near machine
   // epsilon at every size tested, including the largest.
   constexpr double doubleErrorBudget = 1.0e-9;
   if(!(worstDoubleErrorOverall < doubleErrorBudget)) {
      std::printf("FAIL: double accumulator error %.6e exceeded budget "
                  "%.6e at some tested N -- FP64 accumulation regressed\n",
                  worstDoubleErrorOverall, doubleErrorBudget);
      ok = false;
   }

   // The fix must be a genuine, large improvement at the largest N, not a
   // marginal difference.
   if(!(largestNFloatMaxError > 1000.0 * worstDoubleErrorOverall)) {
      std::printf("FAIL: at the largest N, float error %.6e was not >=1000x "
                  "the worst double error %.6e -- fix does not demonstrate a "
                  "material improvement\n",
                  largestNFloatMaxError, worstDoubleErrorOverall);
      ok = false;
   }

   if(!ok) return 1;
   std::puts("ENERGY-FP32-ACCUM: PASS -- FP32 accumulator error grows with "
             "N; FP64 accumulator (this slice's fix) stays near machine "
             "epsilon at every tested size");
   return 0;
}
