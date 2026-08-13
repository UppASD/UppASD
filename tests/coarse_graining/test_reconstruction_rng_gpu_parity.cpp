// RCG-06D (F-20) RECONSTRUCTION-RNG-SPATIAL-STATS GPU companion.
//
// Standalone CUDA/HIP executable, deliberately NOT linked against asdlib
// (matching RCG-06A's stack_overflow_fault_injection/RCG-06B's
// energy_fp32_accum precedent), that computes the same
// (global_seed,channel,ensemble,epoch,block) sweep as
// test_reconstruction_rng_spatial_stats.f90's CPU driver, using the literal
// device generator gpuAdaptiveRuntime.cpp's production kernels use
// (adaptiveTupleSeed/adaptiveNextUniform, gpuAdaptiveReconstructionRng.hpp
// -- shared, not reimplemented).
//
// run_reconstruction_rng_spatial_stats.py compares this CSV output against
// the CPU driver's, keyed by (global_seed,channel,ensemble,epoch,block), as
// the regression test for this slice's two CPU/GPU fixes to
// gpuAdaptiveReconstructionRng.hpp (see that header for the full account of
// both): the generator multiplier (16807 -> 48271, matching CPU) and the
// seed formula's spurious "+1" on epoch (removed; block/channel/ensemble
// keep their "+1", which is a genuine 0-based/1-based index-convention
// conversion, not a bug). Reverting either fix in the shared header
// reproduces a large, easily distinguished mismatch in this comparison.

#include "gpuAdaptiveReconstructionRng.hpp"

#include <cstdint>
#include <cstdio>
#include <vector>

#if defined(CUDA_V)
#include <cuda_runtime.h>
#define GPU_CHECK(call) do { const auto error = (call); if (error != cudaSuccess) { \
   std::fprintf(stderr, "%s\n", cudaGetErrorString(error)); return 1; } } while (false)
#define GPU_NO_DEVICE cudaErrorNoDevice
#define GPU_GET_COUNT cudaGetDeviceCount
#define GPU_MALLOC cudaMalloc
#define GPU_MEMCPY_H2D(dst, src, bytes) cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice)
#define GPU_MEMCPY_D2H(dst, src, bytes) cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToHost)
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
#define GPU_MEMCPY_H2D(dst, src, bytes) hipMemcpy(dst, src, bytes, hipMemcpyHostToDevice)
#define GPU_MEMCPY_D2H(dst, src, bytes) hipMemcpy(dst, src, bytes, hipMemcpyDeviceToHost)
#define GPU_FREE hipFree
#define GPU_DEVICE_SYNC hipDeviceSynchronize
#define GPU_GET_LAST_ERROR hipGetLastError
#else
#error "This fixture requires the configured CUDA or HIP backend"
#endif

namespace {

// Must match test_reconstruction_rng_spatial_stats.f90's sweep exactly, so
// run_reconstruction_rng_spatial_stats.py compares like-for-like rows.
constexpr int nblocks = 4096;
constexpr int nseeds = 3;
constexpr int nchannels = 2;
constexpr int nensembles = 2;
constexpr int nepochs = 2;
constexpr std::uint64_t globalSeeds[nseeds] = {1ULL, 99173ULL, 20260812ULL};
constexpr int channelValues[nchannels] = {1, 2};
constexpr int ensembleValues[nensembles] = {1, 2};
constexpr unsigned int epochValues[nepochs] = {0u, 7u};
constexpr std::uint64_t totalRows =
   static_cast<std::uint64_t>(nseeds) * nchannels * nensembles * nepochs * nblocks;

__global__ void sweepKernel(std::uint64_t* seedOut, int* channelOut, int* ensembleOut,
                            unsigned int* epochOut, int* blockOut, real* u1Out, real* u2Out,
                            const std::uint64_t* seeds, const int* chans, const int* ens,
                            const unsigned int* eps) {
   const std::uint64_t index = static_cast<std::uint64_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   if(index >= totalRows) return;
   std::uint64_t remaining = index;
   const int block = static_cast<int>(remaining % nblocks) + 1;
   remaining /= nblocks;
   const int epochIndex = static_cast<int>(remaining % nepochs);
   remaining /= nepochs;
   const int ensembleIndex = static_cast<int>(remaining % nensembles);
   remaining /= nensembles;
   const int channelIndex = static_cast<int>(remaining % nchannels);
   remaining /= nchannels;
   const int seedIndex = static_cast<int>(remaining);

   const std::uint64_t globalSeed = seeds[seedIndex];
   const int channel = chans[channelIndex];
   const int ensemble = ens[ensembleIndex];
   const unsigned int epoch = eps[epochIndex];

   std::uint64_t rng = adaptiveTupleSeed(globalSeed, static_cast<std::size_t>(block - 1),
      static_cast<std::size_t>(channel - 1), static_cast<std::size_t>(ensemble - 1), epoch);
   const real u1 = adaptiveNextUniform(rng);
   const real u2 = adaptiveNextUniform(rng);

   seedOut[index] = globalSeed;
   channelOut[index] = channel;
   ensembleOut[index] = ensemble;
   epochOut[index] = epoch;
   blockOut[index] = block;
   u1Out[index] = u1;
   u2Out[index] = u2;
}

} // namespace

int main() {
   int count = 0;
   const auto countStatus = GPU_GET_COUNT(&count);
   if(countStatus == GPU_NO_DEVICE || count == 0) {
      std::puts("RECONSTRUCTION-RNG-GPU-PARITY unavailable: no backend device");
      return 77;
   }
   GPU_CHECK(countStatus);

   std::uint64_t* deviceSeeds = nullptr;
   int* deviceChannels = nullptr;
   int* deviceEnsembles = nullptr;
   unsigned int* deviceEpochs = nullptr;
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceSeeds), sizeof(globalSeeds)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceChannels), sizeof(channelValues)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceEnsembles), sizeof(ensembleValues)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceEpochs), sizeof(epochValues)));
   GPU_CHECK(GPU_MEMCPY_H2D(deviceSeeds, globalSeeds, sizeof(globalSeeds)));
   GPU_CHECK(GPU_MEMCPY_H2D(deviceChannels, channelValues, sizeof(channelValues)));
   GPU_CHECK(GPU_MEMCPY_H2D(deviceEnsembles, ensembleValues, sizeof(ensembleValues)));
   GPU_CHECK(GPU_MEMCPY_H2D(deviceEpochs, epochValues, sizeof(epochValues)));

   std::uint64_t* deviceSeedOut = nullptr;
   int* deviceChannelOut = nullptr;
   int* deviceEnsembleOut = nullptr;
   unsigned int* deviceEpochOut = nullptr;
   int* deviceBlockOut = nullptr;
   real* deviceU1Out = nullptr;
   real* deviceU2Out = nullptr;
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceSeedOut), totalRows * sizeof(std::uint64_t)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceChannelOut), totalRows * sizeof(int)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceEnsembleOut), totalRows * sizeof(int)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceEpochOut), totalRows * sizeof(unsigned int)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceBlockOut), totalRows * sizeof(int)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceU1Out), totalRows * sizeof(real)));
   GPU_CHECK(GPU_MALLOC(reinterpret_cast<void**>(&deviceU2Out), totalRows * sizeof(real)));

   constexpr unsigned int threadsPerBlock = 256;
   const unsigned int gridBlocks =
      static_cast<unsigned int>((totalRows + threadsPerBlock - 1) / threadsPerBlock);
   sweepKernel<<<gridBlocks, threadsPerBlock>>>(deviceSeedOut, deviceChannelOut, deviceEnsembleOut,
      deviceEpochOut, deviceBlockOut, deviceU1Out, deviceU2Out, deviceSeeds, deviceChannels,
      deviceEnsembles, deviceEpochs);
   GPU_CHECK(GPU_GET_LAST_ERROR());
   GPU_CHECK(GPU_DEVICE_SYNC());

   std::vector<std::uint64_t> seedOut(totalRows);
   std::vector<int> channelOut(totalRows), ensembleOut(totalRows), blockOut(totalRows);
   std::vector<unsigned int> epochOut(totalRows);
   std::vector<real> u1Out(totalRows), u2Out(totalRows);
   GPU_CHECK(GPU_MEMCPY_D2H(seedOut.data(), deviceSeedOut, totalRows * sizeof(std::uint64_t)));
   GPU_CHECK(GPU_MEMCPY_D2H(channelOut.data(), deviceChannelOut, totalRows * sizeof(int)));
   GPU_CHECK(GPU_MEMCPY_D2H(ensembleOut.data(), deviceEnsembleOut, totalRows * sizeof(int)));
   GPU_CHECK(GPU_MEMCPY_D2H(epochOut.data(), deviceEpochOut, totalRows * sizeof(unsigned int)));
   GPU_CHECK(GPU_MEMCPY_D2H(blockOut.data(), deviceBlockOut, totalRows * sizeof(int)));
   GPU_CHECK(GPU_MEMCPY_D2H(u1Out.data(), deviceU1Out, totalRows * sizeof(real)));
   GPU_CHECK(GPU_MEMCPY_D2H(u2Out.data(), deviceU2Out, totalRows * sizeof(real)));

   GPU_CHECK(GPU_FREE(deviceSeeds));
   GPU_CHECK(GPU_FREE(deviceChannels));
   GPU_CHECK(GPU_FREE(deviceEnsembles));
   GPU_CHECK(GPU_FREE(deviceEpochs));
   GPU_CHECK(GPU_FREE(deviceSeedOut));
   GPU_CHECK(GPU_FREE(deviceChannelOut));
   GPU_CHECK(GPU_FREE(deviceEnsembleOut));
   GPU_CHECK(GPU_FREE(deviceEpochOut));
   GPU_CHECK(GPU_FREE(deviceBlockOut));
   GPU_CHECK(GPU_FREE(deviceU1Out));
   GPU_CHECK(GPU_FREE(deviceU2Out));

   std::puts("# RCG-06D reconstruction RNG spatial statistics raw draws (GPU)");
   std::puts("# global_seed channel ensemble epoch block u1 u2");
   for(std::uint64_t i = 0; i < totalRows; ++i) {
      std::printf("%llu %d %d %u %d %.16e %.16e\n",
         static_cast<unsigned long long>(seedOut[i]), channelOut[i], ensembleOut[i], epochOut[i],
         blockOut[i], static_cast<double>(u1Out[i]), static_cast<double>(u2Out[i]));
   }
   return 0;
}
