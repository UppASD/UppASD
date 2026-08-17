// CGP-06A: isolate the production thermal-RNG cost from the rest of the
// Depondt predictor.
//
// GpuDepondtIntegrator::evolveFirst()/evolveFirst(activeList, count) both call
// GpuThermfield::randomize() unconditionally on the full (N,M) site/ensemble
// layout, before either the ordinary or active-atom-scoped predictor kernel
// runs (see source/gpu_files/gpuDepondtIntegrator.cpp, evolveFirst()'s own
// comment: "Keep the production generator and its full (site,ensemble)
// layout"). This harness exercises the real GpuThermfield class exactly as
// production does -- same initiate()/initiateConstants()/randomize() calls,
// same curand/hiprand generator type -- so the measured cost is the real
// production cost, not a synthetic proxy. It never modifies production code
// or behaviour; it is read-only instrumentation of an existing class.
//
// One (N, M, temperature) configuration per process invocation, deliberately.
// Project memory ("CGP-05 host-barrier removal") records that repeated
// construct/destroy of GpuThermfield-owning objects together with the
// process-global ParallelizationHelperInstance singleton, in one process, can
// destabilize curand's stream binding -- exactly the multi-instance churn
// pattern an earlier version of this file used, and which segfaulted inside
// cudaEventRecord on the very first construction. Production only ever
// constructs one GpuThermfield per process; this harness now matches that,
// and a driver script sweeps sizes across separate process invocations
// instead.
//
// Usage: gpu_thermfield_rng_benchmark --atoms N [--ensembles M]
//        [--temperature T] [--seed S] [--warmup W] [--iterations I]
//        [--repetitions R] [--subphase] [--active-fraction F]
//
// Prints one "thermfield-rng" result line with the randomize() wall-clock
// median/MAD. --subphase additionally reports the GPU-event breakdown of
// randomize() into its own "RNG" (curand/hiprand generate) and "loop"
// (SetupField write) categories, read from the same
// GlobalStopwatchPool("GPU thermfield") production already populates when
// phase timing is enabled -- not new instrumentation, just read out here.
//
// CGP-06B: --active-fraction F (0 <= F <= 1) switches the timed call from
// randomize(mmom) to the real, production active-scoped overload
// randomize(mmom, oneBasedAtoms, activeCount), against a compact list of
// round(F*atoms) one-based ids built once before timing starts. This
// exercises the actual Strategy 1 kernel geometry (gpuActiveAtomSiteCall),
// not the CGP-06A "re-run at a smaller full lattice" projection -- the RNG
// generate sub-phase still draws the full 3*N*M values every call (by
// design: the generator sequence must stay identical to the unscoped path),
// only the SetupField write is scoped, so the real active-fraction win is
// necessarily smaller than that superseded projection implied.

#include "gpuThermfield.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"
#include "real_type.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool backendAvailable() {
   int count = 0;
#if defined(CUDA_V)
   const auto status = cudaGetDeviceCount(&count);
   return status == cudaSuccess && count > 0 && cudaFree(nullptr) == cudaSuccess;
#elif defined(HIP_V)
   const auto status = hipGetDeviceCount(&count);
   return status == hipSuccess && count > 0 && hipFree(nullptr) == hipSuccess;
#endif
}

void gpuCheck(GPU_ERROR_T status, const char* operation) {
   if(status != GPU_SUCCESS)
      throw std::runtime_error(std::string(operation) + ": " + GPU_GET_ERROR_STRING(status));
}

double median(std::vector<double> values) {
   if(values.empty()) return 0.0;
   std::sort(values.begin(), values.end());
   const std::size_t middle = values.size() / 2;
   return values.size() % 2 ? values[middle] : 0.5 * (values[middle - 1] + values[middle]);
}

double medianAbsoluteDeviation(const std::vector<double>& values) {
   const double center = median(values);
   std::vector<double> deviations(values.size());
   std::transform(values.begin(), values.end(), deviations.begin(),
                  [center](double value) { return std::abs(value - center); });
   return median(std::move(deviations));
}

struct Options {
   std::size_t atoms = 0;
   std::size_t ensembles = 1;
   real temperature = real(0.0);
   unsigned long long seed = 918273ULL;
   unsigned int warmup = 5;
   unsigned int iterations = 30;
   unsigned int repetitions = 7;
   bool subphase = false;
   bool activeScoped = false;
   double activeFraction = 1.0;
};

std::size_t parseSize(const char* value, const char* option) {
   char* end = nullptr;
   const unsigned long long parsed = std::strtoull(value, &end, 10);
   if(!value[0] || !end || *end || parsed == 0)
      throw std::runtime_error(std::string("invalid ") + option);
   return static_cast<std::size_t>(parsed);
}

double parseDouble(const char* value, const char* option) {
   char* end = nullptr;
   const double parsed = std::strtod(value, &end);
   if(!value[0] || !end || *end)
      throw std::runtime_error(std::string("invalid ") + option);
   return parsed;
}

Options parse(int argc, char** argv) {
   Options result;
   for(int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if(option == "--atoms" && index + 1 < argc)
         result.atoms = parseSize(argv[++index], "--atoms");
      else if(option == "--ensembles" && index + 1 < argc)
         result.ensembles = parseSize(argv[++index], "--ensembles");
      else if(option == "--temperature" && index + 1 < argc)
         result.temperature = static_cast<real>(parseDouble(argv[++index], "--temperature"));
      else if(option == "--seed" && index + 1 < argc)
         result.seed = std::strtoull(argv[++index], nullptr, 10);
      else if(option == "--warmup" && index + 1 < argc)
         result.warmup = static_cast<unsigned int>(parseSize(argv[++index], "--warmup"));
      else if(option == "--iterations" && index + 1 < argc)
         result.iterations = static_cast<unsigned int>(parseSize(argv[++index], "--iterations"));
      else if(option == "--repetitions" && index + 1 < argc)
         result.repetitions = static_cast<unsigned int>(parseSize(argv[++index], "--repetitions"));
      else if(option == "--subphase")
         result.subphase = true;
      else if(option == "--active-fraction" && index + 1 < argc) {
         result.activeScoped = true;
         result.activeFraction = parseDouble(argv[++index], "--active-fraction");
      } else
         throw std::runtime_error("unknown or incomplete option: " + option);
   }
   if(result.atoms == 0) throw std::runtime_error("--atoms is required and must be > 0");
   if(result.activeFraction < 0.0 || result.activeFraction > 1.0)
      throw std::runtime_error("--active-fraction must be in [0,1]");
   return result;
}

struct TimingSample {
   std::vector<double> wallUs;
   double medianUs = 0.0;
   double madUs = 0.0;
   void finish() {
      medianUs = median(wallUs);
      madUs = medianAbsoluteDeviation(wallUs);
   }
};

template <typename Body>
TimingSample timeSteadyState(const Options& options, Body body) {
   for(unsigned int iteration = 0; iteration < options.warmup; ++iteration) body();
   gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "warm-up synchronization");
   TimingSample sample;
   for(unsigned int repetition = 0; repetition < options.repetitions; ++repetition) {
      const auto begin = std::chrono::steady_clock::now();
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) body();
      gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "steady-state synchronization");
      sample.wallUs.push_back(1.0e6 *
         std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count() /
         options.iterations);
   }
   sample.finish();
   return sample;
}

} // namespace

int main(int argc, char** argv) {
   if(!backendAvailable()) {
      std::puts("THERMFIELD-RNG-BENCHMARK unavailable: no backend device");
      return 77;
   }

   try {
      const Options options = parse(argc, argv);

      ParallelizationHelperInstance.initiate(static_cast<unsigned int>(options.atoms),
                                             static_cast<unsigned int>(options.ensembles), 0);
      GpuThermfield field;
      if(!field.initiate(options.atoms, options.ensembles, GPU_RAND_RNG_PSEUDO_DEFAULT,
                         options.seed))
         throw std::runtime_error("GpuThermfield::initiate failed");
      GpuTensor<real, 2> mmom;
      mmom.Allocate(options.atoms, options.ensembles);
      {
         std::vector<real> moments(options.atoms * options.ensembles, real(1));
         gpuCheck(GPU_MEMCPY(mmom.data(), moments.data(), moments.size() * sizeof(real),
                             GPU_MEMCPY_HOST_TO_DEVICE),
                  "mmom upload");
      }
      Tensor<real, 1> temperatureHost;
      temperatureHost.AllocateHost(options.atoms);
      for(std::size_t site = 0; site < options.atoms; ++site)
         temperatureHost(site) = options.temperature;
      if(!field.initiateConstants(temperatureHost, real(0.001), real(1.0), real(1.0), real(1.0),
                                  real(0.1)))
         throw std::runtime_error("GpuThermfield::initiateConstants failed");

      // Active list, only built/uploaded when requested: round(F*atoms)
      // one-based ids 1..activeCount. Cost depends only on count (Part C),
      // not on which physical sites are selected, so a contiguous prefix is
      // representative.
      GpuTensor<int, 1> activeList;
      unsigned int activeCount = 0;
      if(options.activeScoped) {
         activeCount = static_cast<unsigned int>(
            std::llround(options.activeFraction * static_cast<double>(options.atoms)));
         if(activeCount > 0) {
            std::vector<int> oneBasedHost(activeCount);
            for(unsigned int i = 0; i < activeCount; ++i) oneBasedHost[i] = static_cast<int>(i + 1);
            activeList.Allocate(static_cast<index_t>(activeCount));
            Tensor<int, 1> view(oneBasedHost.data(), static_cast<index_t>(activeCount));
            activeList.copy_sync(view);
         }
      }
      const int* activeData = activeCount > 0 ? activeList.data() : nullptr;

      const auto sample = options.activeScoped
         ? timeSteadyState(options, [&]() { field.randomize(mmom, activeData, activeCount); })
         : timeSteadyState(options, [&]() { field.randomize(mmom); });
      if(options.activeScoped) {
         std::printf(
            "thermfield-rng-active atoms=%zu ensembles=%zu temperature=%.6f "
            "requested_fraction=%.6f active_count=%u "
            "randomize_wall_us=%.6f randomize_mad_us=%.6f warmup=%u iterations=%u "
            "repetitions=%u\n",
            options.atoms, options.ensembles, static_cast<double>(options.temperature),
            options.activeFraction, activeCount,
            sample.medianUs, sample.madUs, options.warmup, options.iterations,
            options.repetitions);
      } else {
         std::printf(
            "thermfield-rng atoms=%zu ensembles=%zu temperature=%.6f "
            "randomize_wall_us=%.6f randomize_mad_us=%.6f warmup=%u iterations=%u "
            "repetitions=%u\n",
            options.atoms, options.ensembles, static_cast<double>(options.temperature),
            sample.medianUs, sample.madUs, options.warmup, options.iterations,
            options.repetitions);
      }

      if(options.subphase) {
         const auto call = [&]() {
            if(options.activeScoped) field.randomize(mmom, activeData, activeCount);
            else field.randomize(mmom);
         };
         Stopwatch::setTimingMode('Y');
         call();
         gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "sub-phase warm-up sync");
         Stopwatch& watch = GlobalStopwatchPool::get("GPU thermfield");
         watch.reset();
         const unsigned int phaseIterations = options.iterations;
         for(unsigned int i = 0; i < phaseIterations; ++i) call();
         gpuCheck(GPU_DEVICE_SYNCHRONIZE(), "sub-phase sync");
         for(const auto& phaseSample : watch.samples()) {
            std::printf(
               "thermfield-subphase atoms=%zu ensembles=%zu name=%s gpu_us_per_call=%.6f "
               "wall_us_per_call=%.6f\n",
               options.atoms, options.ensembles, phaseSample.name.c_str(),
               1000.0 * phaseSample.gpu_ms / phaseIterations,
               1000.0 * phaseSample.wall_ms / phaseIterations);
         }
         Stopwatch::setTimingMode('N');
      }

      activeList.Free();
      mmom.Free();
      temperatureHost.FreeHost();
      ParallelizationHelperInstance.free();
      std::puts("THERMFIELD-RNG-BENCHMARK passed");
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "%s\n", error.what());
      ParallelizationHelperInstance.free();
      return 1;
   }
}
