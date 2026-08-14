// RCG-09 / finding F-10.
//
// THIS IS NOT A PRODUCTION PERFORMANCE BASELINE AND MUST NEVER BE CITED AS ONE.
//
// This program answers exactly one question: does merely having an adaptive
// coarse-graining runtime object in scope, feature-disabled and never
// initialized, change the cost or the device inventory of unrelated GPU work?
// The "work" it times is a synthetic fused-multiply-add loop with no physical
// meaning -- not UppASD's Hamiltonian, not UppASD's integrator, not UppASD's
// memory access pattern.
//
// It lived inside the adaptive benchmark until RCG-09, where its output was
// read as though it were the atomistic path adaptive coarse graining has to
// beat.  It is not.  A real production baseline is what
// `gpu_adaptive_runtime_benchmark` measures (PERF-ATOMISTIC-PROD): UppASD's
// own GpuHamiltonianCalculations and GpuDepondtIntegrator, on the same
// geometry, at the same precision.  This file is retained only as the
// inactive-runtime overhead control required by the blueprint's
// "preserve the feature-off path" principle, and is separated so the two can
// never again be confused in one artifact.

#include "gpuAdaptiveRuntime.hpp"
#include "measurement/memoryMeasurement.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
   std::size_t elements = 8192;
   unsigned int iterations = 100;
   unsigned int rounds = 64;
   unsigned int repetitions = 7;
   double noiseFraction = 0.03;
   bool requireAcceptance = false;
};

std::size_t parseSize(const char* value, const char* option) {
   char* end = nullptr;
   const unsigned long long parsed = std::strtoull(value, &end, 10);
   if(!value[0] || !end || *end || parsed == 0)
      throw std::runtime_error(std::string("invalid ") + option);
   return static_cast<std::size_t>(parsed);
}

double parsePercent(const char* value, const char* option) {
   char* end = nullptr;
   const double parsed = std::strtod(value, &end);
   if(!value[0] || !end || *end || !std::isfinite(parsed) ||
      parsed <= 0.0 || parsed >= 100.0)
      throw std::runtime_error(std::string("invalid ") + option);
   return parsed / 100.0;
}

Options parse(int argc, char** argv) {
   Options result;
   for(int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if(option == "--elements" && index + 1 < argc)
         result.elements = parseSize(argv[++index], "--elements");
      else if(option == "--iterations" && index + 1 < argc)
         result.iterations = static_cast<unsigned int>(parseSize(argv[++index], "--iterations"));
      else if(option == "--rounds" && index + 1 < argc)
         result.rounds = static_cast<unsigned int>(parseSize(argv[++index], "--rounds"));
      else if(option == "--repetitions" && index + 1 < argc)
         result.repetitions = static_cast<unsigned int>(parseSize(argv[++index], "--repetitions"));
      else if(option == "--noise-percent" && index + 1 < argc)
         result.noiseFraction = parsePercent(argv[++index], "--noise-percent");
      else if(option == "--require-acceptance")
         result.requireAcceptance = true;
      else if(option == "--help") {
         std::puts(
            "usage: gpu_adaptive_inactive_overhead_microbenchmark "
            "[--elements N] [--iterations N] [--rounds N] [--repetitions N] "
            "[--noise-percent P] [--require-acceptance]\n"
            "\n"
            "Synthetic inactive-runtime overhead control.  NOT a production\n"
            "performance baseline; see gpu_adaptive_runtime_benchmark for that.");
         std::exit(0);
      } else {
         throw std::runtime_error("unknown or incomplete option: " + option);
      }
   }
   if(result.repetitions < 3)
      throw std::runtime_error("control requires at least 3 repetitions");
   return result;
}

void gpuCheck(GPU_ERROR_T status, const char* operation) {
   if(status != GPU_SUCCESS)
      throw std::runtime_error(std::string(operation) + ": " +
                               GPU_GET_ERROR_STRING(status));
}

double median(std::vector<double> values) {
   if(values.empty()) return 0.0;
   std::sort(values.begin(), values.end());
   const std::size_t middle = values.size() / 2;
   return values.size() % 2 ? values[middle] :
      0.5 * (values[middle - 1] + values[middle]);
}

double medianAbsoluteDeviation(const std::vector<double>& values) {
   const double center = median(values);
   std::vector<double> deviations(values.size());
   std::transform(values.begin(), values.end(), deviations.begin(),
                  [center](double value) { return std::abs(value - center); });
   return median(std::move(deviations));
}

// Synthetic arithmetic only.  Deliberately not named after any physics: it is
// a load, `rounds` fused multiply-adds, and a store.
__global__ void syntheticFmaWork(real* values, std::size_t count,
                                 unsigned int rounds) {
   const std::size_t index =
      static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
   if(index >= count) return;
   real value = values[index];
   for(unsigned int round = 0; round < rounds; ++round)
      value = value * static_cast<real>(1.0000001) + static_cast<real>(1.0e-7);
   values[index] = value;
}

void launchSyntheticWork(real* values, std::size_t count, unsigned int rounds,
                         GPU_STREAM_T stream,
                         const GpuAdaptiveRuntime* inactiveRuntime) {
   // The inactive runtime is dereferenced only for its state, never used to
   // dispatch: this is what "feature-disabled costs nothing" has to mean.
   if(inactiveRuntime && inactiveRuntime->ready())
      throw std::runtime_error("inactive runtime must not be ready");
   constexpr unsigned int threads = 256;
   const unsigned int blocks =
      static_cast<unsigned int>((count + threads - 1) / threads);
#if defined(CUDA_V)
   syntheticFmaWork<<<blocks, threads, 0, stream>>>(values, count, rounds);
#else
   hipLaunchKernelGGL(syntheticFmaWork, dim3(blocks), dim3(threads), 0,
                      stream, values, count, rounds);
#endif
}

double timeBatch(real* values, std::size_t count, const Options& options,
                 GPU_STREAM_T stream,
                 const GpuAdaptiveRuntime* inactiveRuntime) {
   GPU_EVENT_T begin{}, end{};
   gpuCheck(GPU_EVENT_CREATE(&begin), "begin-event creation");
   try {
      gpuCheck(GPU_EVENT_CREATE(&end), "end-event creation");
      gpuCheck(GPU_EVENT_RECORD(begin, stream), "begin-event record");
      for(unsigned int iteration = 0; iteration < options.iterations; ++iteration)
         launchSyntheticWork(values, count, options.rounds, stream, inactiveRuntime);
      gpuCheck(GPU_EVENT_RECORD(end, stream), "end-event record");
      gpuCheck(GPU_EVENT_SYNCHRONIZE(end), "event synchronization");
      float elapsed = 0.0f;
      gpuCheck(GPU_EVENT_ELAPSED_TIME(&elapsed, begin, end), "elapsed time");
      GPU_EVENT_DESTROY(begin);
      GPU_EVENT_DESTROY(end);
      return 1000.0 * static_cast<double>(elapsed) / options.iterations;
   } catch(...) {
      if(begin) GPU_EVENT_DESTROY(begin);
      if(end) GPU_EVENT_DESTROY(end);
      throw;
   }
}

} // namespace

int main(int argc, char** argv) {
   try {
      const Options options = parse(argc, argv);
      std::printf(
         "inactive-runtime-overhead-microbenchmark precision=%s backend=%s "
         "elements=%zu iterations=%u rounds=%u repetitions=%u\n",
         sizeof(real) == sizeof(double) ? "fp64" : "fp32",
#if defined(CUDA_V)
         "CUDA",
#else
         "HIP",
#endif
         options.elements, options.iterations, options.rounds,
         options.repetitions);
      std::puts(
         "inactive-runtime-overhead-microbenchmark note=SYNTHETIC_FMA_LOOP "
         "not_a_production_baseline=true "
         "production_baseline_is=gpu_adaptive_runtime_benchmark:PERF-ATOMISTIC-PROD");

      GpuTensor<real, 1> values;
      values.Allocate(static_cast<index_t>(options.elements));
      values.zeros();
      GPU_STREAM_T stream{};
      gpuCheck(GPU_STREAM_CREATE(&stream), "stream creation");
      std::vector<double> baseline, candidate;
      const auto inventory = TensorMemoryTracker::current_device();
      {
         GpuAdaptiveRuntime inactive;
         if(inactive.ready() || inactive.allocatedBytes() != 0 ||
            TensorMemoryTracker::current_device() != inventory)
            throw std::runtime_error("inactive runtime changed the device inventory");
         // Alternating order suppresses thermal and contention drift, which
         // matters on a shared device (RCG-08-FU3).
         for(unsigned int repetition = 0; repetition < options.repetitions; ++repetition) {
            if(repetition % 2 == 0) {
               baseline.push_back(timeBatch(values.data(), options.elements,
                                            options, stream, nullptr));
               candidate.push_back(timeBatch(values.data(), options.elements,
                                             options, stream, &inactive));
            } else {
               candidate.push_back(timeBatch(values.data(), options.elements,
                                             options, stream, &inactive));
               baseline.push_back(timeBatch(values.data(), options.elements,
                                            options, stream, nullptr));
            }
         }
      }
      gpuCheck(GPU_STREAM_DESTROY(stream), "stream destruction");
      values.Free();

      const double baselineMedian = median(baseline);
      const double candidateMedian = median(candidate);
      const double delta = (candidateMedian - baselineMedian) / baselineMedian;
      const double noise = 3.0 * (medianAbsoluteDeviation(baseline) +
                                  medianAbsoluteDeviation(candidate));
      const double budget = std::max(options.noiseFraction * baselineMedian, noise);
      const bool accepted = std::abs(candidateMedian - baselineMedian) <= budget;
      std::printf(
         "inactive-runtime-overhead baseline_us=%.6f candidate_us=%.6f "
         "delta_percent=%.6f baseline_mad_us=%.6f candidate_mad_us=%.6f "
         "budget_us=%.6f inventory_delta_bytes=0 result=%s\n",
         baselineMedian, candidateMedian, 100.0 * delta,
         medianAbsoluteDeviation(baseline), medianAbsoluteDeviation(candidate),
         budget, accepted ? "PASS" : "NOISY_OR_REGRESSED");
      for(std::size_t sample = 0; sample < baseline.size(); ++sample)
         std::printf("inactive-runtime-overhead-sample index=%zu baseline_us=%.6f "
                     "candidate_us=%.6f\n", sample, baseline[sample],
                     candidate[sample]);

      if(options.requireAcceptance && !accepted) return 2;
      return 0;
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU inactive-runtime overhead microbenchmark: %s\n",
                   error.what());
      return 1;
   }
}
