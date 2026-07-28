#include "gpuDipoleConvolution.hpp"
#include "dipoleOpenKernel.hpp"
#include "measurement/memoryMeasurement.h"
#include "stopwatch.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
   std::array<std::size_t, 3> grid{{32, 32, 32}};
   unsigned int basis = 1;
   unsigned int ensembles = 1;
   unsigned int warmup = 5;
   unsigned int iterations = 50;
   bool open = false;
   std::array<std::size_t, 3> fft_grid{};
};

std::size_t parseSize(const char* value, const char* option) {
   char* end = nullptr;
   const unsigned long long parsed = std::strtoull(value, &end, 10);
   if(!value[0] || !end || *end || parsed == 0) throw std::runtime_error(std::string("invalid ") + option);
   return static_cast<std::size_t>(parsed);
}

Options parse(int argc, char** argv) {
   Options options;
   for(int index = 1; index < argc; ++index) {
      const std::string option = argv[index];
      if(option == "--grid" && index + 3 < argc) {
         options.grid = {parseSize(argv[++index], "--grid"), parseSize(argv[++index], "--grid"),
                         parseSize(argv[++index], "--grid")};
      } else if(option == "--basis" && index + 1 < argc) {
         options.basis = static_cast<unsigned int>(parseSize(argv[++index], "--basis"));
      } else if(option == "--ensembles" && index + 1 < argc) {
         options.ensembles = static_cast<unsigned int>(parseSize(argv[++index], "--ensembles"));
      } else if(option == "--warmup" && index + 1 < argc) {
         options.warmup = static_cast<unsigned int>(parseSize(argv[++index], "--warmup"));
      } else if(option == "--iterations" && index + 1 < argc) {
         options.iterations = static_cast<unsigned int>(parseSize(argv[++index], "--iterations"));
      } else if(option == "--open") {
         options.open = true;
      } else if(option == "--fft-grid" && index + 3 < argc) {
         options.fft_grid = {parseSize(argv[++index], "--fft-grid"), parseSize(argv[++index], "--fft-grid"),
                             parseSize(argv[++index], "--fft-grid")};
      } else if(option == "--help") {
         std::puts("usage: dipole_gpu_fft_benchmark [--open] [--grid N1 N2 N3] [--fft-grid P1 P2 P3] [--basis N] [--ensembles N] [--warmup N] [--iterations N]");
         std::exit(0);
      } else {
         throw std::runtime_error("unknown or incomplete benchmark option: " + option);
      }
   }
   return options;
}

GpuDipoleConvolutionDescriptor descriptor(const Options& options) {
   static std::array<real, 3> c1, c2, c3;
   c1 = {static_cast<real>(1), static_cast<real>(0), static_cast<real>(0)};
   c2 = {static_cast<real>(0), static_cast<real>(1), static_cast<real>(0)};
   c3 = {static_cast<real>(0), static_cast<real>(0), static_cast<real>(1)};
   GpuDipoleConvolutionDescriptor result{};
   result.atomistic_grid = {options.grid[0], options.grid[1], options.grid[2]};
   result.macro_grid = result.atomistic_grid;
   result.basis = options.basis;
   result.ensembles = options.ensembles;
   result.boundary = options.open ? GpuDipoleBoundaryMode::Open : GpuDipoleBoundaryMode::Periodic3D;
   result.discretization = GpuDipoleDiscretization::MacrospinGrid;
   result.fft_grid = {options.fft_grid[0], options.fft_grid[1], options.fft_grid[2]};
   result.c1 = c1.data(); result.c2 = c2.data(); result.c3 = c3.data();
   result.basis_offsets.resize(options.basis);
   for(unsigned int basis = 0; basis < options.basis; ++basis) {
      // Keep every benchmark channel at a distinct physical point.  This is
      // immaterial for NA=1 and prevents a multi-basis sweep from accidentally
      // asking the finite builder to accept overlapping geometry.
      result.basis_offsets[basis] = {0.137 * basis, 0.071 * basis, 0.053 * basis};
   }
   result.alat = 1.0;
   result.tolerance = 1.0e-10;
   result.field_prefactor = 1.0;
   if(!result.valid()) throw std::runtime_error("benchmark descriptor is invalid");
   return result;
}

DipoleOpenGeometry openGeometry(const Options& options, const GpuDipoleConvolutionDescriptor& desc) {
   const auto active = desc.activeGrid();
   const auto fft = desc.fftGrid();
   DipoleOpenGeometry geometry{};
   geometry.atomistic_grid = {options.grid[0], options.grid[1], options.grid[2]};
   geometry.active_grid = {active.n1, active.n2, active.n3};
   geometry.fft_grid = {fft.n1, fft.n2, fft.n3};
   geometry.primitive_vectors = {{1.0, 0.0, 0.0,
                                   0.0, 1.0, 0.0,
                                   0.0, 0.0, 1.0}};
   geometry.basis = options.basis;
   geometry.basis_offsets = desc.basis_offsets;
   return geometry;
}

void printPhaseSamples(unsigned int iterations) {
   const auto samples = GlobalStopwatchPool::get("GPU dipole FFT").samples();
   for(const auto& sample : samples) {
      if(sample.name == "-") continue;
      std::printf("phase name=%s wall_us=%.6f gpu_us=%.6f\n", sample.name.c_str(),
                  1000.0 * static_cast<double>(sample.wall_ms) / iterations,
                  1000.0 * static_cast<double>(sample.gpu_ms) / iterations);
   }
}

} // namespace

int main(int argc, char** argv) {
   try {
      const Options options = parse(argc, argv);
      const auto desc = descriptor(options);
      const auto layout = desc.fftLayout();
      const std::size_t moments_count = 3 * layout.active_macros * options.ensembles;
      std::vector<real> moments(moments_count);
      for(std::size_t index = 0; index < moments.size(); ++index)
         moments[index] = static_cast<real>(0.125 + 0.0001 * static_cast<double>(index % 997));

      GPU_STREAM_T stream{};
      if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("benchmark stream creation failed");
      GpuTensor<real, 3> device_moments;
      GpuTensor<unsigned int, 1> device_map;
      GpuTensor<real, 3> device_beff, device_eneff;
      GpuTensor<real, 2> device_energy;
      GpuDipoleConvolution solver;
      GPU_EVENT_T begin{}, end{};
      try {
         device_moments.Allocate(3L, static_cast<long int>(layout.active_macros),
                                 static_cast<long int>(options.ensembles));
         device_map.Allocate(static_cast<long int>(layout.active_macros));
         device_beff.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(options.ensembles));
         device_eneff.Allocate(3L, static_cast<long int>(layout.active_macros), static_cast<long int>(options.ensembles));
         device_energy.Allocate(static_cast<long int>(options.ensembles), 7L);
         std::vector<unsigned int> map(layout.active_macros);
         std::iota(map.begin(), map.end(), 1U);
         if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real),
                       GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS ||
            GPU_MEMCPY(device_map.data(), map.data(), map.size() * sizeof(unsigned int),
                       GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS) {
            throw std::runtime_error("benchmark input upload failed");
         }
         device_beff.zeros_async(stream);
         device_eneff.zeros_async(stream);
         device_energy.zeros_async(stream);
         const auto builder_begin = std::chrono::steady_clock::now();
         DipoleOpenKernelResult open_kernel{};
         if(options.open) open_kernel = buildOpenDipoleDisplacementKernel(openGeometry(options, desc));
         const double host_builder_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - builder_begin).count();
         const auto setup_begin = std::chrono::steady_clock::now();
         if(!solver.initiate(desc, stream)) throw std::runtime_error("benchmark FFT plan initiation failed");
         if(options.open) {
            solver.uploadOpenKernel(open_kernel.kernel);
         } else {
            // The periodic benchmark deliberately stays a plumbing benchmark:
            // OPEN_FFT always uses its finite production builder above.
            std::vector<double> kernel(layout.fft_cells * layout.kernel_batches);
            for(std::size_t index = 0; index < kernel.size(); ++index)
               kernel[index] = 0.001 * static_cast<double>(1 + index % 37);
            solver.uploadCompleteKernelForTesting(kernel, false);
         }
         const double setup_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - setup_begin).count();
         for(unsigned int iteration = 0; iteration < options.warmup; ++iteration) {
            solver.evaluate(device_moments);
            solver.addFieldsToAtoms(device_beff, device_eneff, device_map.data(), layout.active_macros);
            solver.accumulateEnergy(device_energy, layout.active_macros);
         }
         if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("benchmark warmup synchronization failed");
         // The production code's event-backed phase timers are disabled during
         // warmup so each printed phase is an exact per-steady-state evaluation.
         Stopwatch::setTimingMode('Y');
         GlobalStopwatchPool::get("GPU dipole FFT").reset();
         if(GPU_EVENT_CREATE(&begin) != GPU_SUCCESS || GPU_EVENT_CREATE(&end) != GPU_SUCCESS)
            throw std::runtime_error("benchmark event creation failed");
         GPU_EVENT_RECORD(begin, stream);
         for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) {
            solver.evaluate(device_moments);
            solver.addFieldsToAtoms(device_beff, device_eneff, device_map.data(), layout.active_macros);
            solver.accumulateEnergy(device_energy, layout.active_macros);
         }
         GPU_EVENT_RECORD(end, stream);
         if(GPU_EVENT_SYNCHRONIZE(end) != GPU_SUCCESS) throw std::runtime_error("benchmark event synchronization failed");
         float elapsed_ms = 0.0f;
         if(GPU_EVENT_ELAPSED_TIME(&elapsed_ms, begin, end) != GPU_SUCCESS)
            throw std::runtime_error("benchmark event timing failed");
         const double per_evaluation_us = 1000.0 * static_cast<double>(elapsed_ms) / options.iterations;
         const auto fft = desc.fftGrid();
         std::printf("dipole-benchmark mode=%s precision=%s grid=%zux%zux%zu fft_grid=%zux%zux%zu basis=%u ensembles=%u warmup=%u iterations=%u "
                     "host_builder_ms=%.6f setup_ms=%.6f steady_ms=%.6f steady_us=%.6f persistent_bytes=%zu construction_bytes=%zu workspace_bytes=%zu tracker_current_bytes=%lld tracker_peak_bytes=%lld\n",
                     options.open ? "OPEN_FFT" : "EWALD3D_FFT-plumbing",
                     sizeof(real) == sizeof(double) ? "fp64" : "fp32", options.grid[0], options.grid[1], options.grid[2],
                     fft.n1, fft.n2, fft.n3, options.basis, options.ensembles, options.warmup, options.iterations,
                     host_builder_ms, setup_ms, static_cast<double>(elapsed_ms), per_evaluation_us,
                     GpuDipoleConvolution::estimatePersistentBytes(desc),
                     GpuDipoleConvolution::estimateConstructionBytes(desc),
                     GpuDipoleConvolution::estimateWorkspaceBytes(desc),
                     static_cast<long long>(TensorMemoryTracker::current_device()),
                     static_cast<long long>(TensorMemoryTracker::peak_device()));
         GPU_EVENT_DESTROY(begin); GPU_EVENT_DESTROY(end);
         begin = {}; end = {};
         solver.release();
         printPhaseSamples(options.iterations);
         Stopwatch::setTimingMode('N');
         device_energy.Free(); device_eneff.Free(); device_beff.Free(); device_map.Free(); device_moments.Free();
         if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("benchmark stream destruction failed");
      } catch(...) {
         if(begin) GPU_EVENT_DESTROY(begin);
         if(end) GPU_EVENT_DESTROY(end);
         solver.release();
         if(!device_energy.empty()) device_energy.Free();
         if(!device_eneff.empty()) device_eneff.Free();
         if(!device_beff.empty()) device_beff.Free();
         if(!device_map.empty()) device_map.Free();
         if(!device_moments.empty()) device_moments.Free();
         GPU_STREAM_DESTROY(stream);
         throw;
      }
   } catch(const std::exception& error) {
      std::fprintf(stderr, "FAIL GPU dipole FFT benchmark: %s\n", error.what());
      return 1;
   }
   return 0;
}
