#include "gpuDipoleConvolution.hpp"
#include "tensor.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
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
      } else if(option == "--help") {
         std::puts("usage: dipole_gpu_fft_benchmark [--grid N1 N2 N3] [--basis N] [--ensembles N] [--warmup N] [--iterations N]");
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
   result.boundary = GpuDipoleBoundaryMode::Periodic3D;
   result.discretization = GpuDipoleDiscretization::MacrospinGrid;
   result.c1 = c1.data(); result.c2 = c2.data(); result.c3 = c3.data();
   result.basis_offsets.assign(options.basis, {0.0, 0.0, 0.0});
   result.alat = 1.0;
   result.tolerance = 1.0e-10;
   result.field_prefactor = 1.0;
   if(!result.valid()) throw std::runtime_error("benchmark descriptor is invalid");
   return result;
}

} // namespace

int main(int argc, char** argv) {
   try {
      const Options options = parse(argc, argv);
      const auto desc = descriptor(options);
      const auto layout = desc.fftLayout();
      const std::size_t moments_count = 3 * layout.real_cells * options.basis * options.ensembles;
      const std::size_t kernel_count = layout.real_cells * layout.kernel_batches;
      std::vector<real> moments(moments_count);
      std::vector<double> kernel(kernel_count);
      for(std::size_t index = 0; index < moments.size(); ++index)
         moments[index] = static_cast<real>(0.125 + 0.0001 * static_cast<double>(index % 997));
      // A non-physical tensor is intentional: this benchmark measures the
      // runtime FFT/operator pipeline without charging a representative sweep
      // for the one-time Ewald construction.  Physics validation remains in
      // dipole_gpu_fft_tests.
      for(std::size_t index = 0; index < kernel.size(); ++index)
         kernel[index] = 0.001 * static_cast<double>(1 + index % 37);

      GPU_STREAM_T stream{};
      if(GPU_STREAM_CREATE(&stream) != GPU_SUCCESS) throw std::runtime_error("benchmark stream creation failed");
      GpuTensor<real, 3> device_moments;
      GpuDipoleConvolution solver;
      GPU_EVENT_T begin{}, end{};
      try {
         device_moments.Allocate(3L, static_cast<long int>(layout.real_cells * options.basis),
                                 static_cast<long int>(options.ensembles));
         if(GPU_MEMCPY(device_moments.data(), moments.data(), moments.size() * sizeof(real),
                       GPU_MEMCPY_HOST_TO_DEVICE) != GPU_SUCCESS) {
            throw std::runtime_error("benchmark moment upload failed");
         }
         const auto setup_begin = std::chrono::steady_clock::now();
         if(!solver.initiate(desc, stream)) throw std::runtime_error("benchmark FFT plan initiation failed");
         solver.uploadCompleteKernelForTesting(kernel, false);
         const double setup_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - setup_begin).count();
         for(unsigned int iteration = 0; iteration < options.warmup; ++iteration) solver.evaluate(device_moments);
         if(GPU_STREAM_SYNC(stream) != GPU_SUCCESS) throw std::runtime_error("benchmark warmup synchronization failed");
         if(GPU_EVENT_CREATE(&begin) != GPU_SUCCESS || GPU_EVENT_CREATE(&end) != GPU_SUCCESS)
            throw std::runtime_error("benchmark event creation failed");
         GPU_EVENT_RECORD(begin, stream);
         for(unsigned int iteration = 0; iteration < options.iterations; ++iteration) solver.evaluate(device_moments);
         GPU_EVENT_RECORD(end, stream);
         if(GPU_EVENT_SYNCHRONIZE(end) != GPU_SUCCESS) throw std::runtime_error("benchmark event synchronization failed");
         float elapsed_ms = 0.0f;
         if(GPU_EVENT_ELAPSED_TIME(&elapsed_ms, begin, end) != GPU_SUCCESS)
            throw std::runtime_error("benchmark event timing failed");
         const double per_evaluation_us = 1000.0 * static_cast<double>(elapsed_ms) / options.iterations;
         std::printf("dipole-benchmark precision=%s grid=%zux%zux%zu basis=%u ensembles=%u warmup=%u iterations=%u "
                     "setup_ms=%.6f eval_ms=%.6f eval_us=%.6f persistent_bytes=%zu peak_bytes=%zu\n",
                     sizeof(real) == sizeof(double) ? "fp64" : "fp32", options.grid[0], options.grid[1], options.grid[2],
                     options.basis, options.ensembles, options.warmup, options.iterations, setup_ms,
                     static_cast<double>(elapsed_ms), per_evaluation_us,
                     GpuDipoleConvolution::estimatePersistentBytes(desc), GpuDipoleConvolution::estimateBytes(desc));
         GPU_EVENT_DESTROY(begin); GPU_EVENT_DESTROY(end);
         solver.release(); device_moments.Free();
         if(GPU_STREAM_DESTROY(stream) != GPU_SUCCESS) throw std::runtime_error("benchmark stream destruction failed");
      } catch(...) {
         if(begin) GPU_EVENT_DESTROY(begin);
         if(end) GPU_EVENT_DESTROY(end);
         solver.release();
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
