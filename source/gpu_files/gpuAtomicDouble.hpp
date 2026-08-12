#pragma once

// RCG-06B (F-11): a portable atomicAdd(double*, double).
//
// No minimum CUDA/HIP compute capability is documented or enforced anywhere
// in this codebase -- CMakeLists.txt defaults CMAKE_CUDA_ARCHITECTURES to
// "native" (whatever GPU is present at build time) unless a caller overrides
// it explicitly, and no minimum is checked at configure time. Native
// atomicAdd(double*, double) requires compute capability >= 6.0 (Pascal+),
// which is therefore not guaranteed on every target this codebase can be
// built for. Use the standard CAS-loop double-atomicAdd idiom unconditionally
// instead of the native intrinsic, so energy accumulation is correct
// everywhere. This reuses the same __double_as_longlong/__longlong_as_double
// /atomicCAS(unsigned long long*) idiom
// gpuAdaptiveRuntime.cpp's atomicMaxSelector already relies on, rather than
// inventing a new one.
//
// Shared between gpuAdaptiveRuntime.cpp (the production kernels) and
// test_energy_fp32_accum.cpp (the standalone RCG-06B ENERGY-FP32-ACCUM
// microbenchmark) so the fixture measures the literal accumulator the
// production kernels use, not a reimplementation that could silently drift.
__device__ inline void atomicAddEnergyTerm(double* address, double value) {
   auto* bits = reinterpret_cast<unsigned long long*>(address);
   unsigned long long observed = *bits;
   unsigned long long assumed;
   do {
      assumed = observed;
      const double sum = value +
         __longlong_as_double(static_cast<long long>(assumed));
      observed = atomicCAS(
         bits, assumed,
         static_cast<unsigned long long>(__double_as_longlong(sum)));
   } while(observed != assumed);
}
