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

// CGP-00: the deterministic block-local tree reduction RCG-09B introduced
// (one partial per launch block, fixed binary-tree accumulation order),
// generalized over the accumulator type instead of hardcoded to `double`.
// Production instantiates this at `Acc = double` only, so
// `reduceAdaptiveEnergyBlockT<double>` is byte-for-byte the same codegen
// RCG-09B shipped; the CGP-00 hierarchical-precision evidence fixture
// (test_energy_hierarchical_precision.cpp) additionally instantiates
// `Acc = float`, sharing this literal primitive rather than reimplementing
// it, for the same anti-drift reason `atomicAddEnergyTerm` above is shared.
template <typename Acc>
__device__ inline void reduceAdaptiveEnergyBlockT(
   Acc value, Acc* partial, Acc* shared) {
   shared[threadIdx.x] = value;
   __syncthreads();
   for(unsigned int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
      if(threadIdx.x < stride)
         shared[threadIdx.x] += shared[threadIdx.x + stride];
      __syncthreads();
   }
   if(threadIdx.x == 0) partial[blockIdx.x] = shared[0];
}
