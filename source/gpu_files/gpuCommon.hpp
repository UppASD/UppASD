#pragma once

#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#if defined (HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif
#include "real_type.h"
#include "tensor.hpp"
using ParallelizationHelper = GpuParallelizationHelper;

// Class wrapper
class GpuCommon {
public:
   // AddTo parallelization helper
   class AddTo : public ParallelizationHelper::Element {
      real* a;
      const real* b;

   public:
      AddTo(GpuTensor<real, 3>& A, const GpuTensor<real, 3>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int element) {
         a[element] += b[element];
      }
   };

   // Add parallelization helper
   class Add : public ParallelizationHelper::Element {
      real* a;

   public:
      const real* b;
      const real* c;

   public:
      Add(GpuTensor<real, 3>& A, const GpuTensor<real, 3>& B, const GpuTensor<real, 3>& C) {
         a = A.data();
         b = B.data();
         c = C.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] + c[element];
      }
   };

   // Avg parallelization helper
   class Avg : public ParallelizationHelper::Element {
      real* a;
      const real* b;

   public:
      Avg(GpuTensor<real, 3>& A, const GpuTensor<real, 3>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = real(0.5) * (a[element] + b[element]);
      }
   };

   // ScalarMult parallelization helper
   class ScalarMult : public ParallelizationHelper::Element {
      real* a;
      const real* b;
      const real* c;

   public:
      ScalarMult(GpuTensor<real, 3>& A, const GpuTensor<real, 3>& B, const GpuTensor<real, 2>& C) {
         a = A.data();
         b = B.data();
         c = C.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] * c[element / 3];
      }
   };

   // Inv parallelization helper
   class Inv : public ParallelizationHelper::Atom {
      real* a;
      const real* b;

   public:
      Inv(GpuTensor<real, 2>& A, const GpuTensor<real, 2>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int atom) {
         a[atom] = real(1.0) / b[atom];
      }
   };
};

