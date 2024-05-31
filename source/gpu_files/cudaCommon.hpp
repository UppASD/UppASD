#pragma once

#include <cuda_runtime.h>
#include "cudaParallelizationHelper.hpp"
#include "real_type.h"
#include "tensor.cuh"

// Class wrapper
class CudaCommon {
public:
   // AddTo parallelization helper
   class AddTo : public CudaParallelizationHelper::Element {
      real* a;
      const real* b;

   public:
      AddTo(CudaTensor<real, 3>& A, const CudaTensor<real, 3>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int element) {
         a[element] += b[element];
      }
   };

   // Add parallelization helper
   class Add : public CudaParallelizationHelper::Element {
      real* a;

   public:
      const real* b;
      const real* c;

   public:
      Add(CudaTensor<real, 3>& A, const CudaTensor<real, 3>& B, const CudaTensor<real, 3>& C) {
         a = A.data();
         b = B.data();
         c = C.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] + c[element];
      }
   };

   // Avg parallelization helper
   class Avg : public CudaParallelizationHelper::Element {
      real* a;
      const real* b;

   public:
      Avg(CudaTensor<real, 3>& A, const CudaTensor<real, 3>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = real(0.5) * (a[element] + b[element]);
      }
   };

   // ScalarMult parallelization helper
   class ScalarMult : public CudaParallelizationHelper::Element {
      real* a;
      const real* b;
      const real* c;

   public:
      ScalarMult(CudaTensor<real, 3>& A, const CudaTensor<real, 3>& B, const CudaTensor<real, 2>& C) {
         a = A.data();
         b = B.data();
         c = C.data();
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] * c[element / 3];
      }
   };

   // Inv parallelization helper
   class Inv : public CudaParallelizationHelper::Atom {
      real* a;
      const real* b;

   public:
      Inv(CudaTensor<real, 2>& A, const CudaTensor<real, 2>& B) {
         a = A.data();
         b = B.data();
      }

      __device__ void each(unsigned int atom) {
         a[atom] = real(1.0) / b[atom];
      }
   };
};

