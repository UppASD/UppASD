#pragma once

#include <cuda_runtime.h>

#include "cudaMatrix.hpp"
#include "cudaParallelizationHelper.hpp"
#include "real_type.h"

// Class wrapper
class CudaCommon {
public:
   // AddTo parallelization helper
   class AddTo : public CudaParallelizationHelper::Element {
      real *a;
      const real *b;

   public:
      AddTo(cudaMatrix<real, 3, 3>& A, const cudaMatrix<real, 3, 3>& B) {
         a = A;
         b = B;
      }

      __device__ void each(unsigned int element) {
         a[element] += b[element];
      }
   };

   // Add parallelization helper
   class Add : public CudaParallelizationHelper::Element {
      real *a;

   public:
      const real *b;
      const real *c;

   public:
      Add(cudaMatrix<real, 3, 3>& A, const cudaMatrix<real, 3, 3>& B, const cudaMatrix<real, 3, 3>& C) {
         a = A;
         b = B;
         c = C;
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] + c[element];
      }
   };

   // Avg parallelization helper
   class Avg : public CudaParallelizationHelper::Element {
      real *a;
      const real *b;

   public:
      Avg(cudaMatrix<real, 3, 3>& A, const cudaMatrix<real, 3, 3>& B) {
         a = A;
         b = B;
      }

      __device__ void each(unsigned int element) {
         a[element] = real(0.5) * (a[element] + b[element]);
      }
   };

   // ScalarMult parallelization helper
   class ScalarMult : public CudaParallelizationHelper::Element {
      real *a;
      const real *b;
      const real *c;

   public:
      ScalarMult(cudaMatrix<real, 3, 3>& A, const cudaMatrix<real, 3, 3>& B, const cudaMatrix<real, 2>& C) {
         a = A;
         b = B;
         c = C;
      }

      __device__ void each(unsigned int element) {
         a[element] = b[element] * c[element / 3];
      }
   };

   // Inv parallelization helper
   class Inv : public CudaParallelizationHelper::Atom {
      real *a;
      const real *b;

   public:
      Inv(cudaMatrix<real, 2> &A, const cudaMatrix<real, 2> &B) {
         a = A;
         b = B;
      }

      __device__ void each(unsigned int atom) {
         a[atom] = real(1.0) / b[atom];
      }
   };
};

