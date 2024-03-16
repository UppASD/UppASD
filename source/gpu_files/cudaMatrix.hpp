// Matrix class for device matrices stored in
// column-major format
//
// Niklas Fejes 2012-2013
//

#pragma once

#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "hostMatrix.hpp"
#include "matrix.hpp"
#include "real_type.h"

template <typename T, std::size_t D = 1, std::size_t I = 0, std::size_t J = 0, std::size_t K = 0, std::size_t L = 0>
class cudaMatrix : public matrix<T, D, I, J, K, L> {
private:
   // Initiate
   bool init(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
      // Calculate size
      std::size_t size = i * j * k * l * sizeof(T);

      // Deallocate data if already allocated
      if(this->data != nullptr) {
         cudaFree(this->data);
         this->data = nullptr;
      }

      // Allocate new  memory
      cudaError e = cudaMalloc((void**)&this->data, size);

      // Error?
      if(e != cudaSuccess) {
#ifdef DEBUG
         std::printf(
             "cudaMatrix::initiate: cudaMalloc returned %d (%s)"
             " when trying to allocate %ld bytes\n",
             e,
             cudaGetErrorString(e),
             size);
         __MAT_ERR();
#endif
         return false;
      }

      // Set dimensions
      for(int n = 0; n < D; n++) {
         switch(n) {
            case 0: this->dim_size[n] = i; break;
            case 1: this->dim_size[n] = j; break;
            case 2: this->dim_size[n] = k; break;
            case 3: this->dim_size[n] = l; break;
         }
      }
      return true;
   }

public:
   // Constructors
   cudaMatrix() {
   }

   cudaMatrix(const hostMatrix<T, D, I, J, K, L>& m) {
      clone(m);
   }

   //	cudaMatrix(std::size_t i, std::size_t j = 1, std::size_t k = 1, std::size_t l = 1) {
   //		initiate(i, j, k, l);
   //	}

   // Destructor
   ~cudaMatrix() {
      free();
   }

   // Initiate
   bool initiate(std::size_t i) {
      __MAT_TEST_DIM(1);
      return init(i, 1, 1, 1);
   }

   bool initiate(std::size_t i, std::size_t j) {
      __MAT_TEST_DIM(2);
      return init(i, j, 1, 1);
   }

   bool initiate(std::size_t i, std::size_t j, std::size_t k) {
      __MAT_TEST_DIM(3);
      return init(i, j, k, 1);
   }

   bool initiate(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
      __MAT_TEST_DIM(4);
      return init(i, j, k, l);
   }

   bool initiate(const matrix<T, D, I, J, K, L>& m) {
      return init((D < 1 ? 1 : m.dimension_size(0)),
                  (D < 2 ? 1 : m.dimension_size(1)),
                  (D < 3 ? 1 : m.dimension_size(2)),
                  (D < 4 ? 1 : m.dimension_size(3)));
   }

   // Free
   void free() {
      if(this->data != nullptr) {
         cudaFree(this->data);
      }
      this->data = nullptr;
      for(std::size_t i = 0; i < D; i++) {
         this->dim_size[i] = 0;
      }
   }

   // Swap pointers
   inline void swap(cudaMatrix<T, D, I, J, K, L>& m) {
#ifdef DEBUG
      for(int n = 0; n < D; n++) {
         if(this->dim_size[n] != m.dim_size[n]) {
            std::printf("Warning: swapping pointers between matrices with different sizes\n");
            __MAT_ERR();
            break;
         }
      }
#endif
      T* tmp = this->data;
      this->data = m.data;
      m.data = tmp;
   }

   // Copy data
   inline bool memcopy(const hostMatrix<T, D, I, J, K, L>& m) {
#ifdef DEBUG
      for(int n = 0; n < D; n++) {
         if(this->dim_size[n] != m.dimension_size(n)) {
            std::printf("Warning: copying data between matrices with different sizes\n");
            __MAT_ERR();
            break;
         }
      }
#endif

      // Copy data
      cudaError_t e = cudaMemcpy(this->data, m.get_data(), this->data_size(), cudaMemcpyHostToDevice);

      // Error?
      if(e != cudaSuccess) {
         std::printf("cudaMatrix::memcopy: cudaMemcpy returned %d (%s)\n", e, cudaGetErrorString(e));
      }

      return (e == cudaSuccess);
   }

   // Is this really ok? Will the compiler know which one to use if this version is defaulting on second arg
   // Thomas Nystrand
   inline bool memcopy(const cudaMatrix<T, D, I, J, K, L>& m, cudaStream_t stream = 0) {
#ifdef DEBUG
      for(int n = 0; n < D; n++) {
         if(this->dim_size[n] != m.dim_size[n]) {
            std::printf("Warning: copying data between device matrices with different sizes\n");
            __MAT_ERR();
            break;
         }
      }
#endif

      // Copy data
      cudaError_t e
          = cudaMemcpyAsync(this->data, m.data, this->data_size(), cudaMemcpyDeviceToDevice, stream);

      // Error?
      if(e != cudaSuccess) {
         std::printf("cudaMatrix::memcopy: cudaMemcpyAsync returned %d (%s)\n", e, cudaGetErrorString(e));
      }

      return (e == cudaSuccess);
   }

   inline bool memcopyTo(hostMatrix<T, D, I, J, K, L>& m) const {
#ifdef DEBUG
      for(int n = 0; n < D; n++) {
         if(this->dim_size[n] != m.dimension_size(n)) {
            std::printf("Warning: copying data between matrices with different sizes\n");
            __MAT_ERR();
            break;
         }
      }
#endif

      // Copy data
      cudaError_t e = cudaMemcpyAsync(m.get_data(), this->data, this->data_size(), cudaMemcpyDeviceToHost);

      // Error?
      if(e != cudaSuccess) {
         std::printf("cudaMatrix::memcopyTo: cudaMemcpyAsync returned %d (%s)\n", e, cudaGetErrorString(e));
      }

      return (e == cudaSuccess);
   }

   // Set to zero
   inline void zero() {
      cudaMemset(this->data, 0, this->data_size());
   }

   // Clone
   bool clone(const hostMatrix<T, D, I, J, K, L>& m) {
      // Zero matrix?
      if(m.data_size() == 0) {
         free();
         for(int n = 0; n < D; n++) {
            this->dim_size[n] = m.dimension_size(n);
         }
         return true;
      }

      // Can we keep old data?
      if(this->data && m.size() == this->size()) {
         // Update dimensions
         for(int i = 0; i < D; i++) {
            this->dim_size[i] = m.dimension_size(i);
         }
         // Copy data
         return memcopy(m);
      }

      // Free old
      free();

      // Initiate new
      if(!init(m.dimension_size(0),
               (D < 2 ? 1 : m.dimension_size(1)),
               (D < 3 ? 1 : m.dimension_size(2)),
               (D < 4 ? 1 : m.dimension_size(3)))) {
         return false;
      }

      // Copy the data
      return memcopy(m);
   }

   // Read
   void read(const T* d) {
      // Get memory size
      std::size_t size = this->data_size();

      // Invalid copy?
      if(d == nullptr || this->data == nullptr) {
         std::printf("cudaMatrix::read: Invalid read attempt: (%p,%p,%ld)\n", this->data, d, size);
         return;
      }

      // Copy
      cudaError_t e = cudaMemcpy(this->data, d, size, cudaMemcpyHostToDevice);

      // Error?
      if(e != cudaSuccess) {
         std::printf("cudaMatrix::read: cudaMemcpy returned %d (%s)\n", e, cudaGetErrorString(e));
      }
   }

   void write(T* d) const {
      // Get memory size
      std::size_t size = this->data_size();

      // Invalid copy?
      if(d == nullptr || this->data == nullptr) {
         std::printf("cudaMatrix::write: Invalid write attempt: (%p,%p,%ld)\n", this->data, d, size);
         return;
      }

      cudaError_t e = cudaMemcpy(d, this->data, size, cudaMemcpyDeviceToHost);

      if(e != cudaSuccess) {
         std::printf("cudaMatrix::write: cudaMemcpy returned %d (%s)\n", e, cudaGetErrorString(e));
      }
   }

   void writeAsync(T* d, cudaStream_t stream = 0) const {
      // Get memory size
      std::size_t size = this->data_size();

      // Invalid copy?
      if(d == nullptr || this->data == nullptr) {
         std::printf("cudaMatrix::write: Invalid write attempt: (%p,%p,%ld)\n", this->data, d, size);
         return;
      }

      cudaError_t e = cudaMemcpyAsync(d, this->data, size, cudaMemcpyDeviceToHost, stream);

      if(e != cudaSuccess) {
         std::printf("cudaMatrix::write: cudaMemcpy returned %d (%s)\n", e, cudaGetErrorString(e));
      }
   }
};

