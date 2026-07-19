// Arkadijs Slobodkins
// Uppsala, 2024

#pragma once


#if defined(HIP_V)
  #include <hip/hip_runtime.h>
#elif defined(CUDA_V)
  #include <cuda_runtime.h>
#endif
#include "gpu_wrappers.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>


using index_t = long int;


// Bounds checks are valuable for host-side tensor manipulation, but device
// indexing is on the innermost kernel path.  Do not emit device asserts there.
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  #define GPU_HOST_ASSERT(condition) ((void)0)
#else
  #define GPU_HOST_ASSERT(condition) assert(condition)
#endif


#define ASSERT_GPU(gpuCall)                            \
   {                                                     \
      GPU_ERROR_T error = gpuCall;                      \
      if(error != GPU_SUCCESS) {                         \
          throw std::runtime_error(GPU_GET_ERROR_STRING(error)); \
      }                                                  \
   }

#define ASSERT_CURAND(randCall)                         \
   { const auto status = randCall;                       \
     if(status != GPU_RAND_STATUS_SUCCESS)               \
       throw std::runtime_error("GPU random-generator error: " + std::to_string(static_cast<int>(status))); }


template <index_t dim>
struct Extents {
   static_assert(dim > 0, "");

private:
   index_t x[dim]{};

public:
   Extents() = default;
   Extents(const Extents&) = default;
   Extents& operator=(const Extents&) = default;

   template <typename... Ints>
   __host__ __device__ Extents(Ints... ext) : x{static_cast<index_t>(ext)...} {
      static_assert(static_cast<index_t>(sizeof...(ext)) == dim, "");
   }

   __host__ __device__ index_t& operator[](index_t d) {
      GPU_HOST_ASSERT(d > -1 && d < dim);
      return x[d];
   }

   __host__ __device__ const index_t& operator[](index_t d) const {
      GPU_HOST_ASSERT(d > -1 && d < dim);
      return x[d];
   }

   __host__ __device__ index_t size() const {
      index_t p = 1;
      for(index_t d = 0; d < dim; ++d) {
         p *= x[d];
      }
      return p;
   }

   friend std::ostream& operator<<(std::ostream& os, const Extents<dim>& ext) {
      for(index_t d = 0; d < dim; ++d) {
         os << ext[d] << std::endl;
      }
      return os;
   }
};


template <typename T, index_t dim>
class IndexBase {
private:
   Extents<dim> ext_{};

protected:
   template <typename... Ints>
   __host__ __device__ void SetExtents(Ints... ext) {
      ext_ = {ext...};
   }

   __host__ __device__ void SetExtents(const Extents<dim>& ext) {
      ext_ = ext;
   }

public:
   template <typename... Ints>
   __host__ __device__ IndexBase(Ints... ext) {
      SetExtents(ext...);
   }


   __host__ __device__ IndexBase(const Extents<dim>& ext) {
      SetExtents(ext);
   }


   IndexBase() = default;
   IndexBase(const IndexBase&) = default;
   IndexBase& operator=(const IndexBase&) = default;


   ////////////////////////////////////////////////////////////////////////////////////////////////
   __host__ __device__ static constexpr index_t dimension() {
      return dim;
   }


   __host__ __device__ index_t size() const {
      return ext_.size();
   }


   __host__ __device__ index_t bytes() const {
      return ext_.size() * sizeof(T);
   }


   __host__ __device__ index_t extent(index_t d) const {
      return ext_[d];
   }


   __host__ __device__ Extents<dim> extents() const {
      return ext_;
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////
   __host__ __device__ index_t index_of(index_t i) const {
      GPU_HOST_ASSERT(valid_index(i, 0));
      return i;
   }


   __host__ __device__ index_t index_of(index_t i, index_t j) const {
      GPU_HOST_ASSERT(valid_index(j, 1));
      return index_of(i) + j * ext_[0];
   }


   __host__ __device__ index_t index_of(index_t i, index_t j, index_t k) const {
      GPU_HOST_ASSERT(valid_index(k, 2));
      return index_of(i, j) + k * ext_[0] * ext_[1];
   }


   __host__ __device__ index_t index_of(index_t i, index_t j, index_t k, index_t l) const {
      GPU_HOST_ASSERT(valid_index(l, 3));
      return index_of(i, j, k) + l * ext_[0] * ext_[1] * ext_[2];
   }


   __host__ __device__ index_t index_of(index_t i, index_t j, index_t k, index_t l,
                                        index_t m) const {
      GPU_HOST_ASSERT(valid_index(m, 4));
      return index_of(i, j, k, l) + m * ext_[0] * ext_[1] * ext_[2] * ext_[3];
   }


   __host__ __device__ index_t index_of(index_t i, index_t j, index_t k, index_t l, index_t m,
                                        index_t n) const {
      GPU_HOST_ASSERT(valid_index(n, 5));
      return index_of(i, j, k, l, m) + n * ext_[0] * ext_[1] * ext_[2] * ext_[3] * ext_[4];
   }


   __host__ __device__ bool valid_index(index_t i, index_t d) const {
      return i > -1 && i < ext_[d];
   }
};
