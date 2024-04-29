// Arkadijs Slobodkins
// Uppsala, 2024

#pragma once

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <type_traits>

#include "common.cuh"


#define TENSOR_STATIC_ASSERT_DIMENSION() \
   static_assert(sizeof...(Ints) == dim, "index dimension must be equal to the dimension of the array")


template <typename T, index_t dim>
class Tensor;


template <typename T, index_t dim>
class CudaTensor;


template <typename TensorType1, typename TensorType2>
__host__ __device__ bool same_extents(const TensorType1& A, const TensorType2& B) {
   assert(A.dimension() == B.dimension());
   for(index_t d = 0; d < A.dimension(); ++d) {
      if(A.extent(d) != B.extent(d)) {
         return false;
      }
   }
   return true;
}


template <typename T, index_t dim>
class Tensor : private IndexBase<T, dim> {
public:
   using size_type = index_t;
   using value_type = T;

   using IndexBase<T, dim>::dimension;
   using IndexBase<T, dim>::size;
   using IndexBase<T, dim>::extent;
   using IndexBase<T, dim>::extents;

   ////////////////////////////////////////////////////////////////////////////////////////////////
   template <typename... Ints>
   Tensor(T* data, Ints... ext) : IndexBase<T, dim>{ext...},
                                  data_{data} {
   }


   Tensor(T* data, const Extents<dim>& ext) : IndexBase<T, dim>{ext}, data_{data} {
   }


   Tensor() = default;


   // deleted since we currently don't need these operations
   Tensor(const Tensor&) = delete;
   Tensor& operator=(const Tensor&) = delete;


   template <typename... Ints>
   void set(T* data, Ints... ext) {
      data_ = data;
      IndexBase<T, dim>::SetExtents(ext...);
   }


   void set(T* data, const Extents<dim>& ext) {
      data_ = data;
      IndexBase<T, dim>::SetExtents(ext);
   }


   T* data() {
      return size() == 0 ? nullptr : data_;
   }


   const T* data() const {
      return size() == 0 ? nullptr : data_;
   }


   bool empty() const {
      return !size();
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////
   template <typename... Ints>
   T& operator()(Ints... indexes) {
      TENSOR_STATIC_ASSERT_DIMENSION();
      return data_[IndexBase<T, dim>::index_of(indexes...)];
   }


   template <typename... Ints>
   const T& operator()(Ints... indexes) const {
      TENSOR_STATIC_ASSERT_DIMENSION();
      return data_[IndexBase<T, dim>::index_of(indexes...)];
   }


   T& operator[](index_t i) {
      assert(i > -1 && i < size());
      return data_[i];
   }


   const T& operator[](index_t i) const {
      assert(i > -1 && i < size());
      return data_[i];
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////
   void copy_sync(const CudaTensor<T, dim>& A);


private:
   T* data_{};
};


template <typename T>
std::ostream& operator<<(std::ostream& os, const Tensor<T, 1>& A) {
   os << std::fixed << std::setprecision(7);
   for(index_t i = 0; i < A.extent(0); ++i) {
      std::cout << A(i) << std::endl;
   }
   return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const Tensor<T, 2>& A) {
   os << std::fixed << std::setprecision(7);

   for(index_t i = 0; i < A.extent(0); ++i) {
      for(index_t j = 0; j < A.extent(1); ++j) {
         os << A(i, j) << " ";
      }
      os << std::endl;
   }
   return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const Tensor<T, 3>& A) {
   os << std::fixed << std::setprecision(7);

   for(index_t i = 0; i < A.extent(0); ++i) {
      os << "A(" << i << ", :, :) = " << std::endl;
      for(index_t j = 0; j < A.extent(1); ++j) {
         os << "  ";
         for(index_t k = 0; k < A.extent(2); ++k) {
            os << A(i, j, k) << " ";
         }
         os << std::endl;
      }
      os << std::endl;
   }
   return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, const Tensor<T, 4>& A) {
   os << std::fixed << std::setprecision(7);

   for(index_t i = 0; i < A.extent(0); ++i) {
      for(index_t j = 0; j < A.extent(1); ++j) {
         os << "A(" << i << ", " << j << ", :, :) = " << std::endl;
         for(index_t k = 0; k < A.extent(2); ++k) {
            os << "  ";
            for(index_t l = 0; l < A.extent(3); ++l) {
               os << A(i, j, k, l) << " ";
            }
            os << std::endl;
         }
         os << std::endl;
      }
   }
   return os;
}


template <typename T, index_t dim>
class CudaTensor : private IndexBase<T, dim> {
public:
   using size_type = index_t;
   using value_type = T;

   using IndexBase<T, dim>::dimension;
   using IndexBase<T, dim>::size;
   using IndexBase<T, dim>::extent;
   using IndexBase<T, dim>::extents;


   // can only be created on the host!
   __host__ CudaTensor() : IndexBase<T, dim>{}, data_{} {
   }


   // shallow copy of data
   CudaTensor(const CudaTensor&) = default;
   CudaTensor& operator=(const CudaTensor&) = delete;


   template <typename... Ints>
   __host__ void Allocate(Ints... ext) {
      IndexBase<T, dim>::SetExtents(ext...);
      ASSERT_CUDA_SUCCESS(cudaMalloc(&data_, size() * sizeof(T)));
   }


   // for example if t is Tensor and ct is CudaTensor, can be used as follows: ct.Allocate(t.extents());
   __host__ void Allocate(const Extents<dim>& ext) {
      IndexBase<T, dim>::SetExtents(ext);
      ASSERT_CUDA_SUCCESS(cudaMalloc(&data_, size() * sizeof(T)));
   }


   __host__ void Free() {
      ASSERT_CUDA_SUCCESS(cudaFree(data_));
      IndexBase<T, dim>::SetExtents(Extents<dim>{});
   }


   __host__ __device__ T* data() {
      return size() == 0 ? nullptr : data_;
   }


   __host__ __device__ const T* data() const {
      return size() == 0 ? nullptr : data_;
   }


   __host__ __device__ bool empty() const {
      return !size();
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////
   template <typename... Ints>
   __device__ T& operator()(Ints... indexes) {
      TENSOR_STATIC_ASSERT_DIMENSION();
      return data_[IndexBase<T, dim>::index_of(indexes...)];
   }


   template <typename... Ints>
   __device__ const T& operator()(Ints... indexes) const {
      TENSOR_STATIC_ASSERT_DIMENSION();
      return data_[IndexBase<T, dim>::index_of(indexes...)];
   }


   __device__ T& operator[](index_t i) {
      assert(i > -1 && i < size());
      return data_[i];
   }


   __device__ const T& operator[](index_t i) const {
      assert(i > -1 && i < size());
      return data_[i];
   }


   ////////////////////////////////////////////////////////////////////////////////////////////////
   __host__ void copy_sync(const Tensor<T, dim>& A) {
      assert(same_extents(*this, A));
      ASSERT_CUDA_SUCCESS(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyHostToDevice));
   }


   __host__ void copy_sync(const CudaTensor<T, dim>& A) {
      assert(same_extents(*this, A));
      ASSERT_CUDA_SUCCESS(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToDevice));
   }


private:
   T* data_{};
};


template <typename T, index_t dim>
void Tensor<T, dim>::copy_sync(const CudaTensor<T, dim>& A) {
   assert(same_extents(*this, A));
   ASSERT_CUDA_SUCCESS(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToHost));
}
