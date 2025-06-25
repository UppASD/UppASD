// Arkadijs Slobodkins
// Uppsala, 2024

#pragma once

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>

#include "base.cuh"


#define TENSOR_STATIC_ASSERT_DIMENSION() \
   static_assert(sizeof...(Ints) == dim, \
                 "index dimension must be equal to the dimension of the array")


template <typename T, index_t dim>
class Tensor;


template <typename T, index_t dim>
class CudaTensor;


template <typename T>
using Vector = Tensor<T, 1>;

template <typename T>
using Matrix = Tensor<T, 2>;

template <typename T>
using CudaVector = CudaTensor<T, 1>;

template <typename T>
using CudaMatrix = CudaTensor<T, 2>;


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
   using IndexBase<T, dim>::bytes;
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
   // shallow copy of data
   Tensor(const Tensor&) = default;
   Tensor& operator=(const Tensor&) = default;


   void AllocateHost(const Extents<dim>& ext) {
      IndexBase<T, dim>::SetExtents(ext);
      ASSERT_CUDA(cudaMallocHost(&data_, size() * sizeof(T)));
   }


   template <typename... Ints>
   void AllocateHost(Ints... ext) {
      AllocateHost(Extents<dim>{ext...});
   }


   void FreeHost() {
      IndexBase<T, dim>::SetExtents(Extents<dim>{});
      ASSERT_CUDA(cudaFreeHost(data_));
   }


   void set(T* data, const Extents<dim>& ext) {
      IndexBase<T, dim>::SetExtents(ext);
      data_ = data;
   }


   template <typename... Ints>
   void set(T* data, Ints... ext) {
      set(data, Extents<dim>{ext...});
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
   void copy_sync(const Tensor& A) {
      assert(same_extents(*this, A));
      ASSERT_CUDA(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyHostToHost));
   }


   void copy_sync(const CudaTensor<T, dim>& A);


   // Tensor::copy_async(Tensor) is not supported since it is not an asynchronous operation
   void copy_async(const CudaTensor<T, dim>& A, cudaStream_t stream = 0);


   void swap(Tensor<T, dim>& A) noexcept {
      std::swap(data_, A.data_);
      auto tmp_ext = this->extents();
      this->SetExtents(A.extents());
      A.SetExtents(tmp_ext);
   }


   T* begin() {
      return data();
   }


   const T* begin() const {
      return data();
   }


   T* end() {
      return data() + size();
   }


   const T* end() const {
      return data() + size();
   }


   const T* cbegin() const {
      return data();
   }


   const T* cend() const {
      return data() + size();
   }


   void transpose();


   void zeros() {
      std::fill(begin(), end(), T{});
   }


   void resize(const Extents<dim>& ext_new, bool preserve = true) {
      if(preserve) {
         Extents<dim> ext_old = this->extents();
         Tensor<T, dim> old;
         this->swap(old);
         this->AllocateHost(ext_new);

         if constexpr(dim == 1) {
            this->copy_min_dim1(old);
         } else if constexpr(dim == 2) {
            this->copy_min_dim2(old);
         } else if constexpr(dim == 3) {
            this->copy_min_dim3(old);
         } else if constexpr(dim == 4) {
            this->copy_min_dim4(old);
         } else if constexpr(dim == 5) {
            this->copy_min_dim5(old);
         } else if constexpr(dim == 6) {
            this->copy_min_dim6(old);
         }

      } else {
         this->FreeHost();
         this->AllocateHost(ext_new);
      }
   }


private:
   void copy_min_dim1(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         (*this)(i) = A(i);
   }

   void copy_min_dim2(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         for(index_t j = 0; j < std::min(this->extent(1), A.extent(1)); ++j)
            (*this)(i, j) = A(i, j);
   }

   void copy_min_dim3(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         for(index_t j = 0; j < std::min(this->extent(1), A.extent(1)); ++j)
            for(index_t k = 0; k < std::min(this->extent(2), A.extent(2)); ++k)
               (*this)(i, j, k) = A(i, j, k);
   }

   void copy_min_dim4(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         for(index_t j = 0; j < std::min(this->extent(1), A.extent(1)); ++j)
            for(index_t k = 0; k < std::min(this->extent(2), A.extent(2)); ++k)
               for(index_t u = 0; u < std::min(this->extent(3), A.extent(3)); ++u)
                  (*this)(i, j, k, u) = A(i, j, k, u);
   }

   void copy_min_dim5(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         for(index_t j = 0; j < std::min(this->extent(1), A.extent(1)); ++j)
            for(index_t k = 0; k < std::min(this->extent(2), A.extent(2)); ++k)
               for(index_t u = 0; u < std::min(this->extent(3), A.extent(3)); ++u)
                  for(index_t w = 0; w < std::min(this->extent(4), A.extent(4)); ++w)
                     (*this)(i, j, k, u, w) = A(i, j, k, u, w);
   }

   void copy_min_dim6(const Tensor<T, dim>& A) {
      for(index_t i = 0; i < std::min(this->extent(0), A.extent(0)); ++i)
         for(index_t j = 0; j < std::min(this->extent(1), A.extent(1)); ++j)
            for(index_t k = 0; k < std::min(this->extent(2), A.extent(2)); ++k)
               for(index_t u = 0; u < std::min(this->extent(3), A.extent(3)); ++u)
                  for(index_t w = 0; w < std::min(this->extent(4), A.extent(4)); ++w)
                     for(index_t z = 0; z < std::min(this->extent(5), A.extent(5)); ++z)
                        (*this)(i, j, k, u, w, z) = A(i, j, k, u, w, z);
   }

   T* data_{};
};


template <typename T, index_t dim>
void Tensor<T, dim>::transpose() {
   // statically assert since partial template specialization is not allowed
   static_assert(dim == 2, "");

   std::unique_ptr<T[]> TR(new T[size()]);
   for(index_t i = 0; i < extent(0); ++i) {
      for(index_t j = 0; j < extent(1); ++j) {
         TR[i * extent(1) + j] = (*this)(i, j);
      }
   }

   this->SetExtents(extent(1), extent(0));
   std::copy(TR.get(), TR.get() + size(), data());
}


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
   using IndexBase<T, dim>::bytes;
   using IndexBase<T, dim>::extent;
   using IndexBase<T, dim>::extents;


   // can only be created on the host!
   __host__ CudaTensor() : IndexBase<T, dim>{}, data_{} {
   }


   // shallow copy of data
   CudaTensor(const CudaTensor&) = default;
   CudaTensor& operator=(const CudaTensor&) = default;


   // for example if t is Tensor and ct is CudaTensor, can be used as follows:
   // ct.Allocate(t.extents());
   __host__ void Allocate(const Extents<dim>& ext) {
      IndexBase<T, dim>::SetExtents(ext);
      ASSERT_CUDA(cudaMalloc(&data_, size() * sizeof(T)));
   }


   template <typename... Ints>
   __host__ void Allocate(Ints... ext) {
      Allocate(Extents<dim>{ext...});
   }


   __host__ void Free() {
      IndexBase<T, dim>::SetExtents(Extents<dim>{});
      ASSERT_CUDA(cudaFree(data_));
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
      ASSERT_CUDA(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyHostToDevice));
   }


   __host__ void copy_sync(const CudaTensor<T, dim>& A) {
      assert(same_extents(*this, A));
      ASSERT_CUDA(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToDevice));
   }


   __host__ void copy_async(const Tensor<T, dim>& A, cudaStream_t stream = 0) {
      assert(same_extents(*this, A));
      ASSERT_CUDA(cudaMemcpyAsync(
          data(), A.data(), size() * sizeof(T), cudaMemcpyHostToDevice, stream));
   }


   __host__ void copy_async(const CudaTensor<T, dim>& A, cudaStream_t stream = 0) {
      assert(same_extents(*this, A));
      ASSERT_CUDA(cudaMemcpyAsync(
          data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToDevice, stream));
   }


   __host__ void swap(CudaTensor<T, dim>& A) noexcept {
      std::swap(data_, A.data_);
      auto tmp_ext = this->extents();
      this->SetExtents(A.extents());
      A.SetExtents(tmp_ext);
   }


   void zeros() {
      ASSERT_CUDA(cudaMemset(data(), 0, size() * sizeof(T)));
   }


private:
   T* data_{};
};


template <typename T, index_t dim>
void Tensor<T, dim>::copy_sync(const CudaTensor<T, dim>& A) {
   assert(same_extents(*this, A));
   ASSERT_CUDA(cudaMemcpy(data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToHost));
}


template <typename T, index_t dim>
void Tensor<T, dim>::copy_async(const CudaTensor<T, dim>& A, cudaStream_t stream) {
   assert(same_extents(*this, A));
   ASSERT_CUDA(
       cudaMemcpyAsync(data(), A.data(), size() * sizeof(T), cudaMemcpyDeviceToHost, stream));
}
