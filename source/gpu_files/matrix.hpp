// Matrix interface class for matrices stored in
// column-major format
//
// Niklas Fejes 2012-2013
//

#pragma once

#include "c_headers.hpp"

// Debug definitions
#ifdef MATRIX_ERROR_INTERRUPT
#include <csignal>
#define __MAT_ERR() std::raise(SIGINT)
#else
#define __MAT_ERR() ((void)0)
#endif

#define __MAT_TEST_DIM(_n_)                                                                       \
   if(D != _n_) {                                                                                 \
      std::printf("Warning: wrong number of indexes (%ld, should be %ld)\n", (long)_n_, (long)D); \
      __MAT_ERR();                                                                                \
   }
#define __MAT_TEST_IOB(_idx_, _i_, _s_)                                                             \
   if(_i_ >= _s_) {                                                                                 \
      std::printf("Warning: index %d out of bounds (i=%ld >= %ld)\n", _idx_, (long)_i_, (long)_s_); \
      __MAT_ERR();                                                                                  \
   }

template <typename T, std::size_t D = 1, std::size_t I = 0, std::size_t J = 0, std::size_t K = 0,
          std::size_t L = 0>
class matrix {
protected:
   // Data fields
   T* data;
   std::size_t dim_size[D];

   // Default constructor
   matrix() {
      data = nullptr;
      for(std::size_t i = 0; i < D; i++) {
         dim_size[i] = 0;
      }
   }

   // Get the offset of an index
   inline std::size_t index(std::size_t i, std::size_t j = 0, std::size_t k = 0, std::size_t l = 0) const {
#ifdef DEBUG
      if(data == nullptr) {
         std::printf("Error: trying to access uninitialized data\n");
         __MAT_ERR();
      }
      if((D == 1 && j + k + l != 0) || (D == 2 && k + l != 0) || (D == 3 && l != 0)) {
         std::printf("Warning: attempting to read with more than %ld indexes from %ld dimension matrix\n",
                     (long)D,
                     (long)D);
         __MAT_ERR();
      }
      for(int n = 0; n < D; n++) {
         switch(n) {
            case 0: __MAT_TEST_IOB(1, i, dim_size[n]); break;
            case 1: __MAT_TEST_IOB(2, j, dim_size[n]); break;
            case 2: __MAT_TEST_IOB(3, k, dim_size[n]); break;
            case 3: __MAT_TEST_IOB(4, l, dim_size[n]); break;
         }
      }
#endif
      std::size_t ind = 0;

      // Causes NVCC to throw false warnings
      // if (D >= 4) ind = (ind + l) * (K != 0 ? K : dim_size[2]);
      // if (D >= 3) ind = (ind + k) * (J != 0 ? J : dim_size[1]);
      // if (D >= 2) ind = (ind + j) * (I != 0 ? I : dim_size[0]);

      // Should be optimized to the same code as above
      for(int n = (int)D - 2; n >= 0; n--) {
         switch(n) {
            case 2: ind = (ind + l) * (K != 0 ? K : dim_size[n]); break;
            case 1: ind = (ind + k) * (J != 0 ? J : dim_size[n]); break;
            case 0: ind = (ind + j) * (I != 0 ? I : dim_size[n]); break;
         }
      }
      return ind + i;
   }

public:
   // Get the size of a dimension
   inline std::size_t dimension_size(std::size_t d) const {
      if(d >= D) {
#ifdef DEBUG
         std::printf("Warning: dimension out of bound (d=%ld, max=%ld)\n", d, D);
         __MAT_ERR();
#endif
         return 1;
      }
      return dim_size[d];
   }

   // Size of data (in bytes)
   inline std::size_t data_size() const {
      std::size_t size = sizeof(T);
      for(std::size_t i = 0; i < D; i++) {
         size *= dim_size[i];
      }
      return size;
   }

   // Size of data (number of elements)
   inline std::size_t size() const {
      std::size_t size = 1;
      for(std::size_t i = 0; i < D; i++) {
         size *= dim_size[i];
      }
      return size;
   }

   // Data member access
   inline const T* get_data() const {
      return data;
   }

   inline T* get_data() {
      return data;
   }

   inline bool has_data() const {
      return data != nullptr;
   }

   // parenthesis-operator
   inline T& operator()(std::size_t i) {
      __MAT_TEST_DIM(1);
      return data[index(i, 0, 0, 0)];
   }

   inline const T& operator()(std::size_t i) const {
      __MAT_TEST_DIM(1);
      return data[index(i, 0, 0, 0)];
   }

   inline T& operator()(std::size_t i, std::size_t j) {
      __MAT_TEST_DIM(2);
      return data[index(i, j, 0, 0)];
   }

   inline const T& operator()(std::size_t i, std::size_t j) const {
      __MAT_TEST_DIM(2);
      return data[index(i, j, 0, 0)];
   }

   inline T& operator()(std::size_t i, std::size_t j, std::size_t k) {
      __MAT_TEST_DIM(3);
      return data[index(i, j, k, 0)];
   }

   inline const T& operator()(std::size_t i, std::size_t j, std::size_t k) const {
      __MAT_TEST_DIM(3);
      return data[index(i, j, k, 0)];
   }

   inline T& operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
      __MAT_TEST_DIM(4);
      return data[index(i, j, k, l)];
   }

   inline const T& operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
      __MAT_TEST_DIM(4);
      return data[index(i, j, k, l)];
   }

   // Allow cast to pointer to type
   operator T*() {
      return data;
   }

   operator T*() const {
      return data;
   }

   operator const T*() const {
      return data;
   }

   // Print matrix info to stdout
   void print_info(const char* name) const {
      std::printf("%s: [data: %p, dims: %ld, elementSize: %ld, fixed size:[%ld,%ld,%ld,%ld]]\n",
                  name,
                  data,
                  (long)D,
                  (long)sizeof(T),
                  (long)I,
                  (long)J,
                  (long)K,
                  (long)L);
      for(std::size_t i = 0; i < D; i++) {
         std::printf("    dim %ld: %ld\n", (long)i + 1, (long)dim_size[i]);
      }
   }
};

