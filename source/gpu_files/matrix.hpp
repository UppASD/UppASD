// Matrix interface class for matrices stored in
// column-major format
//
// Niklas Fejes 2012-2013
//

#pragma once

#include "c_headers.hpp"
#include "real_type.h"

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

template <typename T, usd_int D = 1, usd_int I = 0, usd_int J = 0, usd_int K = 0,
          usd_int L = 0>
class matrix {
protected:
   // Data fields
   T* data;
   usd_int dim_size[D];

   // Default constructor
   matrix() {
      data = nullptr;
      for(usd_int i = 0; i < D; i++) {
         dim_size[i] = 0;
      }
   }

   // Get the offset of an index
   inline usd_int index(usd_int i, usd_int j = 0, usd_int k = 0, usd_int l = 0) const {
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
      usd_int ind = 0;

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
   inline usd_int dimension_size(usd_int d) const {
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
   inline usd_int data_size() const {
      usd_int size = sizeof(T);
      for(usd_int i = 0; i < D; i++) {
         size *= dim_size[i];
      }
      return size;
   }

   // Size of data (number of elements)
   inline usd_int size() const {
      usd_int size = 1;
      for(usd_int i = 0; i < D; i++) {
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
   inline T& operator()(usd_int i) {
      __MAT_TEST_DIM(1);
      return data[index(i, 0, 0, 0)];
   }

   inline const T& operator()(usd_int i) const {
      __MAT_TEST_DIM(1);
      return data[index(i, 0, 0, 0)];
   }

   inline T& operator()(usd_int i, usd_int j) {
      __MAT_TEST_DIM(2);
      return data[index(i, j, 0, 0)];
   }

   inline const T& operator()(usd_int i, usd_int j) const {
      __MAT_TEST_DIM(2);
      return data[index(i, j, 0, 0)];
   }

   inline T& operator()(usd_int i, usd_int j, usd_int k) {
      __MAT_TEST_DIM(3);
      return data[index(i, j, k, 0)];
   }

   inline const T& operator()(usd_int i, usd_int j, usd_int k) const {
      __MAT_TEST_DIM(3);
      return data[index(i, j, k, 0)];
   }

   inline T& operator()(usd_int i, usd_int j, usd_int k, usd_int l) {
      __MAT_TEST_DIM(4);
      return data[index(i, j, k, l)];
   }

   inline const T& operator()(usd_int i, usd_int j, usd_int k, usd_int l) const {
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
      for(usd_int i = 0; i < D; i++) {
         std::printf("    dim %ld: %ld\n", (long)i + 1, (long)dim_size[i]);
      }
   }
};

