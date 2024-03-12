#pragma once

#include <cstdio>

#include "hostMatrix.hpp"
#include "matrix.hpp"

template <typename T, size_t D = 1, size_t I = 0, size_t J = 0, size_t K = 0, size_t L = 0>
class fortMatrix : public hostMatrix<T, D, I, J, K, L> {
private:
   // Set data
   void setData(T *d, size_t i, size_t j, size_t k, size_t l) {
      this->data = d;
      for(int n = 0; n < D; n++) {
         switch(n) {
            case 0: this->dim_size[n] = i; break;
            case 1: this->dim_size[n] = j; break;
            case 2: this->dim_size[n] = k; break;
            case 3: this->dim_size[n] = l; break;
         }
      }
   }

public:
   // Constructors
   fortMatrix() {
   }

   fortMatrix(const hostMatrix<T, D, I, J, K, L> &m) {
      clone(m);
   }

   fortMatrix(const fortMatrix<T, D, I, J, K, L> &m) {
      clone(m);
   }

   fortMatrix(T *d, size_t i, size_t j = 1, size_t k = 1, size_t l = 1) {
      set(d, i, j, k, l);
   }

   // Destructor
   ~fortMatrix() {
      // Prevent parent destructor from attempting to free fortran data
      setData(nullptr, 1, 1, 1, 1);
   }

   // Set
   void set(T *d, size_t i) {
      __MAT_TEST_DIM(1);
      setData(d, i, 1, 1, 1);
   }

   void set(T *d, size_t i, size_t j) {
      __MAT_TEST_DIM(2);
      setData(d, i, j, 1, 1);
   }

   void set(T *d, size_t i, size_t j, size_t k) {
      __MAT_TEST_DIM(3);
      setData(d, i, j, k, 1);
   }

   void set(T *d, size_t i, size_t j, size_t k, size_t l) {
      __MAT_TEST_DIM(4);
      setData(d, i, j, k, l);
   }
};


