#pragma once

#include "hostMatrix.hpp"
#include "matrix.hpp"
#include "real_type.h"

template <typename T, usd_int D = 1, usd_int I = 0, usd_int J = 0, usd_int K = 0, usd_int L = 0>
class fortMatrix : public hostMatrix<T, D, I, J, K, L> {
private:
   // Set data
   void setData(T *d, usd_int i, usd_int j, usd_int k, usd_int l) {
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

   fortMatrix(T *d, usd_int i, usd_int j = 1, usd_int k = 1, usd_int l = 1) {
      set(d, i, j, k, l);
   }

   // Destructor
   ~fortMatrix() {
      // Prevent parent destructor from attempting to free fortran data
      setData(nullptr, 1, 1, 1, 1);
   }

   // Set
   void set(T *d, usd_int i) {
      __MAT_TEST_DIM(1);
      setData(d, i, 1, 1, 1);
   }

   void set(T *d, usd_int i, usd_int j) {
      __MAT_TEST_DIM(2);
      setData(d, i, j, 1, 1);
   }

   void set(T *d, usd_int i, usd_int j, usd_int k) {
      __MAT_TEST_DIM(3);
      setData(d, i, j, k, 1);
   }

   void set(T *d, usd_int i, usd_int j, usd_int k, usd_int l) {
      __MAT_TEST_DIM(4);
      setData(d, i, j, k, l);
   }
};

