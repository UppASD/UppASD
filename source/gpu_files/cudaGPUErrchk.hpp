#pragma once

#include <cuda_runtime.h>

#include "c_headers.hpp"

#define gpuErrchk(ans) \
   { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
   if(code != cudaSuccess) {
      std::fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if(abort) {
         std::exit(code);
      }
   }
}

