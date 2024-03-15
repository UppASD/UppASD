#pragma once

/*  Helper file for setting up the grid for kernels.
 *  The maximum size of the x dimension is limited, so
 *  if there are more atoms than that in the system
 *  we need to split up the indexing over x and y.
 *
 */

#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "real_type.h"

template <usd_int threads, bool big>
class GridHelper {
   usd_int maxGridSize1;

public:
   // Types for template parameters
   enum Type { Small = false, Big = true };

   GridHelper() {
      int device;
      struct cudaDeviceProp prop;
      cudaGetDevice(&device);
      cudaGetDeviceProperties(&prop, device);
      maxGridSize1 = prop.maxGridSize[0];
   }

   // Grid dimensions
   inline bool dim1d(dim3 *block, dim3 *grid, usd_int X) {
      if(big) {
         // Total number of blocks
         usd_int x = (X + threads - 1) / threads;
         usd_int y = 1;

         // The block's first dimension can't be larger than maxGridSize1
         while(x > maxGridSize1) {
            y *= 2;
            x = (x + 1) / 2;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(x, y, 1);

         return true;
      } else {
         // Total number of blocks
         usd_int x = (X + threads - 1) / threads;

         // Small if x <= maxGridSize1
         if(x > maxGridSize1) {
            std::fprintf(stderr, "Error: Too many atoms in system, bigger grid must be used!\n");
            std::exit(EXIT_FAILURE);
            return false;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(x, 1, 1);

         return true;
      }
   }

   inline bool dim2d(dim3 *block, dim3 *grid, usd_int X, usd_int Y) {
      if(big) {
         // Total number of blocks
         usd_int x = (X * Y + threads - 1) / threads;
         usd_int y = 1;

         // The block dimension can't be larger than maxGridSize1
         while(x > maxGridSize1) {
            y *= 2;
            x = (x + 1) / 2;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(x, y, 1);

         return true;
      } else {
         // Total number of blocks
         usd_int x = (X + threads - 1) / threads;

         // Small if x <= maxGridSize1
         if(x > maxGridSize1) {
            std::fprintf(stderr, "Error: Too many atoms in system, bigger grid must be used!\n");
            std::exit(EXIT_FAILURE);
            return false;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(x, Y, 1);

         return true;
      }
   }

   inline bool dim3d(dim3 *block, dim3 *grid, usd_int X, usd_int Y, usd_int Z) {
      if(big) {
         // Total number of blocks
         usd_int x = (X * Y * Z + threads - 1) / threads;
         usd_int y = 1;

         // The block dimension can't be larger than maxGridSize1
         while(x > maxGridSize1) {
            y *= 2;
            x = (x + 1) / 2;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(x, y, 1);

         return true;
      } else {
         // Total number of blocks
         usd_int xy = (X * Y + threads - 1) / threads;

         // Small if y <= maxGridSize1
         if(xy > maxGridSize1) {
            std::fprintf(stderr, "Error: Too many atoms in system, bigger grid must be used!\n");
            std::exit(EXIT_FAILURE);
            return false;
         }

         // Create dims
         *block = dim3(threads, 1, 1);
         *grid = dim3(xy, Z, 1);

         return true;
      }
   }

   // 1D-index
   inline static __device__ bool index1d(usd_int *x, usd_int X) {
      if(big) {
         *x = (blockIdx.x + gridDim.x * blockIdx.y) * threads + threadIdx.x;
         return (*x < X);
      } else {
         *x = blockIdx.x * threads + threadIdx.x;
         return (*x < X);
      }
   }

   // 2D-index
   inline static __device__ bool index2d(usd_int *x, usd_int *y, usd_int X, usd_int Y) {
      if(big) {
         usd_int xy = (blockIdx.x + gridDim.x * blockIdx.y) * threads + threadIdx.x;
         *x = xy % X;
         *y = xy / X;

         // Index x is guaranteed to be less than X (due xy % X)
         return (*y < Y);
      } else {
         *x = blockIdx.x * threads + threadIdx.x;
         *y = blockIdx.y;

         // Index y is guaranteed to be less than Y
         return (*x < X);
      }
   }

   // 3D-index
   inline static __device__ bool index3d(usd_int *x, usd_int *y, usd_int *z, usd_int X, usd_int Y,
                                         usd_int Z) {
      if(big) {
         usd_int xyz = (blockIdx.x + gridDim.x * blockIdx.y) * threads + threadIdx.x;
         usd_int yz = xyz / X;
         *x = xyz % X;
         *y = yz % Y;
         *z = yz / Y;

         // Index x and y is guaranteed to be less than X and Y
         return (*z < Z);
      } else {
         usd_int xy = blockIdx.x * threads + threadIdx.x;
         *x = xy % X;
         *y = xy / X;
         *z = blockIdx.y;

         // Index y and z is guaranteed to be less than Y and Z
         return (*y < Y);
      }
   }
};

