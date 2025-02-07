#pragma once

#include <cuda_runtime.h>
#include <curand.h>
#include "cudaParallelizationHelper.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"

// CUDA thermfield class
// Optimized for simulations with constant temperatures and timestep.

class CudaThermfield {
private:
   // Generator
   curandGenerator_t gen;

   // Initiation flags
   bool dataInitiated;
   bool constantsInitiated;

   // Data
   CudaTensor<real, 3> field;
   CudaTensor<real, 1> sigmaFactor;  // = sqrt(Dp*temperature(i))

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   CudaParallelizationHelper& parallel;

public:
   // Parallelization helpers
   class SetupSigmaFactor;
   class SetupField;

   // Constructor / destructor
   CudaThermfield();
   ~CudaThermfield();

   // Initiate
   bool initiate(std::size_t N, std::size_t M, curandRngType_t rngType = CURAND_RNG_PSEUDO_DEFAULT,
                 unsigned long long seed = 0);
   bool initiateConstants(const Tensor<real, 1>& temperature, real timestep, real gamma, real k_bolt,
                          real mub, real damping);

   // Initiated?
   inline bool initiated() {
      return dataInitiated && constantsInitiated;
   }

   // Get field
   inline const CudaTensor<real, 3>& getField() {
      return field;
   }

   // Randomize
   void randomize(const CudaTensor<real, 2>& mmom);
};

