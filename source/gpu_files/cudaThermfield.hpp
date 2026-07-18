#pragma once

#include <cuda_runtime.h>
#include <curand.h>
#include "gpuParallelizationHelper.hpp"
#include "tensor.hpp"
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
   GpuTensor<real, 3> field;
   GpuVector<real> randomValues;
   GpuTensor<real, 1> sigmaFactor;  // = sqrt(Dp*temperature(i))

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
  GpuParallelizationHelper& parallel;

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
   void resetConstants(const Tensor<real, 1>& temperature, real timestep, real gamma, real k_bolt,
                          real mub, real damping);

   // Initiated?
   inline bool initiated() {
      return dataInitiated && constantsInitiated;
   }

   // Get field
   inline const GpuTensor<real, 3>& getField() {
      return field;
   }

   // Randomize
   void randomize(const GpuTensor<real, 2>& mmom);
};
