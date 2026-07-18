#pragma once

#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"

class GpuThermfield {
private:
   GPU_RAND_GENERATOR_T gen;
   bool dataInitiated = false;
   bool constantsInitiated = false;
   GpuTensor<real, 3> field;
   GpuVector<real> randomValues;
   GpuTensor<real, 1> sigmaFactor;
   StopwatchDeviceSync stopwatch;
   GpuParallelizationHelper& parallel;
public:
   class SetupSigmaFactor;
   class SetupField;
   GpuThermfield();
   ~GpuThermfield();
   bool initiate(std::size_t N, std::size_t M,
                 GPU_RAND_RNGTYPE_T rngType = GPU_RAND_RNG_PSEUDO_DEFAULT,
                 unsigned long long seed = 0);
   bool initiateConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                          real k_bolt, real mub, real damping);
   void resetConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                       real k_bolt, real mub, real damping);
   bool initiated() const { return dataInitiated && constantsInitiated; }
   const GpuTensor<real, 3>& getField() const { return field; }
   void randomize(const GpuTensor<real, 2>& mmom);
};
