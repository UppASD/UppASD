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
   // Projected device bytes for a given system size; must mirror initiate().
   static std::size_t estimateBytes(std::size_t N, std::size_t M);
   const GpuTensor<real, 3>& getField() const { return field; }
   void randomize(const GpuTensor<real, 2>& mmom);
   // Strategy 1 (CGP-06B): identical generator sequence and identical
   // randomValues draw to randomize(mmom) above -- state advances exactly as
   // today so invariants 1-5 (CGP-06A Part B) hold unchanged -- but the
   // field write is scoped to the compact active list instead of every
   // (site,ensemble) pair. oneBasedAtoms/activeCount use the same one-based
   // compact-list contract as GpuParallelizationHelper::gpuActiveAtomCall.
   void randomize(const GpuTensor<real, 2>& mmom, const int* oneBasedAtoms,
                  unsigned int activeCount);
};
