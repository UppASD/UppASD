#pragma once

#include "c_headers.hpp"

#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "gpuStructures.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#include "hipThermfield.hpp"
using GpuThermfield = HipThermfield;
#elif defined(CUDA_V)
#include <curand.h>
#include "cudaThermfield.hpp"
using GpuThermfield = CudaThermfield;
#endif
using ParallelizationHelper = GpuParallelizationHelper;

class GpuDepondtIntegrator {
private:
   // System parameters
   real gamma;
   real k_bolt;
   real mub;
   real damping;
   real dp_factor;

   // Integrator parameters
   char stt;  // Spin transfer torque mode
   real timestep;

   // Class local matrices
   GpuTensor<real, 3> mrod;
   GpuTensor<real, 3> blocal;
   GpuTensor<real, 3> bdup;

   // Thermfield
   GpuThermfield thermfield;

   // Timer
   StopwatchDeviceSync stopwatch;

   // Parallelization helper
   ParallelizationHelper& parallel;

   // Algorithm
   void rotate(const GpuTensor<real, 3>& emom, real delta_t);
   void buildbeff(const GpuTensor<real, 3>& emom, const GpuTensor<real, 3>& btorque);

public:
   // Parallelization helpers
   class Rotate;
   class BuildEffectiveField;

   // Constructor
   GpuDepondtIntegrator();

   // Destructor
   ~GpuDepondtIntegrator();

   // Initiator
   bool initiate(const SimulationParameters SimParam);

   // Set up constants
   bool initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature);
   void resetConstants(const Tensor<real, 1> temperature);

   // Releaser
   void release();

   // Algorithm
   void evolveFirst(deviceLattice& gpuLattice);

   void evolveSecond(deviceLattice& gpuLattice);
};

