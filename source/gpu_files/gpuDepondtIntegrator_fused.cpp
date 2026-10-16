#include "c_headers.hpp"
#include "gpuCommon.hpp"
#include "gpuDepondtIntegrator_fused.hpp"
#include "tensor.hpp"
#include "printDebug.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "gpuStructures.hpp"
#include "depondt_kernels.hpp"

#include "gpu_wrappers.h"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <cuda.h>
#include <curand.h>
#endif


// ConstructoruDepondtIntegrator
GpuDepondtIntegrator_fused::GpuDepondtIntegrator_fused(const SimulationParameters SimParam, deviceLattice& p_gpuLattice, const Tensor<real, 1> temperature)
    : stopwatch(GlobalStopwatchPool::get("Cuda Depondt integrator")),
      gpuLattice(p_gpuLattice)
{
   gamma = SimParam.gamma;
   k_bolt = SimParam.k_bolt;
   mub = SimParam.mub;
   damping = SimParam.damping;
   timestep = SimParam.delta_t;

   //dp_factor = (2.0 * damping * k_bolt) / (gamma * mub * (1 + damping * damping));
    dp_factor = (2.0 * damping * k_bolt) / (timestep * gamma * mub *(1.0 + damping * damping));


   // Integrator parameters
   stt = SimParam.stt;

   applyTorque = ((stt !='N') ? true : false);

   N = SimParam.N;
   M = SimParam.M;

   tasks3NM = 3*N*M;
   tasksNM = N*M;


   threads = {256, 1, 1};
   blocks = {(N + threads.x - 1) / threads.x, 1, 1};

   thermfield.Allocate(3, N, M);
   sigmaFactor.Allocate(SimParam.N);
   d_state.Allocate(N * M);
   dataInitiated = true;

   unsigned long long seed;
   if(SimParam.randomSeed == 0ULL) {
         srand(time(NULL));
         seed = (unsigned long long)rand();
      }
   else{
      seed = (unsigned long long) SimParam.randomSeed;
   }
   const dim3 rngBlocks((tasksNM + threads.x - 1) / threads.x, 1, 1);
   sigmaFactor.copy_sync(temperature);


   InitializeDeviceRngStates<<<rngBlocks, threads>>>(d_state.data(), seed, taskNM);
   InitializeThermalSigmaFactor<<<blocks, threads>>>(sigmaFactor,  N, dp_factor);

   const GPU_ERROR_T err = GPU_DEVICE_SYNCHRONIZE();
    if (err != GPU_SUCCESS) {
        std::fprintf(
            stderr,
            "Failed to reset Depondt thermal constants: %s\n",
            GPU_GET_ERROR_STRING(err));

         constantsInitiated = true;
    }

    constantsInitiated = true;

}

// Destructor
GpuDepondtIntegrator_fused::~GpuDepondtIntegrator_fused() {
   release();
}
void GpuDepondtIntegrator_fused::resetConstants(const Tensor<real, 1> temperature, real phaseTimestep, real phaseDamping){
    timestep = phaseTimestep;
    damping = phaseDamping;

    dp_factor = (2.0 * damping * k_bolt) / (timestep * gamma * mub *(1.0 + damping * damping));

    const dim3 sigmaThreads(256, 1, 1);
    const dim3 sigmaBlocks((N + sigmaThreads.x - 1) / sigmaThreads.x,  1, 1);

    sigmaFactor.copy_sync(temperature);
    InitializeThermalSigmaFactor<<<sigmaBlocks, sigmaThreads>>>(sigmaFactor, N, dp_factor);

    const GPU_ERROR_T err = GPU_DEVICE_SYNCHRONIZE();
    if (err != GPU_SUCCESS) {
        std::fprintf(
            stderr,
            "Failed to reset Depondt thermal constants: %s\n",
            GPU_GET_ERROR_STRING(err));
         constantsInitiated = false;

    }

    constantsInitiated = true;
}



// Releaser
void GpuDepondtIntegrator_fused::release() {
   thermfield.Free();
   sigmaFactor.Free();
   d_state.Free();
   dataInitiated = false;
}

// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field
// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
void GpuDepondtIntegrator_fused::evolveFirst(deviceLattice& gpuLattice) {
   threads = {256, 1, 1};
   blocks = {(N + threads.x - 1)/threads.x, M, 1};

   DepondtEvolveFirst<<<blocks, threads>>>(thermfield, sigmaFactor, d_state, gpuLattice.beff, gpuLattice.b2eff, gpuLattice.emom, gpuLattice.emom2, gpuLattice.mmom, 
                                 gpuLattice.emomM, gpuLattice.btorque, N, M, timestep, gamma, damping, applyTorque);

}

// Second step of Depond solver, calculates the corrected effective field from
// the predicted effective fields. Rotates the moments in the corrected field
void GpuDepondtIntegrator_fused::evolveSecond(deviceLattice& gpuLattice) {
   threads = {256, 1, 1};
   blocks = {(N + threads.x - 1)/threads.x, M, 1};
   DepondtEvolveSecond<<<blocks, threads>>>(thermfield, gpuLattice.beff, gpuLattice.b2eff, gpuLattice.emom,gpuLattice.emom2, gpuLattice.btorque,
                                            N, M, timestep,gamma,damping,applyTorque);

}

