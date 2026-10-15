#include "c_headers.hpp"
#include "gpuCommon.hpp"
#include "gpuDepondtIntegrator.hpp"
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
GpuDepondtIntegrator_fused::GpuDepondtIntegrator_fused(const SimulationParameters SimParam, deviceLattice& p_gpuLattice)
    : stopwatch(GlobalStopwatchPool::get("Cuda Depondt integrator")),
      gpuLattice(p_gpuLattice),
{
   gamma = SimParam.gamma;
   k_bolt = SimParam.k_bolt;
   mub = SimParam.mub;
   damping = SimParam.damping;
   dp_factor = (2.0 * damping * k_bolt) / (gamma * mub * (1 + damping * damping));

   // Integrator parameters
   stt = SimParam.stt;
   timestep = SimParam.delta_t;

   applyTorque = ((stt =='Y') ? true : false);

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

   InitializeDeviceRngStates<<<taskNM, 1>>>(d_state.data(), seed, taskNM);
   InitializeThermalSigmaFactor<<<blocks, threads>>>(sigmaFactor, gpuLattice.temperature, N, dp_factor);

   if(GPU_DEVICE_SYNCHRONIZE() != GPU_SUCCESS) {
      release();
      return false;
   }
}

// Destructor
GpuDepondtIntegrator_fused::~GpuDepondtIntegrator_fused() {
   release();
}




void GpuDepondtIntegrator::resetConstants(const Tensor<real, 1> temperature, real phaseTimestep, real phaseDamping) {
   // Update per-phase constants before resetting thermfield state.
   timestep = phaseTimestep;
   damping = phaseDamping;
   dp_factor = (2.0 * damping * k_bolt) / (gamma * mub * (1 + damping * damping));
   thermfield.resetConstants(temperature, timestep, gamma, k_bolt, mub, damping); 
}

// Releaser
void GpuDepondtIntegrator::release() {
   thermfield.Free();
   sigmaFactor.Free();
   d_state.Free();
   dataInitiated = false;
}

// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field
// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
void GpuDepondtIntegrator::evolveFirst(deviceLattice& gpuLattice) {
   // Timing
   stopwatch.skip();
   threads = {256, 1, 1};
   blocks = {(N + threads.x - 1)/threads.x, M, 1};

   DepondtEvolveFirst<<<blocks, threads>>>(thermfield, sigmaFactor, d_state, gpuLattice.beff, gpuLattice.emom, gpuLattice.emom2, gpuLattice.mmom, 
                                 gpuLattice.emomM, gpuLattice.btorque, N, M, applyTorque);


///////////////////////////
   // copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
   parallel.gpuElementCall(GpuCommon::ScalarMult(gpuLattice.emomM, mrod, gpuLattice.mmom));
   gpuLattice.emom2.swap(gpuLattice.emom);  // Previous emom will not be needed
   gpuLattice.emom.swap(mrod);   // Previous mrod will not be needed
   gpuLattice.b2eff.swap(bdup);  // Previous bdup will not be needed
   stopwatch.add("copy");
}

// Second step of Depond solver, calculates the corrected effective field from
// the predicted effective fields. Rotates the moments in the corrected field
void GpuDepondtIntegrator::evolveSecond(deviceLattice& gpuLattice) {
   // Timing
   stopwatch.skip();

   // Construct local field
   parallel.gpuElementCall(GpuCommon::Add(blocal, gpuLattice.beff, thermfield.getField()));
   stopwatch.add("localfield");

   // Construct effective field (including damping term)
   buildbeff(gpuLattice.emom, gpuLattice.btorque);
   stopwatch.add("buildbeff");

   // Corrected field
   parallel.gpuElementCall(GpuCommon::Avg(bdup, gpuLattice.b2eff));
   gpuLattice.emom.swap(gpuLattice.emom2);  // Vaild as emom2 wont be used after its coming swap
   stopwatch.add("corrfield");

   // Final rotation
   rotate(gpuLattice.emom, timestep);
   stopwatch.add("rotate");

   // Swap
   gpuLattice.emom2.swap(mrod);  // Vaild as mrod wont be needed after this
   stopwatch.add("copy");
}

void GpuDepondtIntegrator::rotate(const GpuTensor<real, 3>& emom, real delta_t) {
   parallel.gpuAtomCall(Rotate(mrod, emom, bdup, timestep, gamma, damping));
}

// Constructs the effective field (including damping term)
void GpuDepondtIntegrator::buildbeff(const GpuTensor<real, 3>& emom,
                                      const GpuTensor<real, 3>& btorque) {
   parallel.gpuAtomCall(BuildEffectiveField(bdup, blocal, emom, damping));

   // TODO untested
   if(stt != 'N') {
      parallel.gpuElementCall(GpuCommon::AddTo(bdup, btorque));
   }
}

