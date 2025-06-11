#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>

#include "c_headers.hpp"
#include "gpuParallelizationHelper.hpp"
#include "hipThermfield.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

// The neighbour list setup helper
class HipThermfield::SetupSigmaFactor : public GpuParallelizationHelper::Site {
private:
   real* sigma_factor;
   real dp;

public:
   SetupSigmaFactor(GpuTensor<real, 1>& p1, real p2) {
      sigma_factor = p1.data();
      dp = p2;
   }

   __device__ void each(unsigned int site) {
      sigma_factor[site] = sqrt(dp * sigma_factor[site]);
   }
};

class HipThermfield::SetupField : public GpuParallelizationHelper::AtomSite {
private:
   real* field;
   const real* sigma_factor;
   const real* mmom;

public:
   SetupField(GpuTensor<real, 3>&  p1, const GpuTensor<real, 1>& p2, const GpuTensor<real, 2>& p3) {
      field = p1.data();
      sigma_factor = p2.data();
      mmom = p3.data();
   }

   __device__ void each(unsigned int atom, unsigned int site) {
      real sigma = sigma_factor[site] * rsqrt(mmom[atom]);
      field[atom * 3 + 0] *= sigma;
      field[atom * 3 + 1] *= sigma;
      field[atom * 3 + 2] *= sigma;
   }
};

////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

HipThermfield::HipThermfield()
    : stopwatch(GlobalStopwatchPool::get("Hip thermfield")),
      parallel(ParallelizationHelperInstance) {
   constantsInitiated = false;
   dataInitiated = false;
}

HipThermfield::~HipThermfield() {
   if(dataInitiated) {
      hiprandDestroyGenerator(gen);
   }
}

bool HipThermfield::initiate(std::size_t N, std::size_t M, hiprandRngType_t rngType,
                              unsigned long long seed) {
   if(dataInitiated) {
      std::fprintf(stderr, "Warning: attempt to initiate already initiated HipThermfield\n");
      return true;
   }

   stopwatch.skip();
   field.Allocate(3, N, M);
   sigmaFactor.Allocate(N);
   if(!field.empty() && !sigmaFactor.empty()) {
      if(hiprandCreateGenerator(&gen, rngType) == HIPRAND_STATUS_SUCCESS) {
         if(seed == 0ULL) {
            seed = time(nullptr);
         }
         hiprandSetPseudoRandomGeneratorSeed(gen, seed);
         hiprandSetStream(gen, parallel.getWorkStream());
         dataInitiated = true;
      } else {
         field.Free();
         sigmaFactor.Free();
      }
   }
   stopwatch.add("initiate");
   return dataInitiated;
}

bool HipThermfield::initiateConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                                       real k_bolt, real mub, real damping) {
   // Timing
   stopwatch.skip();

   // Initiated?
   if(!dataInitiated) {
      return false;
   }

   // Damping parameter
   real dp = (2.0 * damping * k_bolt) / (timestep * gamma * mub * (1 + damping * damping));

   // Set up sigmaFactor
   sigmaFactor.copy_sync(temperature);

   // sF = sqrt(dp*sF) ( = sqrt(dp*temp))
   parallel.gpuSiteCall(SetupSigmaFactor(sigmaFactor, dp));
   stopwatch.add("initiate constants");

   constantsInitiated = true;
   return true;
}

void HipThermfield::resetConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                                       real k_bolt, real mub, real damping) {
   // Set up sigmaFactor
   sigmaFactor.copy_sync(temperature);
   real dp = (2.0 * damping * k_bolt) / (timestep * gamma * mub * (1 + damping * damping));
   // sF = sqrt(dp*sF) ( = sqrt(dp*temp))
   parallel.gpuSiteCall(SetupSigmaFactor(sigmaFactor, dp));
   stopwatch.add("initiate constants");
}

void  HipThermfield::randomize(const GpuTensor<real, 2>& mmom) {
   // Initiated?
   if(!initiated()) {
      return;
   }

   // Timing
   stopwatch.skip();

// Generate random vector
#ifdef SINGLE_PREC
   hiprandGenerateNormal(gen, field.data(), field.size(), 0.0, 1.0);
#else
   hiprandGenerateNormalDouble(gen, field.data(), field.size(), 0.0, 1.0);
#endif
   stopwatch.add("RNG");

   // Expand thermal field
   parallel.gpuAtomSiteCall(SetupField(field, sigmaFactor, mmom));
   stopwatch.add("loop");
}

