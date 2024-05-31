#include <cuda.h>
#include <curand.h>

#include "c_headers.hpp"
#include "cudaParallelizationHelper.hpp"
#include "cudaThermfield.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

// The neighbour list setup helper
class CudaThermfield::SetupSigmaFactor : public CudaParallelizationHelper::Site {
private:
   real* sigma_factor;
   real dp;

public:
   SetupSigmaFactor(CudaTensor<real, 1>& p1, real p2) {
      sigma_factor = p1.data();
      dp = p2;
   }

   __device__ void each(unsigned int site) {
      sigma_factor[site] = sqrt(dp * sigma_factor[site]);
   }
};

class CudaThermfield::SetupField : public CudaParallelizationHelper::AtomSite {
private:
   real* field;
   const real* sigma_factor;
   const real* mmom;

public:
   SetupField(CudaTensor<real, 3>&  p1, const CudaTensor<real, 1>& p2, const CudaTensor<real, 2>& p3) {
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

CudaThermfield::CudaThermfield()
    : stopwatch(GlobalStopwatchPool::get("Cuda thermfield")),
      parallel(CudaParallelizationHelper::def) {
   constantsInitiated = false;
   dataInitiated = false;
}

CudaThermfield::~CudaThermfield() {
   if(dataInitiated) {
      curandDestroyGenerator(gen);
   }
}

bool CudaThermfield::initiate(std::size_t N, std::size_t M, curandRngType_t rngType,
                              unsigned long long seed) {
   if(dataInitiated) {
      std::fprintf(stderr, "Warning: attempt to initiate already initiated CudaThermfield\n");
      return true;
   }

   stopwatch.skip();
   field.Allocate(3, N, M);
   sigmaFactor.Allocate(N);
   if(!field.empty() && !sigmaFactor.empty()) {
      if(curandCreateGenerator(&gen, rngType) == CURAND_STATUS_SUCCESS) {
         if(seed == 0ULL) {
            seed = time(nullptr);
         }
         curandSetPseudoRandomGeneratorSeed(gen, seed);
         curandSetStream(gen, parallel.getWorkStream());
         dataInitiated = true;
      } else {
         field.Free();
         sigmaFactor.Free();
      }
   }
   stopwatch.add("initiate");
   return dataInitiated;
}

bool CudaThermfield::initiateConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
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
   parallel.cudaSiteCall(SetupSigmaFactor(sigmaFactor, dp));
   stopwatch.add("initiate constants");

   constantsInitiated = true;
   return true;
}

void CudaThermfield::randomize(const CudaTensor<real, 2>& mmom) {
   // Initiated?
   if(!initiated()) {
      return;
   }

   // Timing
   stopwatch.skip();

// Generate random vector
#ifdef SINGLE_PREC
   curandGenerateNormal(gen, field.data(), field.size(), 0.0, 1.0);
#else
   curandGenerateNormalDouble(gen, field.data(), field.size(), 0.0, 1.0);
#endif
   stopwatch.add("RNG");

   // Expand thermal field
   parallel.cudaAtomSiteCall(SetupField(field, sigmaFactor, mmom));
   stopwatch.add("loop");
}

