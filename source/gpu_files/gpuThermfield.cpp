#include "gpuThermfield.hpp"
#include "base.hpp"
#include "stopwatchPool.hpp"

class GpuThermfield::SetupSigmaFactor : public GpuParallelizationHelper::Site {
   real* sigmaFactor; real dp;
public:
   SetupSigmaFactor(GpuTensor<real, 1>& sigma, real factor) : sigmaFactor(sigma.data()), dp(factor) {}
   __device__ void each(unsigned int site) { sigmaFactor[site] = sqrt(dp * sigmaFactor[site]); }
};

class GpuThermfield::SetupField : public GpuParallelizationHelper::AtomSite {
   real* field; const real* randomValues; const real* sigmaFactor; const real* mmom;
public:
   SetupField(GpuTensor<real, 3>& fieldIn, const GpuVector<real>& randomIn,
              const GpuTensor<real, 1>& sigmaIn, const GpuTensor<real, 2>& mmomIn)
      : field(fieldIn.data()), randomValues(randomIn.data()), sigmaFactor(sigmaIn.data()), mmom(mmomIn.data()) {}
   __device__ void each(unsigned int atom, unsigned int site) {
      const real sigma = sigmaFactor[site] * rsqrt(mmom[atom]);
      field[atom * 3] = randomValues[atom * 3] * sigma;
      field[atom * 3 + 1] = randomValues[atom * 3 + 1] * sigma;
      field[atom * 3 + 2] = randomValues[atom * 3 + 2] * sigma;
   }
};

GpuThermfield::GpuThermfield()
   : stopwatch(GlobalStopwatchPool::get("GPU thermfield"), ParallelizationHelperInstance.getWorkStream()), parallel(ParallelizationHelperInstance) {}

GpuThermfield::~GpuThermfield() {
   if(dataInitiated) GPU_RAND_DESTROY_GEN(gen);
   // Tensor storage has explicit ownership; release all three buffers in the
   // same destructor that tears down the RNG so a completed SD phase leaves
   // no hidden device allocation behind.
   field.Free();
   sigmaFactor.Free();
   randomValues.Free();
}

// Mirror of the allocations in initiate() below. Keep the two in sync.
std::size_t GpuThermfield::estimateBytes(std::size_t N, std::size_t M) {
   const std::size_t field_count = 3 * N * M;                       // field(3,N,M)
   const std::size_t random_count = field_count + (field_count & 1);// randomValues, padded even
   const std::size_t sigma_count = N;                               // sigmaFactor(N)
   return (field_count + random_count + sigma_count) * sizeof(real);
}

bool GpuThermfield::initiate(std::size_t N, std::size_t M, GPU_RAND_RNGTYPE_T rngType, unsigned long long seed) {
   if(dataInitiated) return true;
   stopwatch.skip();
   field.Allocate(3, N, M);
   randomValues.Allocate(static_cast<long int>(field.size() + (field.size() & 1)));
   sigmaFactor.Allocate(N);
   if(field.empty() || randomValues.empty() || sigmaFactor.empty()) return false;
   if(GPU_RAND_CREATE_GEN(&gen, rngType) != GPU_RAND_STATUS_SUCCESS) return false;
   if(seed == 0ULL) seed = time(nullptr);
   try {
      ASSERT_CURAND(GPU_RAND_SET_PSEUDO_RANDOM_GENERATOR_SEED(gen, seed));
      ASSERT_CURAND(GPU_RAND_SET_STREAM(gen, parallel.getWorkStream()));
   } catch(...) {
      GPU_RAND_DESTROY_GEN(gen);
      field.Free(); sigmaFactor.Free(); randomValues.Free();
      throw;
   }
   dataInitiated = true;
   stopwatch.add("initiate");
   return true;
}

bool GpuThermfield::initiateConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                                      real k_bolt, real mub, real damping) {
   if(!dataInitiated) return false;
   resetConstants(temperature, timestep, gamma, k_bolt, mub, damping);
   constantsInitiated = true;
   return true;
}

void GpuThermfield::resetConstants(const Tensor<real, 1>& temperature, real timestep, real gamma,
                                   real k_bolt, real mub, real damping) {
   stopwatch.skip();
   sigmaFactor.copy_sync(temperature);
   const real dp = (2.0 * damping * k_bolt) / (timestep * gamma * mub * (1 + damping * damping));
   parallel.gpuSiteCall(SetupSigmaFactor(sigmaFactor, dp));
   stopwatch.add("initiate constants");
}

void GpuThermfield::randomize(const GpuTensor<real, 2>& mmom) {
   if(!initiated()) return;
   stopwatch.skip();
   ASSERT_CURAND(GPU_RAND_GENERATE_NORMAL_REAL(gen, randomValues.data(), randomValues.size(),
                                                static_cast<real>(0.0), static_cast<real>(1.0)));
   stopwatch.add("RNG");
   parallel.gpuAtomSiteCall(SetupField(field, randomValues, sigmaFactor, mmom));
   stopwatch.add("loop");
}

void GpuThermfield::randomize(const GpuTensor<real, 2>& mmom, const int* oneBasedAtoms,
                              unsigned int activeCount) {
   if(!initiated()) return;
   stopwatch.skip();
   // Unconditional, full-N,M draw -- identical call to randomize(mmom) above.
   // The generator's sequence must advance exactly as it does in the
   // feature-off/all-fine path (CGP-06A invariants 1-4); only the
   // downstream write is scoped below.
   ASSERT_CURAND(GPU_RAND_GENERATE_NORMAL_REAL(gen, randomValues.data(), randomValues.size(),
                                                static_cast<real>(0.0), static_cast<real>(1.0)));
   stopwatch.add("RNG");
   // Active-list write only: every downstream reader of `field` for adaptive
   // dynamics (GpuDepondtIntegrator's active-list EvolveFirst) is itself
   // restricted to the same compact list, so entries left unwritten here are
   // never read. Empty list (all-coarse step) issues no write kernel at all.
   parallel.gpuActiveAtomSiteCall(SetupField(field, randomValues, sigmaFactor, mmom),
                                  oneBasedAtoms, activeCount);
   stopwatch.add("loop");
}
