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
#include "gpuParallelizationHelper.hpp"

#include "gpu_wrappers.h"
#include <limits>
#include <stdexcept>
using ParallelizationHelper = GpuParallelizationHelper;

// Versine, 1 - cos(angle), evaluated as 2*sin(angle/2)^2.  The direct form
// cancels catastrophically for the small precession angles the Depondt solver
// normally takes; the half-angle form keeps full relative accuracy there.
static __device__ inline real versine(real angle) {
   const real s = sin(angle * real(0.5));
   return real(2.0) * s * s;
}

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

class GpuDepondtIntegrator::BuildEffectiveField : public ParallelizationHelper::Atom {
private:
   real* bdup;
   const real* blocal;
   const real* emom;
   real damping;

public:
   BuildEffectiveField(GpuTensor<real, 3>& p1, const GpuTensor<real, 3>& p2, const GpuTensor<real, 3>& p3, real p4) {
      bdup = p1.data();
      blocal = p2.data();
      emom = p3.data();
      damping = p4;
   }

   __device__ void each(unsigned int atom) {
      real* my_bdup = &bdup[atom * 3];
      const real* my_bloc = &blocal[atom * 3];
      const real* my_emom = &emom[atom * 3];
      my_bdup[0] = my_bloc[0] + damping * (my_emom[1] * my_bloc[2] - my_emom[2] * my_bloc[1]);
      my_bdup[1] = my_bloc[1] + damping * (my_emom[2] * my_bloc[0] - my_emom[0] * my_bloc[2]);
      my_bdup[2] = my_bloc[2] + damping * (my_emom[0] * my_bloc[1] - my_emom[1] * my_bloc[0]);
      //if (atom == 0) printf("%.6lf \n", my_emom[0]); 
   }
};

class GpuDepondtIntegrator::Rotate : public ParallelizationHelper::Atom {
private:
   real* mrod;
   const real* emom;
   const real* bdup;
   real timestep;
   real gamma;
   real damping;

public:
   Rotate(GpuTensor<real, 3>& p1, const GpuTensor<real, 3>& p2, const GpuTensor<real, 3>& p3, real p4, real p5, real p6) {
      mrod = p1.data();
      emom = p2.data();
      bdup = p3.data();
      timestep = p4;
      gamma = p5;
      damping = p6;
   }

   __device__ void each(unsigned int atom) {
      real* my_mrod = &mrod[atom * 3];
      const real* my_emom = &emom[atom * 3];
      const real* my_bdup = &bdup[atom * 3];

      // Get effective field components and size
      real x = my_bdup[0];
      real y = my_bdup[1];
      real z = my_bdup[2];
      real norm = sqrt(x * x + y * y + z * z);

      if(norm == (real)0.0) {
         my_mrod[0] = my_emom[0];
         my_mrod[1] = my_emom[1];
         my_mrod[2] = my_emom[2];
         return;
      }

      // Normalize components
      x /= norm;
      y /= norm;
      z /= norm;

      // Get precession angle
      real angle = norm * timestep * gamma;

      angle *= (real)1.0 / ((real)1.0 + damping * damping);

      // Calculate sin and cos
      real cosv, sinv;
      sincos(angle, &sinv, &cosv);
      real u = versine(angle);

      // Calculate matrix
      real M[3][3];
      M[0][0] = x * x * u + cosv;
      M[0][1] = x * y * u - z * sinv;
      M[0][2] = x * z * u + y * sinv;
      M[1][0] = x * y * u + z * sinv;
      M[1][1] = y * y * u + cosv;
      M[1][2] = y * z * u - x * sinv;
      M[2][0] = x * z * u - y * sinv;
      M[2][1] = z * y * u + x * sinv;
      M[2][2] = z * z * u + cosv;

      // Rotate
      real mx = my_emom[0], my = my_emom[1], mz = my_emom[2];
      my_mrod[0] = mx * M[0][0] + my * M[0][1] + mz * M[0][2];
      my_mrod[1] = mx * M[1][0] + my * M[1][1] + mz * M[1][2];
      my_mrod[2] = mx * M[2][0] + my * M[2][1] + mz * M[2][2];
      // Alternative
      //		real mx = my_emom[0], my = my_emom[1], mz = my_emom[2];
      //		my_mrod[0] = mx*x*x*u+  cosv + my*y*x*u-z*sinv + mz*z*x*u+y*sinv;
      //		my_mrod[1] = mx*x*y*u+z*sinv + my*y*y*u+  cosv + mz*z*y*u-x*sinv;
      //		my_mrod[2] = mx*x*z*u-y*sinv + my*y*z*u+x*sinv + mz*z*z*u+  cosv;
   }
};

// Predictor-only fused pass.  blocal must remain available after this launch:
// evolveSecond uses the same thermal field when constructing its corrected
// effective field.  bdup and mrod are deliberately written to the same
// temporary buffers as the former individual passes, so their swaps below
// retain their original ownership semantics.
class GpuDepondtIntegrator::EvolveFirst : public ParallelizationHelper::Atom {
private:
   real* mrod;
   real* blocal;
   real* bdup;
   real* emomM;
   const real* emom;
   const real* beff;
   const real* thermalField;
   const real* btorque;
   const real* mmom;
   real timestep;
   real gamma;
   real damping;
   bool addStt;

public:
   EvolveFirst(GpuTensor<real, 3>& pMrod,
               GpuTensor<real, 3>& pLocal,
               GpuTensor<real, 3>& pBdup,
               GpuTensor<real, 3>& pEmomM,
               const GpuTensor<real, 3>& pEmom,
               const GpuTensor<real, 3>& pBeff,
               const GpuTensor<real, 3>& pThermalField,
               const GpuTensor<real, 3>& pTorque,
               const GpuTensor<real, 2>& pMom,
               real pTimestep, real pGamma, real pDamping, bool pAddStt)
       : mrod(pMrod.data()), blocal(pLocal.data()), bdup(pBdup.data()),
         emomM(pEmomM.data()), emom(pEmom.data()), beff(pBeff.data()),
         thermalField(pThermalField.data()), btorque(pTorque.data()),
         mmom(pMom.data()), timestep(pTimestep), gamma(pGamma),
         damping(pDamping), addStt(pAddStt) {}

   __device__ void each(unsigned int atom) {
      const unsigned int element = atom * 3;
      const real bx = beff[element] + thermalField[element];
      const real by = beff[element + 1] + thermalField[element + 1];
      const real bz = beff[element + 2] + thermalField[element + 2];
      blocal[element] = bx;
      blocal[element + 1] = by;
      blocal[element + 2] = bz;

      const real mx = emom[element];
      const real my = emom[element + 1];
      const real mz = emom[element + 2];
      real x = bx + damping * (my * bz - mz * by);
      real y = by + damping * (mz * bx - mx * bz);
      real z = bz + damping * (mx * by - my * bx);
      if(addStt) {
         x += btorque[element];
         y += btorque[element + 1];
         z += btorque[element + 2];
      }
      bdup[element] = x;
      bdup[element + 1] = y;
      bdup[element + 2] = z;

      const real norm = sqrt(x * x + y * y + z * z);
      if(norm == real(0.0)) {
         mrod[element] = mx;
         mrod[element + 1] = my;
         mrod[element + 2] = mz;
      } else {
         x /= norm;
         y /= norm;
         z /= norm;
         real angle = norm * timestep * gamma;
         angle *= real(1.0) / (real(1.0) + damping * damping);
         real cosv, sinv;
         sincos(angle, &sinv, &cosv);
         const real u = versine(angle);

         const real r0 = mx * (x * x * u + cosv) +
                         my * (x * y * u - z * sinv) +
                         mz * (x * z * u + y * sinv);
         const real r1 = mx * (x * y * u + z * sinv) +
                         my * (y * y * u + cosv) +
                         mz * (y * z * u - x * sinv);
         const real r2 = mx * (x * z * u - y * sinv) +
                         my * (z * y * u + x * sinv) +
                         mz * (z * z * u + cosv);
         mrod[element] = r0;
         mrod[element + 1] = r1;
         mrod[element + 2] = r2;
      }

      const real moment = mmom[atom];
      emomM[element] = mrod[element] * moment;
      emomM[element + 1] = mrod[element + 1] * moment;
      emomM[element + 2] = mrod[element + 2] * moment;
   }
};

// The full-range solver expresses the predictor-to-corrector transition with
// tensor swaps.  Swapping a whole tensor is not valid for a subset because it
// would exchange inactive atom state too.  These copy-only operations preserve
// the public Depondt buffer contract for the selected physical atoms while the
// actual predictor/corrector arithmetic stays in EvolveFirst,
// BuildEffectiveField, Rotate, and GpuCommon.
class GpuDepondtIntegrator::CommitActivePredictor : public ParallelizationHelper::Atom {
private:
   real* emom;
   real* emom2;
   real* b2eff;
   const real* mrod;
   const real* bdup;
public:
   CommitActivePredictor(GpuTensor<real, 3>& pEmom, GpuTensor<real, 3>& pEmom2,
                         GpuTensor<real, 3>& pB2eff, const GpuTensor<real, 3>& pMrod,
                         const GpuTensor<real, 3>& pBdup)
      : emom(pEmom.data()), emom2(pEmom2.data()), b2eff(pB2eff.data()),
        mrod(pMrod.data()), bdup(pBdup.data()) {}
   __device__ void each(unsigned int atom) {
      const unsigned int element = 3 * atom;
      emom2[element] = emom[element];
      emom2[element + 1] = emom[element + 1];
      emom2[element + 2] = emom[element + 2];
      emom[element] = mrod[element];
      emom[element + 1] = mrod[element + 1];
      emom[element + 2] = mrod[element + 2];
      b2eff[element] = bdup[element];
      b2eff[element + 1] = bdup[element + 1];
      b2eff[element + 2] = bdup[element + 2];
   }
};

class GpuDepondtIntegrator::RestoreActiveInitial : public ParallelizationHelper::Atom {
private:
   real* emom;
   const real* emom2;
public:
   RestoreActiveInitial(GpuTensor<real, 3>& pEmom, const GpuTensor<real, 3>& pEmom2)
      : emom(pEmom.data()), emom2(pEmom2.data()) {}
   __device__ void each(unsigned int atom) {
      const unsigned int element = 3 * atom;
      emom[element] = emom2[element];
      emom[element + 1] = emom2[element + 1];
      emom[element + 2] = emom2[element + 2];
   }
};

class GpuDepondtIntegrator::CommitActiveCorrector : public ParallelizationHelper::Atom {
private:
   real* emom2;
   const real* mrod;
public:
   CommitActiveCorrector(GpuTensor<real, 3>& pEmom2, const GpuTensor<real, 3>& pMrod)
      : emom2(pEmom2.data()), mrod(pMrod.data()) {}
   __device__ void each(unsigned int atom) {
      const unsigned int element = 3 * atom;
      emom2[element] = mrod[element];
      emom2[element + 1] = mrod[element + 1];
      emom2[element + 2] = mrod[element + 2];
   }
};

////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

// ConstructoruDepondtIntegrator
GpuDepondtIntegrator::GpuDepondtIntegrator()
    : stopwatch(GlobalStopwatchPool::get("Cuda Depondt integrator"), ParallelizationHelperInstance.getWorkStream()),
      parallel(ParallelizationHelperInstance) {
}

// Destructor
GpuDepondtIntegrator::~GpuDepondtIntegrator() {
   release();
}

// Mirror of the device allocations in initiate() (mrod/blocal/bdup, each (3,N,M))
// plus the thermfield this integrator owns. Keep in sync with initiate().
std::size_t GpuDepondtIntegrator::estimateBytes(const SimulationParameters& SimParam) {
   const std::size_t nm3 = 3 * SimParam.N * SimParam.M;
   const std::size_t local = 3 * nm3 * sizeof(real);  // mrod + blocal + bdup
   return local + GpuThermfield::estimateBytes(SimParam.N, SimParam.M);
}

// Initiator
bool GpuDepondtIntegrator::initiate(const SimulationParameters SimParam) {
   // Assert that we're not already initialized
   release();

   // Param
   stt = SimParam.stt;
   timestep = SimParam.delta_t;

   // Initiate thermfield
   if(!thermfield.initiate(SimParam.N, SimParam.M, SimParam.rngType, SimParam.randomSeed)) {
      release();
      return false;
   }

   // Allocate device matrices
   mrod.Allocate(static_cast <long int>(3),static_cast <long int>( SimParam.N), static_cast <long int>(SimParam.M));
   blocal.Allocate(static_cast <long int>(3), static_cast <long int>(SimParam.N), static_cast <long int>(SimParam.M));
   bdup.Allocate(static_cast <long int>(3), static_cast <long int>(SimParam.N), static_cast <long int>(SimParam.M));

   // All initialized?
   if(GPU_DEVICE_SYNCHRONIZE() != GPU_SUCCESS) {
      release();
      return false;
   }

   return true;
}

bool GpuDepondtIntegrator::initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature) {
   // Set parameters
   gamma = SimParam.gamma;
   k_bolt = SimParam.k_bolt;
   mub = SimParam.mub;
   damping = SimParam.damping;
   dp_factor = (2.0 * damping * k_bolt) / (gamma * mub * (1 + damping * damping));

   // Initiate thermfield constants
   if(!thermfield.initiateConstants(temperature, timestep, gamma, k_bolt, mub, damping)) {
      return false; //TODO
   }

   return true;
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
   mrod.Free();
   blocal.Free();
   bdup.Free();
}

// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field
// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
void GpuDepondtIntegrator::evolveFirst(deviceLattice& gpuLattice) {
   // Timing
   stopwatch.skip();

   //_dpr;

   // Randomize stochastic field
   thermfield.randomize(gpuLattice.mmom);
   stopwatch.add("thermfield");

   parallel.gpuAtomCall(EvolveFirst(mrod, blocal, bdup, gpuLattice.emomM,
                                    gpuLattice.emom, gpuLattice.beff,
                                    thermfield.getField(), gpuLattice.btorque,
                                    gpuLattice.mmom, timestep, gamma, damping,
                                    stt != 'N'));
   stopwatch.add("predictor");

   gpuLattice.emom2.swap(gpuLattice.emom);
   gpuLattice.emom.swap(mrod);
   gpuLattice.b2eff.swap(bdup);
   stopwatch.add("copy");
}

void GpuDepondtIntegrator::copyThermalFieldFrom(const GpuDepondtIntegrator& source) {
   if(thermfield.getField().size() != source.thermfield.getField().size())
      throw std::invalid_argument("Depondt thermal-field shapes do not match");
   ASSERT_GPU(GPU_MEMCPY(const_cast<real*>(thermfield.getField().data()),
                         source.thermfield.getField().data(),
                         thermfield.getField().size() * sizeof(real),
                         GPU_MEMCPY_DEVICE_TO_DEVICE));
}

namespace {

unsigned int checkedActiveAtomCount(const int* oneBasedAtoms, std::size_t activeAtomCount,
                                    const GpuTensor<real, 3>& emom) {
   if(activeAtomCount != 0 && oneBasedAtoms == nullptr)
      throw std::invalid_argument("Depondt active atom list is null");
   const std::size_t sites = static_cast<std::size_t>(emom.extent(1));
   if(activeAtomCount > sites || activeAtomCount > std::numeric_limits<unsigned int>::max())
      throw std::invalid_argument("Depondt active atom count is outside the production lattice");
   return static_cast<unsigned int>(activeAtomCount);
}

} // namespace

void GpuDepondtIntegrator::evolveFirstWithThermalField(
   deviceLattice& gpuLattice, const int* oneBasedAtoms,
   std::size_t activeAtomCount) {
   const unsigned int count = checkedActiveAtomCount(oneBasedAtoms, activeAtomCount,
                                                      gpuLattice.emom);
   if(count == 0) return;
   parallel.gpuActiveAtomCall(EvolveFirst(mrod, blocal, bdup, gpuLattice.emomM,
                                          gpuLattice.emom, gpuLattice.beff,
                                          thermfield.getField(), gpuLattice.btorque,
                                          gpuLattice.mmom, timestep, gamma, damping,
                                          stt != 'N'), oneBasedAtoms, count);
   parallel.gpuActiveAtomCall(CommitActivePredictor(gpuLattice.emom, gpuLattice.emom2,
                                                     gpuLattice.b2eff, mrod, bdup),
                              oneBasedAtoms, count);
}

void GpuDepondtIntegrator::evolveFirst(deviceLattice& gpuLattice,
                                       const int* oneBasedAtoms,
                                       std::size_t activeAtomCount) {
   const unsigned int count = checkedActiveAtomCount(oneBasedAtoms, activeAtomCount,
                                                      gpuLattice.emom);
   stopwatch.skip();

   // Keep the production generator and its full (site,ensemble) draw -- the
   // curand/hiprand call still consumes exactly 3*N*M values from the shared
   // stream, every step, so the generator sequence is byte-for-byte
   // unchanged from the feature-off path and compact-list order still cannot
   // change a site's random draw (CGP-06A Strategy 1). Only the downstream
   // field *write* is scoped to the active list: every reader of
   // thermalField below is itself restricted to the same oneBasedAtoms/count
   // list, so an unwritten entry is never read.
   thermfield.randomize(gpuLattice.mmom, oneBasedAtoms, count);
   stopwatch.add("thermfield");

   parallel.gpuActiveAtomCall(EvolveFirst(mrod, blocal, bdup, gpuLattice.emomM,
                                          gpuLattice.emom, gpuLattice.beff,
                                          thermfield.getField(), gpuLattice.btorque,
                                          gpuLattice.mmom, timestep, gamma, damping,
                                          stt != 'N'), oneBasedAtoms, count);
   stopwatch.add("predictor");

   parallel.gpuActiveAtomCall(CommitActivePredictor(gpuLattice.emom, gpuLattice.emom2,
                                                     gpuLattice.b2eff, mrod, bdup),
                              oneBasedAtoms, count);
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

void GpuDepondtIntegrator::evolveSecond(deviceLattice& gpuLattice,
                                        const int* oneBasedAtoms,
                                        std::size_t activeAtomCount) {
   const unsigned int count = checkedActiveAtomCount(oneBasedAtoms, activeAtomCount,
                                                      gpuLattice.emom);
   stopwatch.skip();

   parallel.gpuActiveAtomCall(
      GpuCommon::AddAtoms(blocal, gpuLattice.beff, thermfield.getField()),
      oneBasedAtoms, count);
   stopwatch.add("localfield");

   parallel.gpuActiveAtomCall(BuildEffectiveField(bdup, blocal, gpuLattice.emom, damping),
                              oneBasedAtoms, count);
   if(stt != 'N')
      parallel.gpuActiveAtomCall(GpuCommon::AddToAtoms(bdup, gpuLattice.btorque),
                                 oneBasedAtoms, count);
   stopwatch.add("buildbeff");

   parallel.gpuActiveAtomCall(GpuCommon::AvgAtoms(bdup, gpuLattice.b2eff),
                              oneBasedAtoms, count);
   parallel.gpuActiveAtomCall(RestoreActiveInitial(gpuLattice.emom, gpuLattice.emom2),
                              oneBasedAtoms, count);
   stopwatch.add("corrfield");

   parallel.gpuActiveAtomCall(Rotate(mrod, gpuLattice.emom, bdup, timestep, gamma, damping),
                              oneBasedAtoms, count);
   stopwatch.add("rotate");

   parallel.gpuActiveAtomCall(CommitActiveCorrector(gpuLattice.emom2, mrod),
                              oneBasedAtoms, count);
   stopwatch.add("copy");
}

void GpuDepondtIntegrator::rotate(const GpuTensor<real, 3>& emom, real delta_t) {
   parallel.gpuAtomCall(Rotate(mrod, emom, bdup, delta_t, gamma, damping));
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
