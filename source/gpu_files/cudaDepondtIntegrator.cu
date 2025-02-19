#include <cuda.h>
#include <curand.h>

#include "c_headers.hpp"
#include "cudaCommon.hpp"
#include "cudaDepondtIntegrator.hpp"
#include "tensor.cuh"
#include "printDebug.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "cudaStructures.hpp"

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

class CudaDepondtIntegrator::BuildEffectiveField : public CudaParallelizationHelper::Atom {
private:
   real* bdup;
   const real* blocal;
   const real* emom;
   real damping;

public:
   BuildEffectiveField(CudaTensor<real, 3>& p1, const CudaTensor<real, 3>& p2, const CudaTensor<real, 3>& p3, real p4) {
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

class CudaDepondtIntegrator::Rotate : public CudaParallelizationHelper::Atom {
private:
   real* mrod;
   const real* emom;
   const real* bdup;
   real timestep;
   real gamma;
   real damping;

public:
   Rotate(CudaTensor<real, 3>& p1, const CudaTensor<real, 3>& p2, const CudaTensor<real, 3>& p3, real p4, real p5, real p6) {
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
      real u = 1 - cosv;

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

////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

// Constructor
CudaDepondtIntegrator::CudaDepondtIntegrator()
    : stopwatch(GlobalStopwatchPool::get("Cuda Depondt integrator")),
      parallel(CudaParallelizationHelper::def) {
}

// Destructor
CudaDepondtIntegrator::~CudaDepondtIntegrator() {
   release();
}

// Initiator
bool CudaDepondtIntegrator::initiate(const SimulationParameters SimParam) {
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
   mrod.Allocate(3, SimParam.N, SimParam.M);
   blocal.Allocate(3, SimParam.N, SimParam.M);
   bdup.Allocate(3, SimParam.N, SimParam.M);

   // All initialized?
   if(cudaDeviceSynchronize() != cudaSuccess) {
      release();
      return false;
   }

   return true;
}

bool CudaDepondtIntegrator::initiateConstants(const SimulationParameters SimParam, const Tensor<real, 1>temperature) {
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


void CudaDepondtIntegrator::resetConstants(const Tensor<real, 1> temperature) {
   // Set parameters
   thermfield.resetConstants(temperature, timestep, gamma, k_bolt, mub, damping); 
}

// Releaser
void CudaDepondtIntegrator::release() {
   mrod.Free();
   blocal.Free();
   bdup.Free();
}

// First step of Depond solver, calculates the stochastic field and rotates the
// magnetic moments according to the effective field
// Dupont recipe J. Phys.: Condens. Matter 21 (2009) 336005
void CudaDepondtIntegrator::evolveFirst(cudaLattice& gpuLattice) {
   // Timing
   stopwatch.skip();

   //_dpr;

   // Randomize stochastic field
   thermfield.randomize(gpuLattice.mmom);
   stopwatch.add("thermfield");

   CudaCommon::Add cm(blocal, gpuLattice.beff, thermfield.getField());

   //_dpr;
   // Construct local field
   parallel.cudaElementCall(CudaCommon::Add(blocal, gpuLattice.beff, thermfield.getField()));
   stopwatch.add("localfield");

   //_dpr;
   // Construct effective field (including damping term)
   buildbeff(gpuLattice.emom, gpuLattice.btorque);
   stopwatch.add("buildbeff");

   //_dpr;
   // Set up rotation matrices and perform rotations
   rotate(gpuLattice.emom, timestep);
   stopwatch.add("rotate");

   // copy m(t) to emom2 and m(t+dt) to emom for heisge, save b(t)
   parallel.cudaElementCall(CudaCommon::ScalarMult(gpuLattice.emomM, mrod, gpuLattice.mmom));
   gpuLattice.emom2.swap(gpuLattice.emom);  // Previous emom will not be needed
   gpuLattice.emom.swap(mrod);   // Previous mrod will not be needed
   gpuLattice.b2eff.swap(bdup);  // Previous bdup will not be needed
   stopwatch.add("copy");
}

// Second step of Depond solver, calculates the corrected effective field from
// the predicted effective fields. Rotates the moments in the corrected field
void CudaDepondtIntegrator::evolveSecond(cudaLattice& gpuLattice) {
   // Timing
   stopwatch.skip();

   // Construct local field
   parallel.cudaElementCall(CudaCommon::Add(blocal, gpuLattice.beff, thermfield.getField()));
   stopwatch.add("localfield");

   // Construct effective field (including damping term)
   buildbeff(gpuLattice.emom, gpuLattice.btorque);
   stopwatch.add("buildbeff");

   // Corrected field
   parallel.cudaElementCall(CudaCommon::Avg(bdup, gpuLattice.b2eff));
   gpuLattice.emom.swap(gpuLattice.emom2);  // Vaild as emom2 wont be used after its coming swap
   stopwatch.add("corrfield");

   // Final rotation
   rotate(gpuLattice.emom, timestep);
   stopwatch.add("rotate");

   // Swap
   gpuLattice.emom2.swap(mrod);  // Vaild as mrod wont be needed after this
   stopwatch.add("copy");
}

void CudaDepondtIntegrator::rotate(const CudaTensor<real, 3>& emom, real delta_t) {
   parallel.cudaAtomCall(Rotate(mrod, emom, bdup, timestep, gamma, damping));
}

// Constructs the effective field (including damping term)
void CudaDepondtIntegrator::buildbeff(const CudaTensor<real, 3>& emom,
                                      const CudaTensor<real, 3>& btorque) {
   parallel.cudaAtomCall(BuildEffectiveField(bdup, blocal, emom, damping));

   // TODO untested
   if(stt != 'N') {
      parallel.cudaElementCall(CudaCommon::AddTo(bdup, btorque));
   }
}

