// Routines for updating magnetic moment after time evolution

#include "c_headers.hpp"
#include "gpuCommon.hpp"
#include "tensor.hpp"
#include "gpuMomentUpdater.hpp"

#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "gpuStructures.hpp"
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"
#if defined (HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda.h>
#endif

using ParallelizationHelper = GpuParallelizationHelper;

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

class GpuMomentUpdater::Mompar1 : public ParallelizationHelper::Atom {
private:
   real* mmom2;
   const real* mmom0;
   const real* emom2;

public:
   Mompar1(GpuTensor<real, 2>& p1, const GpuTensor<real, 2>& p2, const GpuTensor<real, 3>& p3) {
      mmom2 = p1.data();
      mmom0 = p2.data();
      emom2 = p3.data();
   }

   __device__ void each(unsigned int atom) {
      mmom2[atom] = fmax((real)1e-4, mmom0[atom] * fabs(emom2[atom * 3 + 2]));
   }
};

class GpuMomentUpdater::Mompar2 : public ParallelizationHelper::Atom {
private:
   real* mmom2;
   const real* mmom0;
   const real* emom2;

public:
   Mompar2(GpuTensor<real, 2>& p1, const GpuTensor<real, 2>& p2, const GpuTensor<real, 3>& p3) {
      mmom2 = p1.data();
      mmom0 = p2.data();
      emom2 = p3.data();
   }

   __device__ void each(unsigned int atom) {
      real mz = emom2[atom * 3 + 2];
      mmom2[atom] = fmax((real)1e-4, mmom0[atom] * mz * mz);
   }
};

// mmomi = 1.0 / mmom
// mmomM = emom*  mmom
class GpuMomentUpdater::Copy1 : public ParallelizationHelper::Atom {
private:
   real* mmomi;
   real* emomM;
   const real* mmom;
   const real* emom;

public:
   Copy1(GpuTensor<real, 2>& p1, GpuTensor<real, 3>& p2, const GpuTensor<real, 2>& p3, const GpuTensor<real, 3>& p4) {
      mmomi = p1.data();
      emomM = p2.data();
      mmom = p3.data();
      emom = p4.data();
   }

   __device__ void each(unsigned int atom) {
      real m = mmom[atom];
      mmomi[atom] = 1 / m;

      real* my_emomM = &emomM[atom * 3];
      const real* my_emom = &emom[atom * 3];
      my_emomM[0] = m * my_emom[0];
      my_emomM[1] = m * my_emom[1];
      my_emomM[2] = m * my_emom[2];
   }
};

// mmomi = (mmom < 0.000001) ? 1 : (1.0 / mmom)
// mmomM = emom*  mmom
class GpuMomentUpdater::Copy2 : public ParallelizationHelper::Atom {
private:
   real* mmomi;
   real* emomM;
   const real* mmom;
   const real* emom;

public:
   Copy2(GpuTensor<real, 2>& p1, GpuTensor<real, 3>& p2, const GpuTensor<real, 2>& p3, const GpuTensor<real, 3>& p4) {
      mmomi = p1.data();
      emomM = p2.data();
      mmom = p3.data();
      emom = p4.data();
   }

   __device__ void each(unsigned int atom) {
      real m = mmom[atom];
      mmomi[atom] = (m < (real)0.000001) ? 1 : (1 / m);

      real* my_emomM = &emomM[atom * 3];
      const real* my_emom = &emom[atom * 3];
      my_emomM[0] = m * my_emom[0];
      my_emomM[1] = m * my_emom[1];
      my_emomM[2] = m * my_emom[2];
   }
};

////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

// Constructor
GpuMomentUpdater::GpuMomentUpdater(deviceLattice& gpuLattice, int p8, char p9)
    : mmom(gpuLattice.mmom),
      mmom0(gpuLattice.mmom0),
      mmom2(gpuLattice.mmom2),
      emom(gpuLattice.emom),
      emom2(gpuLattice.emom2),
      emomM(gpuLattice.emomM),
      mmomi(gpuLattice.mmomi),
      mompar(p8),
      initexc(p9),
      stopwatch(GlobalStopwatchPool::get("Cuda moment")),
      parallel(ParallelizationHelperInstance) {

   // Exit if mompar is not supported
   if(mompar == 3) {
      std::fprintf(stderr, "mompar 3 (ptnanowire) not implemented!\n");
      std::exit(EXIT_FAILURE);
   }
}

// Wrapper routine for updating the magnetic moments
void GpuMomentUpdater::update() {
   // Timing
   stopwatch.skip();
   // Calculate
   switch(mompar) {
      case 0:  // mmom2 = mmom
         mmom2.swap(mmom);
         break;
      case 1:  // mmom2 = abs(emom_z)
         parallel.gpuAtomCall(Mompar1(mmom2, mmom0, emom2));
         break;
      case 2:  // mmom2 = abs(emom_z^2)
         parallel.gpuAtomCall(Mompar2(mmom2, mmom0, emom2));
         break;
   }
   stopwatch.add("calculate");

   // Copy
   emom.swap(emom2);
   mmom.swap(mmom2);
   stopwatch.add("copy - A");

   // Invert
   if(initexc != 'I') {
      parallel.gpuAtomCall(Copy1(mmomi, emomM, mmom, emom));
   } else {
      parallel.gpuAtomCall(Copy2(mmomi, emomM, mmom, emom));
   }
   stopwatch.add("copy - B");
}

