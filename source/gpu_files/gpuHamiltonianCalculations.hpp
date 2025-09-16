#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#if defined(HIP_V)
#include <hip/hip_runtime.h>
#elif defined(CUDA_V)
#include <cuda_runtime.h>
#endif
using ParallelizationHelper = GpuParallelizationHelper;



class GpuHamiltonianCalculations {
private:
   // Local matrices
  struct Exchange {
      unsigned int mnn;
      GpuTensor<real, 2> coupling;
      GpuTensor<unsigned int, 1> neighbourCount;
      GpuTensor<unsigned int, 2> neighbourPos;
      GpuTensor<unsigned int, 1> rHam;
   };

  struct DMinteraction {
      unsigned int mnn;
      GpuTensor<real, 3> interaction;
      GpuTensor<unsigned int, 1> neighbourCount;
      GpuTensor<unsigned int, 2> neighbourPos;
   };

   struct TensorialExchange {
      unsigned int mnn;
      GpuTensor<real, 4> tensor;
      GpuTensor<unsigned int, 1> neighbourCount;
      GpuTensor<unsigned int, 2> neighbourPos;
   };

   struct Anisotropy {
      GpuTensor<unsigned int, 1> taniso;
      GpuTensor<real, 2> eaniso;
      GpuTensor<real, 2> kaniso;
      GpuTensor<real, 1> sb;  // Ratio between uniaxial and cubic anisotropie
   };

   struct HamRed{
		GpuTensor<unsigned int, 1> redNeibourCount; //Reduced Hamiltonian -- shared between Jij, DMI anf Jtens
	} ;


   Exchange ex;
   DMinteraction dm;
   TensorialExchange tenEx;
   Anisotropy aniso;
   HamRed redHam;
   GpuTensor<real, 3> external_field;

   bool do_j_tensor = false;
   bool do_dm = false;
   int do_aniso = 0;

   // Initiation flag
   bool initiated;

   // System size
   unsigned int N;
   unsigned int NH;
   unsigned int mnn;
   unsigned int mnndm;

   // Parallelization helper
   ParallelizationHelper& parallel;

public:
   // Parallelization helpers
   class SetupNeighbourList;
   class SetupNeighbourListDM;
   class SetupNeighbourListExchangeTensor;
   class SetupAnisotropy;
   class HeisgeJij;
   class HeisgeJijDM;
   class HeisgeJijAniso;
   class HeisgeJijDMAniso;
   class HeisgeJijTensor;
   class HeisgeJijTensorAniso;
   class HeisgeJijElement;

   // Constructor
   GpuHamiltonianCalculations();

   // Initiate
   bool initiate(const Flag Flags, const SimulationParameters SimParam, deviceHamiltonian& gpuHamiltonian);

   // Initiated
   bool isInitiated() {
      return initiated;
   }
   // Calculate "heisge"
   void heisge(deviceLattice& gpuLattice);
};
