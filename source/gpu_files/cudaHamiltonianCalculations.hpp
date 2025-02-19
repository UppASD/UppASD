#pragma once

#include <cuda_runtime.h>
#include "tensor.cuh"
#include "cudaParallelizationHelper.hpp"
#include "real_type.h"
#include "cudaStructures.hpp"


class CudaHamiltonianCalculations {
private:
   // Local matrices
  struct Exchange {
      unsigned int mnn;
      CudaTensor<real, 2> coupling;
      CudaTensor<unsigned int, 1> neighbourCount;
      CudaTensor<unsigned int, 2> neighbourPos;
      CudaTensor<unsigned int, 1> rHam;
   };

  struct DMinteraction {
      unsigned int mnn;
      CudaTensor<real, 3> interaction;
      CudaTensor<unsigned int, 1> neighbourCount;
      CudaTensor<unsigned int, 2> neighbourPos;
   };

   struct TensorialExchange {
      unsigned int mnn;
      CudaTensor<real, 4> tensor;
      CudaTensor<unsigned int, 1> neighbourCount;
      CudaTensor<unsigned int, 2> neighbourPos;
   };

   struct Anisotropy {
      CudaTensor<unsigned int, 1> taniso;
      CudaTensor<real, 2> eaniso;
      CudaTensor<real, 2> kaniso;
      CudaTensor<real, 1> sb;  // Ratio between uniaxial and cubic anisotropie
   };

   struct HamRed{
		CudaTensor<unsigned int, 1> redNeibourCount; //Reduced Hamiltonian -- shared between Jij, DMI anf Jtens
	} ;


   Exchange ex;
   DMinteraction dm;
   TensorialExchange tenEx;
   Anisotropy aniso;
   HamRed redHam;
   CudaTensor<real, 3> external_field;

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
   CudaParallelizationHelper& parallel;

public:
   // Parallelization helpers
   class SetupNeighbourList;
   class SetupNeighbourListDM;
   class SetupNeighbourListExchangeTensor;
   class SetupAnisotropy;
   class HeisgeJij;
   class HeisJijTensor;
   class HeisgeJijAniso;
   class HeisgeJijElement;

   // Constructor
   CudaHamiltonianCalculations();

   // Initiate
   bool initiate(const Flag Flags, const SimulationParameters SimParam, cudaHamiltonian& gpuHamiltonian);

   // Initiated
   bool isInitiated() {
      return initiated;
   }
   // Calculate "heisge"
   void heisge(cudaLattice& gpuLattice);
};
