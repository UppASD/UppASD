#pragma once

#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"

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

   struct MacroBlockLayout {
      unsigned int enabled = 0;
      unsigned int nblocks = 0;
      unsigned int npairGroups = 0;
      unsigned int nentries = 0;
      unsigned int maxBlockAtoms = 0;
      GpuTensor<unsigned int, 1> atomToBlock;
      GpuTensor<unsigned int, 1> blockAtomOffset;
      GpuTensor<unsigned int, 1> blockAtoms;
      GpuTensor<unsigned int, 1> pairGroupSrc;
      GpuTensor<unsigned int, 1> pairGroupDst;
      GpuTensor<unsigned int, 1> pairGroupEntryOffset;
      GpuTensor<unsigned int, 1> entryAtomI;
      GpuTensor<unsigned int, 1> entryAtomJ;
      GpuTensor<unsigned int, 1> entryIH;
      GpuTensor<unsigned int, 1> entryJslot;
      GpuTensor<unsigned int, 1> entryLocalI;
      GpuTensor<unsigned int, 1> entryLocalJ;
   };


   Exchange ex;
   DMinteraction dm;
   TensorialExchange tenEx;
   Anisotropy aniso;
   HamRed redHam;
   MacroBlockLayout macroblocks;
   GpuTensor<real, 3> external_field;

   bool do_j_tensor = false;
   bool do_dm = false;
   int do_aniso = 0;
   int do_ene = 0;
   bool use_macroblock_backend = false;
   bool use_macroblock_pair_backend = false;

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
   void heisge(deviceLattice& gpuLattice, deviceEnergies& gpuEnergies, bool measure);
};
