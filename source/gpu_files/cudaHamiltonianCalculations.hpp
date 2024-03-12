#pragma once

#include <cuda_runtime.h>

#include "cudaMatrix.hpp"
#include "cudaParallelizationHelper.hpp"
#include "hostMatrix.hpp"
#include "real_type.h"

class CudaHamiltonianCalculations {
private:
   // Local matrices
   typedef struct Exchange {
      unsigned int mnn;
      cudaMatrix<real, 2> coupling;
      cudaMatrix<unsigned int, 1> neighbourCount;
      cudaMatrix<unsigned int, 2> neighbourPos;
   } Exchange;

   typedef struct DMinteraction {
      unsigned int mnn;
      cudaMatrix<real, 3, 3> interaction;
      cudaMatrix<unsigned int, 1> neighbourCount;
      cudaMatrix<unsigned int, 2> neighbourPos;
   } DMinteraction;

   typedef struct TensorialExchange {
      unsigned int mnn;
      cudaMatrix<real, 4, 3, 3> tensor;
      cudaMatrix<unsigned int, 1> neighbourCount;
      cudaMatrix<unsigned int, 2> neighbourPos;
   } TensorialExchange;

   typedef struct Anisotropy {
      cudaMatrix<unsigned int, 1> taniso;
      cudaMatrix<real, 2, 3> eaniso;
      cudaMatrix<real, 2, 2> kaniso;
      cudaMatrix<real, 1> sb;  // Ratio between uniaxial and cubic anisotropie
   } Anisotropy;

   Exchange ex;
   DMinteraction dm;
   TensorialExchange tenEx;
   Anisotropy aniso;

   bool do_j_tensor = false;
   int do_aniso = 0;

   // Initiation flag
   bool initiated;

   // System size
   unsigned int N;

   // Parallelization helper
   CudaParallelizationHelper &parallel;

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
   bool initiate(const hostMatrix<real, 2> &ncoup, const hostMatrix<unsigned int, 2> &nlist,
                 const hostMatrix<unsigned int, 1> &nlistsize, const hostMatrix<real, 3, 3> &dm_ncoup,
                 const hostMatrix<unsigned int, 2> &dm_nlist, const hostMatrix<unsigned int, 1> &dm_nlistsize,
                 const int do_dm, const int do_j_tensor, const hostMatrix<real, 4, 3, 3> j_tensor,
                 const int do_aniso, const hostMatrix<real, 2, 2> kaniso, const hostMatrix<real, 2, 3> eaniso,
                 const hostMatrix<unsigned int, 1> taniso, const hostMatrix<real, 1> sb);

   // Initiated
   bool isInitiated() {
      return initiated;
   }

   // Release
   void release();

   // Calculate "heisge"
   void heisge(cudaMatrix<real, 3, 3> &beff, const cudaMatrix<real, 3, 3> &emomM,
               const cudaMatrix<real, 3, 3> &external_field);
};
