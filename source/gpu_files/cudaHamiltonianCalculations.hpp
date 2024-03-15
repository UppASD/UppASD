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
      usd_int mnn;
      cudaMatrix<real, 2> coupling;
      cudaMatrix<usd_int, 1> neighbourCount;
      cudaMatrix<usd_int, 2> neighbourPos;
   } Exchange;

   typedef struct DMinteraction {
      usd_int mnn;
      cudaMatrix<real, 3, 3> interaction;
      cudaMatrix<usd_int, 1> neighbourCount;
      cudaMatrix<usd_int, 2> neighbourPos;
   } DMinteraction;

   typedef struct TensorialExchange {
      usd_int mnn;
      cudaMatrix<real, 4, 3, 3> tensor;
      cudaMatrix<usd_int, 1> neighbourCount;
      cudaMatrix<usd_int, 2> neighbourPos;
   } TensorialExchange;

   typedef struct Anisotropy {
      cudaMatrix<usd_int, 1> taniso;
      cudaMatrix<real, 2, 3> eaniso;
      cudaMatrix<real, 2, 2> kaniso;
      cudaMatrix<real, 1> sb;  // Ratio between uniaxial and cubic anisotropy
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
   usd_int N;

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
   bool initiate(const hostMatrix<real, 2>& ncoup, const hostMatrix<usd_int, 2>& nlist,
                 const hostMatrix<usd_int, 1>& nlistsize, const hostMatrix<real, 3, 3>& dm_ncoup,
                 const hostMatrix<usd_int, 2>& dm_nlist, const hostMatrix<usd_int, 1>& dm_nlistsize,
                 const int do_dm, const int do_j_tensor, const hostMatrix<real, 4, 3, 3> j_tensor,
                 const int do_aniso, const hostMatrix<real, 2, 2> kaniso, const hostMatrix<real, 2, 3> eaniso,
                 const hostMatrix<usd_int, 1> taniso, const hostMatrix<real, 1> sb);

   // Initiated
   bool isInitiated() {
      return initiated;
   }

   // Release
   void release();

   // Calculate "heisge"
   void heisge(cudaMatrix<real, 3, 3>& beff, const cudaMatrix<real, 3, 3>& emomM,
               const cudaMatrix<real, 3, 3>& external_field);
};
