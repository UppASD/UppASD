#pragma once

#include <curand.h>

#include "cudaMatrix.hpp"
#include "fortMatrix.hpp"
#include "real_type.h"

class CudaMdSimulation {
private:
   char stt;
   int SDEalgh;

   std::size_t rstep;
   std::size_t nstep;
   std::size_t Natom;
   std::size_t Mensemble;
   std::size_t max_no_neigh;

   real delta_t;
   real gamma;
   real k_bolt;
   real mub;
   real damping;

   const real* binderc;
   real* mavg;

   int mompar;
   char initexc;

   bool do_dm;
   bool do_jtensor;
   unsigned int do_aniso;
   std::size_t max_no_dmneigh;

   // Thermfield parameters
   curandRngType_t rngType;
   unsigned long long randomSeed;

   // Matrix class wrappers for Fortran data.
   fortMatrix<real, 2> f_ncoup;
   fortMatrix<unsigned int, 2> f_nlist;
   fortMatrix<unsigned int, 1> f_nlistsize;
   fortMatrix<real, 3, 3> f_dmvect;
   fortMatrix<unsigned int, 2> f_dmlist;
   fortMatrix<unsigned int, 1> f_dmlistsize;
   fortMatrix<real, 3, 3> f_beff;
   fortMatrix<real, 3, 3> f_b2eff;
   fortMatrix<real, 3, 3> f_emomM;
   fortMatrix<real, 3, 3> f_emom;
   fortMatrix<real, 3, 3> f_emom2;
   fortMatrix<real, 3, 3> f_external_field;
   fortMatrix<real, 2> f_mmom;
   fortMatrix<real, 3, 3> f_btorque;
   fortMatrix<real, 1> f_temperature;
   fortMatrix<real, 2> f_mmom0;
   fortMatrix<real, 2> f_mmom2;
   fortMatrix<real, 2> f_mmomi;
   fortMatrix<real, 4, 3, 3> f_j_tensor;

   fortMatrix<real, 2, 2> f_kaniso;
   fortMatrix<real, 2, 3> f_eaniso;
   fortMatrix<unsigned int, 1> f_taniso;
   fortMatrix<real, 1> f_sb;

   cudaMatrix<real, 3, 3> beff;
   cudaMatrix<real, 3, 3> b2eff;
   cudaMatrix<real, 3, 3> emomM;
   cudaMatrix<real, 3, 3> emom;
   cudaMatrix<real, 3, 3> emom2;
   cudaMatrix<real, 3, 3> external_field;
   cudaMatrix<real, 2> mmom;
   cudaMatrix<real, 3, 3> btorque;
   //	cudaMatrix<real,1>   temperature;
   cudaMatrix<real, 2> mmom0;
   cudaMatrix<real, 2> mmom2;
   cudaMatrix<real, 2> mmomi;

   bool isInitiated;

   void printConstants();
   void printMdStatus(std::size_t mstep);

   void initiate_fortran();

public:
   CudaMdSimulation();
   ~CudaMdSimulation();

   void measurementPhase();

   void initiateConstants();
   bool initiateMatrices();

   void release();

   void copyFromFortran();
   void copyToFortran();
};

