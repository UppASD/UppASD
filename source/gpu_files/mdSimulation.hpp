#pragma once

#include "fortMatrix.hpp"
#include "real_type.h"

class MdSimulation {
private:
   char stt;
   int SDEalgh;

   usd_int rstep;
   usd_int nstep;
   usd_int Natom;
   usd_int Mensemble;
   usd_int max_no_neigh;

   real delta_t;
   real gamma;
   real k_bolt;
   real mub;
   real damping;

   const real* binderc;
   real* mavg;

   int mompar;
   char initexc;

   int do_dm;
   usd_int max_no_dmneigh;

   // <real,3,3> -- real precision, 3 dimensions, first dim always 3
   fortMatrix<real, 2> ncoup;
   fortMatrix<unsigned int, 2> nlist;
   fortMatrix<unsigned int, 1> nlistsize;
   fortMatrix<real, 3, 3> dmvect;
   fortMatrix<unsigned int, 2> dmlist;
   fortMatrix<unsigned int, 1> dmlistsize;
   fortMatrix<real, 3, 3> beff;
   fortMatrix<real, 3, 3> b2eff;
   fortMatrix<real, 3, 3> emomM;
   fortMatrix<real, 3, 3> emom;
   fortMatrix<real, 3, 3> emom2;
   fortMatrix<real, 3, 3> external_field;
   fortMatrix<real, 2> mmom;
   fortMatrix<real, 3, 3> btorque;
   fortMatrix<real, 1> temperature;
   fortMatrix<real, 2> mmom0;
   fortMatrix<real, 2> mmom2;
   fortMatrix<real, 2> mmomi;

   void printConstants();

   bool isOwnData;
   void freeOwn();

public:
   ~MdSimulation();

   void measurementPhase();

   void initiateConstants();
   void initiateFortran();
   void initiateOwn();

   void copyFromFortran();
   void copyToFortran();
};


