#pragma once

#include <curand.h>

#include "c_headers.hpp"
#include "tensor.cuh"
#include "real_type.h"

struct Flag {    
   bool do_dm;
   bool do_jtensor;
   unsigned int do_aniso;
   bool do_avrg;
   bool do_cumu;    
   bool do_iphase_now;    
   bool do_mphase_now;    
};

struct SimulationParameters {    
   char stt;
   int SDEalgh;

   std::size_t rstep;
   std::size_t nstep;
   std::size_t N;
   std::size_t NH;
   std::size_t M;
   std::size_t mnn;
   std::size_t mnndm;
   unsigned int ipnphase;
   unsigned int ipmcnphase;
   unsigned int mcnstep;

   real delta_t;
   real gamma;
   real k_bolt;
   real mub;
   real damping;
   real Temp;

   std::size_t avrg_step;  
   std::size_t avrg_buff;
   std::size_t cumu_step;
   std::size_t cumu_buff;  

   int mompar;
   char initexc;
    // Thermfield parameters
   curandRngType_t rngType;
   unsigned long long randomSeed;

   const real* binderc;
   real* mavg;
};


struct hostHamiltonian {
   Tensor<unsigned int, 1>     aHam;                             //reduced Hamiltonian
   Tensor<real, 2>             ncoup;            //Jij
   Tensor<unsigned int, 2>     nlist;
   Tensor<unsigned int, 1>     nlistsize;
   Tensor<real, 3>             dmvect;     //DMI
   Tensor<unsigned int, 2>     dmlist;
   Tensor<unsigned int, 1>     dmlistsize;
   Tensor<real, 4>             j_tensor;
   Tensor<real, 2>             kaniso;
   Tensor<real, 2>             eaniso;
   Tensor<unsigned int, 1>     taniso;
   Tensor<real, 1>             sb;
   Tensor<real, 3>             extfield;

};
   
struct hostLattice {
   Tensor<real, 3> beff;
   Tensor<real, 3> b2eff;
   Tensor<real, 3> emomM;
   Tensor<real, 3> emom;
   Tensor<real, 3>  emom2;
   Tensor<real, 2>  mmom;
   Tensor<real, 2>  mmom0;
   Tensor<real, 2>  mmom2;
   Tensor<real, 2>  mmomi;
   Tensor<real, 3> btorque;
   Tensor<real, 1> temperature;
   Tensor<real, 1> ipTemp;
   Tensor<real, 2> ipTemp_array;
   Tensor<unsigned int, 1> ipnstep;
   Tensor<unsigned int, 1> ipmcnstep;

};

struct hostMeasurables {    
   Tensor<real, 1> mavg_buff;
   Tensor<real, 1> eavg_buff;
   Tensor<real, 1> mcumu_buff;
   Tensor<real, 1> ecumu_buff;
   Tensor<real, 1> binderc; 
};
   
struct cudaHamiltonian {
   CudaTensor<unsigned int, 1>     aHam;                             //reduced Hamiltonian
   CudaTensor<real, 2>             ncoup;            //Jij
   CudaTensor<unsigned int, 2>     nlist;
   CudaTensor<unsigned int, 1>     nlistsize;
   CudaTensor<real, 3>             dmvect;     //DMI
   CudaTensor<unsigned int, 2>     dmlist;
   CudaTensor<unsigned int, 1>     dmlistsize;
   CudaTensor<real, 4>             j_tensor;
   CudaTensor<real, 2>             kaniso;
   CudaTensor<real, 2>             eaniso;
   CudaTensor<unsigned int, 1>     taniso;
   CudaTensor<real, 1>             sb;
   CudaTensor<real, 3>             extfield;
}; 
   
struct cudaLattice {
   CudaTensor<real, 3> beff;
   CudaTensor<real, 3> b2eff;
   CudaTensor<real, 3> eneff;
   CudaTensor<real, 3> emomM;
   CudaTensor<real, 3> emom;
   CudaTensor<real, 3>  emom2;
   CudaTensor<real, 2>  mmom;
   CudaTensor<real, 2>  mmom0;
   CudaTensor<real, 2>  mmom2;
   CudaTensor<real, 2>  mmomi;
   CudaTensor<real, 3> btorque;
   CudaTensor<real, 1> temperature;
   CudaTensor<real, 1> ipTemp;
   CudaTensor<real, 1> ipTemp_array;
} ;

struct cudaMeasurables {    
   CudaTensor<real, 1> mavg_buff;
   CudaTensor<real, 1> eavg_buff;
   CudaTensor<real, 1> mcumu_buff;
   CudaTensor<real, 1> ecumu_buff;
   CudaTensor<real, 1> binderc;
};
   
