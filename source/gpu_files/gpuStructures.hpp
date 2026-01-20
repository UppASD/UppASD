#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"

#include "gpu_wrappers.h"
#include "measurementData.h"
#if defined(HIP_V)
#include <hiprand/hiprand.h>
#elif defined(CUDA_V)
#include <curand.h>
#endif


struct Flag {    
   bool do_dm;
   bool do_jtensor;
   unsigned int do_aniso;
   bool do_avrg;
   bool do_cumu;    
   bool do_iphase_now;    
   bool do_mphase_now;    
   char do_sc;
   bool do_gpu_correlations;
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

   std::size_t sc_sep;
   std::size_t sc_step;
   std::size_t nw;
   std::size_t nq;
   std::size_t sc_max_nstep;
   std::size_t sc_window_fun;



    // Thermfield parameters
#if defined(HIP_V)
   hiprandRngType_t rngType;
#elif defined(CUDA_V) 
   curandRngType_t rngType;
#endif

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

struct hostCorrelations {    
   Tensor<real, 2> coord;
   Tensor<real, 1> r_mid;
   Tensor<real, 2> q;
   Tensor<real, 1> w;
};
   
struct deviceHamiltonian {
   GpuTensor<unsigned int, 1>     aHam;                             //reduced Hamiltonian
   GpuTensor<real, 2>             ncoup;            //Jij
   GpuTensor<unsigned int, 2>     nlist;
   GpuTensor<unsigned int, 1>     nlistsize;
   GpuTensor<real, 3>             dmvect;     //DMI
   GpuTensor<unsigned int, 2>     dmlist;
   GpuTensor<unsigned int, 1>     dmlistsize;
   GpuTensor<real, 4>             j_tensor;
   GpuTensor<real, 2>             kaniso;
   GpuTensor<real, 2>             eaniso;
   GpuTensor<unsigned int, 1>     taniso;
   GpuTensor<real, 1>             sb;
   GpuTensor<real, 3>             extfield;
}; 

struct deviceLattice {
   GpuTensor<real, 3> beff;
   GpuTensor<real, 3> b2eff;
   GpuTensor<real, 3> eneff;
   GpuTensor<real, 3> emomM;
   GpuTensor<real, 3> emom;
   GpuTensor<real, 3>  emom2;
   GpuTensor<real, 2>  mmom;
   GpuTensor<real, 2>  mmom0;
   GpuTensor<real, 2>  mmom2;
   GpuTensor<real, 2>  mmomi;
   GpuTensor<real, 3> btorque;
   GpuTensor<real, 1> temperature;
   GpuTensor<real, 1> ipTemp;
   GpuTensor<real, 1> ipTemp_array;
   GpuVector<EnergyData> energy;
};

struct deviceMeasurables {    
   GpuTensor<real, 1> mavg_buff;
   GpuTensor<real, 1> eavg_buff;
   GpuTensor<real, 1> mcumu_buff;
   GpuTensor<real, 1> ecumu_buff;
   GpuTensor<real, 1> binderc;
};
   
