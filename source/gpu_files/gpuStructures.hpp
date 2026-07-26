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
   bool do_gpu_measurements;
   char do_sc_proj;
   char do_sc_projch;
   int do_ene;
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
   std::size_t N1 = 0;
   std::size_t N2 = 0;
   std::size_t N3 = 0;
   std::size_t NA = 0;
   bool do_gpu_convolution = false;
   char BC1 = '0';
   char BC2 = '0';
   char BC3 = '0';
   const real* C1 = nullptr;
   const real* C2 = nullptr;
   const real* C3 = nullptr;
   const real* Bas = nullptr;
   // Kept in fp64 even in a single-precision device build: it is part of the
   // physical dipole contract, not a device-field scratch value.
   double alat = 0.0;
   double gpu_dipole_tol = 1.0e-10;
   unsigned int ipnphase;
   unsigned int ipmcnphase;
   unsigned int mcnstep;

   real delta_t;
   real gamma;
   real k_bolt;
   real mub;
   real mry;
   real damping;
   real Temp;

   std::size_t avrg_step;  
   std::size_t avrg_buff;
   std::size_t cumu_step;
   std::size_t cumu_buff;  
   std::size_t ene_step;
   std::size_t ene_buff;  

   int mompar;
   char initexc;

   std::size_t sc_sep;
   std::size_t sc_step;
   std::size_t nw;
   std::size_t nq;
   std::size_t sc_max_nstep;
   std::size_t sc_window_fun;
   std::size_t NT;
   std::size_t Nchmax;



    // Thermfield parameters
#if defined(HIP_V)
   hiprandRngType_t rngType;
#elif defined(CUDA_V) 
   curandRngType_t rngType;
#endif

   unsigned long long randomSeed;

   const fortran_real* binderc;
   fortran_real* mavg;
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
   Tensor<unsigned int, 1>     macro_cell_index;
   Tensor<unsigned int, 1>     macro_nlistsize;
   Tensor<real, 2>             macro_center;
   Tensor<real, 2>             macro_min_coord;
   Tensor<real, 2>             macro_max_coord;

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
   Tensor<real, 1> ipdelta_t;
   Tensor<real, 1> iplambda1;
   Tensor<unsigned int, 1> ipmcnstep;

};

struct hostMeasurables {    
   Tensor<real, 1> mavg_buff;
   Tensor<real, 1> eavg_buff;
   Tensor<real, 1> mcumu_buff;
   Tensor<real, 1> ecumu_buff;
   Tensor<real, 1> binderc; 
};

//TODO: do we need?
struct hostEnergies {    
   Tensor<real, 1> tot;
   Tensor<real, 1> xc;
   Tensor<real, 1> dm;
   Tensor<real, 1> ani;
   Tensor<real, 1> ext;
   Tensor<real, 1> pair;
};

struct hostCorrelations {    
   Tensor<real, 2> coord;
   Tensor<real, 1> r_mid;
   Tensor<real, 2> q;
   Tensor<real, 1> w;
   Tensor<cpu_complex, 2> m_k;
   Tensor<cpu_complex, 3> m_kt;
   Tensor<cpu_complex, 3> m_kw;
   Tensor<real, 1> deltat_corr;  // Per-sample delta_t timing array
   Tensor<real, 1> scstep_arr;   // Per-sample sc_step array
   int sc_nsamp;  // Number of samples from GPU correlations
   int sc_tidx;   // Number of time steps accumulated in GPU correlations
   Vector<int> atype;
   Vector<int> achtype;
   Tensor<cpu_complex, 3> m_k_proj;//TODO: check dimentions
   Tensor<cpu_complex, 4> m_kt_proj;
   Tensor<cpu_complex, 4> m_kw_proj;
   Tensor<cpu_complex, 3> m_k_projch;//TODO: check dimentions
   Tensor<cpu_complex, 4> m_kt_projch;
   Tensor<cpu_complex, 4> m_kw_projch;

};
   
struct deviceHamiltonian {
   bool neighbourListsPrepared = false;
   GpuTensor<unsigned int, 1>     aHam;                             //reduced Hamiltonian
   GpuTensor<real, 2>             ncoup;            //Jij
   GpuTensor<unsigned int, 2>     nlist;
   GpuTensor<unsigned int, 1>     nlistsize;
   GpuTensor<real, 3>             dmvect;     //DMI
   GpuTensor<unsigned int, 2>     dmlist;
   GpuTensor<unsigned int, 1>     dmlistsize;
   // Device layout: [reduced site, output axis, input axis, neighbour].
   GpuTensor<real, 4>             j_tensor;
   GpuTensor<real, 2>             kaniso;
   GpuTensor<real, 2>             eaniso;
   GpuTensor<unsigned int, 1>     taniso;
   GpuTensor<real, 1>             sb;
   GpuTensor<real, 3>             extfield;
   GpuTensor<unsigned int, 1>     macro_cell_index;
   GpuTensor<unsigned int, 1>     macro_nlistsize;
   GpuTensor<real, 2>             macro_center;
   GpuTensor<real, 2>             macro_min_coord;
   GpuTensor<real, 2>             macro_max_coord;
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
   //GpuVector<EnergyData> energy;
};

struct deviceMeasurables {    
   GpuTensor<real, 1> mavg_buff;
   GpuTensor<real, 1> eavg_buff;
   GpuTensor<real, 1> mcumu_buff;
   GpuTensor<real, 1> ecumu_buff;
   GpuTensor<real, 1> binderc;
};

struct deviceEmpties {   
 
   GpuTensor<real, 1> real1D;
   GpuTensor<real, 2> real2D;
   GpuTensor<real, 3> real3D;
   GpuTensor<real, 4> real4D;

};

struct deviceEnergies {
   //GpuTensor<real, 1> exchM;
   //GpuTensor<real, 1> aniM;
   //GpuTensor<real, 1> dmM;
   //GpuTensor<real, 1> totalM;
   //GpuTensor<real, 1> tensorM;
   //GpuTensor<real, 1> extM;

   // 0 - exch, 1 - ani, 2 - dm, 3 - tensor, 4 - external, 5 - total,
   // 6 - periodic FFT dipole.  All entries are per atom in T*mu_B before
   // the common measurement conversion to mRy/atom.
   GpuTensor<real,2> energyM;

};

   
