#pragma once

#include "c_headers.hpp"
#include "real_type.h"
//#include <thrust/complex.h>
#include <complex>


// NAME         TYPE     DIMENSION   DESCRIPTION
//
// mompar       int       1          Parametrization of magnetic moment magnitudes (0=no)
// initexc      char      1          Mode of excitation of initial magnetic moments (I=vacancies, R=two magnon
// Raman, F=no)
//
// emom         real     (3,N,M)     Current unit moment vector
// emom2        real     (3,N,M)     Final (or temporary) unit moment vector
// emomM        real     (3,N,M)     Current magnetic moment vector
// mmom         real     (N,M)       Magnitude of magnetic moments
// mmom0        real     (N,M)       Starting magnitude of magnetic moments
// mmom2        real     (N,M)       Temporary value of magnitude of magnetic moments
// mmomi        real     (N,M)       Inverse of magnitude of magnetic moments
//
// mrod         real     (3,N,M)     Rotated magnetic moments
// btherm       real     (3,N,M)     Thermal stochastic field
// bloc         real     (3,N,M)     Local effective field
// bdup         real     (3,N,M)     Resulting effective field
//
// beff         real     (3,N,M)     Total effective field from application of Hamiltonian
// b2eff        real     (3,N,M)     Temporary storage of magnetic field
// btorque      real     (3,N,M)     Spin transfer torque
// emom         real     (3,N,M)     Current unit moment vector
// emom2        real     (3,N,M)     Final (or temporary) unit moment vector
// emomM        real     (3,N,M)     Current magnetic moment vector
// mmom         real     (N,M)       Magnitude of magnetic moments
// delta_t      real      1          Time step
// temperature  real     (N)         Temperature
//
// stt          char      1          Method to handle spin transfer torque
// sb 			real 	 (N)	     Ratio between cubic and uniaxial anisotropy
//
// dxyz_vec
// dxyz_atom
// dxyz_list

class FortranData {
public:
   // Scalars
   static char* stt;
   static int* SDEalgh;

   static unsigned int* rstep;
   static unsigned int* nstep;
   static unsigned int* Natom;
   static unsigned int* nHam;
   static unsigned int* Mensemble;
   static unsigned int* max_no_neigh;
   static unsigned int* ipmcnphase;
   static unsigned int* mcnstep;
   static unsigned int* ipnphase;
   

   static unsigned int* nq;
   static unsigned int* sc_step;
   static unsigned int* sc_sep;
   static char* do_sc;
   static unsigned int* sc_max_nstep;
   static unsigned int* sc_window_fun;
   static unsigned int* nw;

   static real* r_mid;
   static real* q;
   static real* coord;
   static real* w;
   static std::complex<real> * m_k;
   static std::complex<real>* m_kw;
   static std::complex<real>* m_kt;
   static real* deltat_corr;
   static real* scstep_arr;
   static int* sc_nsamp_ptr;
   static int* sc_tidx_ptr;


   //delta_t;

   static real* delta_t;
   static real* gamma;
   static real* k_bolt;
   static real* mub;
   static real* damping;
   static real * Temp;

   static real* binderc;
   static real* mavg;

   static int* mompar;
   static char* initexc;

   static unsigned int* do_dm;
   static unsigned int* max_no_dmneigh;


   static unsigned int*do_jtensor;  // Information on weather the exchange coupling tensor should be used or not
   static unsigned int* do_aniso;  // Information on weather the anisotropy should be used or not


   static char* do_cuda_measurements;           // Do measurements in CUDA (Y/N)
   static char* do_avrg;                        // Measure average magnetization (Y/N)
   static char* do_cumu;                        // Measure Binder cumulant, susceptibility, and specific heat(Y/N)
   static char* do_autocorr;                    // Perform autocorrelation (Y/N)
   static unsigned int* plotenergy;             // Calculate and plot energy (0/1)
   static char* do_skyno;
   static char* do_gpu_correlations;
   static char* do_gpu_macroblocks;


   // Matrices / vectors
   static unsigned int * aHam;

   static real* ncoup;
   static unsigned int* nlist;
   static unsigned int* nlistsize;

   static real* dmvect;
   static unsigned int* dmlist;
   static unsigned int* dmlistsize;

   static real* j_tensor;

   static real* kaniso;
   static real* eaniso;
   static unsigned int* taniso;
   static real* sb;

   static real* beff;
   static real* b2eff;
   static real* emomM;
   static real* emom;
   static real* emom2;
   static real* external_field;
   static real* mmom;
   static real* btorque;
   static real* temperature;
   static real* ipTemp;
   static unsigned int* ipmcnstep;
   static real* ipTemp_array;
   static unsigned int* ipnstep;
   static real* ipdelta_t;
   static real* iplambda1;
   static real* mmom0;
   static real* mmom2;
   static real* mmomi;

   static real* dxyz_vec;
   static int* dxyz_atom;
   static int* dxyz_list;
   static unsigned int* macroblock_enabled;
   static unsigned int* macroblock_nblocks;
   static unsigned int* macroblock_npair_groups;
   static unsigned int* macroblock_nentries;
   static unsigned int* macroblock_max_block_atoms;
   static unsigned int* macroblock_atom_to_block;
   static unsigned int* macroblock_block_atom_offset;
   static unsigned int* macroblock_block_atoms;
   static unsigned int* macroblock_pair_group_src;
   static unsigned int* macroblock_pair_group_dst;
   static unsigned int* macroblock_pair_group_entry_offset;
   static unsigned int* macroblock_entry_atom_i;
   static unsigned int* macroblock_entry_atom_j;
   static unsigned int* macroblock_entry_ih;
   static unsigned int* macroblock_entry_jslot;
   static unsigned int* macroblock_entry_local_i;
   static unsigned int* macroblock_entry_local_j;

   // Input
   static int* gpu_mode;
   static int* gpu_rng;
   static int* gpu_rng_seed;

   static real* mavg_buff;
   static real* mavg2_buff;
   static real* mavg4_buff;
   static real* eavg_buff;
   static real* eavg2_buff;


   static unsigned int* avrg_step;
   static unsigned int* avrg_buff;
   static unsigned int* cumu_step;
   static unsigned int* cumu_buff;
   static unsigned int* eavrg_step;
   static unsigned int* eavrg_buff;
   static unsigned int* skyno_step;
   static unsigned int* skyno_buff;


   // Initiators
    static void setFlagPointers(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy, 
                                char* p_do_avrg, char* p_do_proj_avrg, char* p_do_cumu,
                                unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
                                unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc,
                                char* p_do_gpu_correlations, char* p_do_gpu_macroblocks);

    static void setConstantPointers(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
                                    unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
                                    real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
                                    real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
                                    unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
                                    unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
                                    unsigned int* p_eavrg_step, unsigned int* p_eavrg_buff, unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
                                    unsigned int* p_skyno_step, unsigned int* p_skyno_buff,  unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
                                    unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep);

    static void setHamiltonianPointers(real* p_ncoup, unsigned int* p_nlist, unsigned int* p_nlistsize,
                                       real* p_dmvect, unsigned int* p_dmlist, unsigned int* p_dmlistsize,
                                       real* p_kaniso, real* p_eaniso, unsigned int* p_taniso, real* p_sb,
                                       real* p_j_tensor, unsigned int* p_aHam, 
                                       real* p_external_field, real* p_btorque, real* p_Temp_array, 
                                       real * p_ipTemp, unsigned int * p_ipmcnstep,
                                       real * p_ipTemp_array, unsigned int* p_ipnstep,
                                       real * p_ipdelta_t, real * p_iplambda1);

    static void setLatticePointers(real* p_beff, real* p_b2eff, real* p_emomM, real* p_emom, real* p_emom2, 
                                   real* p_mmom, real* p_mmom0, real* p_mmom2, real* p_mmomi,
                                   real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list);
   static void setMacroBlockPointers(unsigned int* p_enabled, unsigned int* p_nblocks,
                                       unsigned int* p_npair_groups, unsigned int* p_nentries,
                                       unsigned int* p_max_block_atoms,
                                       unsigned int* p_atom_to_block,
                                      unsigned int* p_block_atom_offset,
                                      unsigned int* p_block_atoms,
                                      unsigned int* p_pair_group_src,
                                      unsigned int* p_pair_group_dst,
                                       unsigned int* p_pair_group_entry_offset,
                                       unsigned int* p_entry_atom_i,
                                       unsigned int* p_entry_atom_j,
                                       unsigned int* p_entry_ih,
                                       unsigned int* p_entry_jslot,
                                       unsigned int* p_entry_local_i,
                                       unsigned int* p_entry_local_j);

    static void setMeasurablePointers(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
                                       real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
                                       real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
                                       real* p_eavg_buff, real* p_eavg2_buff, 
                                       real* p_spinwait, real* p_autocorr_buff, real* p_indxb_ac, 
                                       real* p_traj_step, real* p_traj_buff, real* p_traj_atom,
                                       real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj);
   
   static void setCorrelationPointers(real* p_q, real* p_r_mid, real* p_coord, real* p_w, void* p_m_k, 
                                     void* p_m_kw, void* p_m_kt, real* p_deltat_corr, real* p_scstep_arr, int* p_sc_nsamp, int* p_sc_tidx);

    
   
    /*static void setConstantPointers(char* p1, int* p2, unsigned int* p3, unsigned int* p4, unsigned int* p5,
                                   unsigned int* p6, unsigned int* p7, real* p8, real* p9, real* p10,
                                   real* p11, real* p12, real* p13, real* p14, int* p15, char* p16,
                                   unsigned int* p17, unsigned int* p18, unsigned int* p19,
                                   unsigned int* p20, unsigned int* p21, real * p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* ipnphase);
    */

    /*static void setMatrixPointers(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5, real* p6,
                                 real* p7, real* p8, real* p9, real* p10, real* p11, real* p12, real* p13,
                                 real* p14, real* p15, real* p16, unsigned int* p17, unsigned int* p18,
                                 real* p19, real* p20, real* p21, unsigned int* p22, real* p23, unsigned int* p24, 
                                 real* p_ipTemp, unsigned int * p_ipmcnstep, real* ipTemp_array, unsigned int* ipnstep);
    */

   static void setInputDataPointers(int* p1, int* p2, int* p3);
};
