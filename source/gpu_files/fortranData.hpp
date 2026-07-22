#pragma once

#include "c_headers.hpp"
#include "real_type.h"
#include <complex>
#include <thrust/complex.h>

// The Fortran ABI always supplies real(dblprec) storage.  Keep this type
// independent from the precision selected for device-side `real`.
using fortran_real = double;
using fortran_complex = std::complex<fortran_real>;

#define real fortran_real

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
   static unsigned int* NA;
   static unsigned int* Natom_full;
   static unsigned int* N1;
   static unsigned int* N2;
   static unsigned int* N3;
   static char* do_gpu_convolution;
   static char* BC1;
   static char* BC2;
   static char* BC3;
   static real* C1;
   static real* C2;
   static real* C3;
   static real* Bas;
   // Macrocell data for the CV6 dipole backend. cell_index is Fortran
   // one-based; conversion happens in the device aggregation kernel.
   static int* do_dip;
   static unsigned int* num_macro;
   static unsigned int* macro_block_x;
   static unsigned int* macro_block_y;
   static unsigned int* macro_block_z;
   static unsigned int* macro_cell_index;
   static unsigned int* macro_nlistsize;
   static real* macro_center;
   static real* macro_min_coord;
   static real* macro_max_coord;
   static unsigned int* pme_num_macro;
   static unsigned int* pme_macro_grid;
   static unsigned int* pme_cell_index;
   static unsigned int* pme_macro_nlistsize;
   static real* pme_macro_center;
   static real* pme_macro_min_coord;
   static real* pme_macro_max_coord;
   

   static unsigned int* nq;
   static unsigned int* sc_step;
   static unsigned int* sc_sep;
   static char* do_sc;
   static unsigned int* sc_max_nstep;
   static unsigned int* sc_window_fun;
   static unsigned int* nw;
   static char* do_sc_proj;
   static char* do_sc_projch;
   static unsigned int* NT;
   static unsigned int* Nchmax;

   static real* r_mid;
   static real* q;
   static real* coord;
   static real* w;
   static int* atype;
   // Atomic types for lattice measurements.  This is distinct from the
   // correlation-grid type array, which may have a different extent.
   static int* lattice_atype;
   static int* achtype;

   static real* deltat_corr;
   static real* scstep_arr;
   static int* sc_nsamp_ptr;  // Pointer to GPU-computed sample count
   static int* sc_tidx_ptr;   // Pointer to GPU-computed time step index
   static fortran_complex * m_k;
   static fortran_complex * m_kw;
   static fortran_complex * m_kt;
   static fortran_complex* m_k_proj;
   static fortran_complex* m_k_projch;
   static fortran_complex* m_kt_proj;
   static fortran_complex* m_kt_projch;
   static fortran_complex* m_kw_proj;
   static fortran_complex* m_kw_projch;


   //delta_t;

   static real* delta_t;
   static real* gamma;
   static real* k_bolt;
   static real* mub;
   static real* mry;
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
   static unsigned int* do_ene;             // Calculate and plot energy (0/1)
   static char* do_skyno;
   static char* do_gpu_correlations;
   static char* do_gpu_timings;
   static char* real_time_measure;
   static char* do_avrg_proj;
   static char* do_avrg_projch;
   static char* do_cumu_proj;


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

   // Inputindxb_ac, autocorr_buff
   static int* gpu_mode;
   static int* gpu_rng;
   static int* gpu_rng_seed;

   static real* mavg_buff;
   static real* mavg2_buff;
   static real* mavg4_buff;
   static real* eavg_buff;
   static real* eavg2_buff;
   static real* autocorr_buff;
   static real* spinwait;
   static real* indxb_ac; //TODO: get rid of once printed on C
   static unsigned int* spinwaittable;
   static int* achem_ch;
   static int* asite_ch;


   static unsigned int* avrg_step;
   static unsigned int* avrg_buff;
   static unsigned int* cumu_step;
   static unsigned int* cumu_buff;
   static unsigned int* ene_step;
   static unsigned int* ene_buff;
   static unsigned int* skyno_step;
   static unsigned int* skyno_buff;
   static unsigned int* ac_step;
   static unsigned int* ac_buff;
   static unsigned int* nspinwait;

   static unsigned int* do_ralloy;


   // Initiators
    static void setFlagPointers(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy, 
                                char* p_do_avrg, char* p_do_proj_avrg, char* p_do_projch_avrg, char* p_do_cumu,
                                char*p_do_cumu_projo, unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
                                unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc,
                                char* p_do_gpu_correlations, char* p_do_gpu_timings, char* p_real_time_measure,
                                char* p_do_sc_proj, char* p_do_sc_projch, unsigned int* p_do_ralloy);

    static void setConstantPointers(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
                                    unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
                                    real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
                                    real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
                                    unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
                                    unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
                                    unsigned int* p_ene_step, unsigned int* p_ene_buff, unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
                                    unsigned int* p_skyno_step, unsigned int* p_skyno_buff,  unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
                                    unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep,
                                    unsigned int* p_nspinwait, unsigned int* p_ac_step, unsigned int* p_ac_buff,
                                    unsigned int* p_nt, unsigned int* p_Nchmax, real* p_mry,
                                    unsigned int* p_NA, unsigned int* p_Natom_full);

    static void setGpuGeometryPointers(unsigned int* p_N1, unsigned int* p_N2, unsigned int* p_N3,
                                       unsigned int* p_NA, char* p_BC1, char* p_BC2, char* p_BC3,
                                       real* p_C1, real* p_C2, real* p_C3, real* p_Bas,
                                       char* p_do_gpu_convolution);

    static void setMacrocellPointers(int* p_do_dip, unsigned int* p_num_macro,
                                     unsigned int* p_block_x, unsigned int* p_block_y,
                                     unsigned int* p_block_z, unsigned int* p_cell_index,
                                     unsigned int* p_macro_nlistsize, real* p_macro_center,
                                     real* p_macro_min_coord, real* p_macro_max_coord);
   static void setPmeMacrocellPointers(unsigned int* p_num_macro, unsigned int* p_macro_grid,
                                       unsigned int* p_cell_index, unsigned int* p_macro_nlistsize,
                                       real* p_macro_center, real* p_macro_min_coord, real* p_macro_max_coord);
   static void clearPmeMacrocellPointers();
    static void clearMacrocellPointers();

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
                                   real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list, int* p_atype);

    static void setMeasurablePointers(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
                                       real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
                                       real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
                                       real* p_eavg_buff, real* p_eavg2_buff, 
                                       int* p_traj_step, int* p_traj_buff, int* p_traj_atom,
                                       real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj,
                                       unsigned int* p_spinwaitt, real* p_spinwait, real* p_indxb_ac, real* p_autocorr_buff,
                                       int* p_achem_ch, int* p_asite_ch);
   
   static void setCorrelationPointers(real* p_q, real* p_r_mid, real* p_coord, real* p_w, void* p_m_k, 
                                       void* p_m_kw, void* p_m_kt, real* p_deltat_corr, real* p_scstep_arr, 
                                       int* p_sc_nsamp, int* p_sc_tidx, int* p_atype, int* p_achtype, void* p_m_k_proj, 
                                       void* p_m_k_projch, void* p_m_kt_proj, void* p_m_kt_projch, void* p_m_kw_proj, 
                                       void* p_m_kw_projch);

    
   
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

#undef real
