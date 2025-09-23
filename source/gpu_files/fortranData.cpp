#include "fortranData.hpp"

#include "real_type.h"

// Constants
char* FortranData::stt;
int* FortranData::SDEalgh;

unsigned int* FortranData::rstep;
unsigned int* FortranData::nstep;
unsigned int* FortranData::Natom;
unsigned int* FortranData::nHam;
unsigned int* FortranData::Mensemble;
unsigned int* FortranData::max_no_neigh;
unsigned int* FortranData::ipmcnphase;
unsigned int* FortranData::mcnstep;
unsigned int* FortranData::ipnphase;

unsigned int* FortranData::nq;
unsigned int* FortranData::sc_step;
unsigned int* FortranData::sc_sep;
char* FortranData::do_sc;
unsigned int* FortranData::sc_max_nstep;
unsigned int* FortranData::sc_window_fun;
unsigned int* FortranData::nw;

real* FortranData::r_mid;
real* FortranData::q;
real* FortranData::coord;
real* FortranData::w;
std::complex<real>* FortranData::m_k;
std::complex<real>* FortranData::m_kw;
std::complex<real>* FortranData::m_kt;

real* FortranData::delta_t;
real* FortranData::gamma;
real* FortranData::k_bolt;
real* FortranData::mub;
real* FortranData::damping;
real* FortranData::Temp;

real* FortranData::binderc;
real* FortranData::mavg;

int* FortranData::mompar;
char* FortranData::initexc;

unsigned int* FortranData::do_dm;
unsigned int* FortranData::do_jtensor;
unsigned int* FortranData::do_aniso;
unsigned int* FortranData::max_no_dmneigh;

char* FortranData::do_cuda_measurements;
char* FortranData::do_avrg;
char* FortranData::do_cumu;
char* FortranData::do_autocorr;
unsigned int* FortranData::plotenergy;
char* FortranData::do_skyno;

// Matrices
unsigned int * FortranData::aHam;
real* FortranData::ncoup;
unsigned int* FortranData::nlist;
unsigned int* FortranData::nlistsize;
real* FortranData::dmvect;
unsigned int* FortranData::dmlist;
unsigned int* FortranData::dmlistsize;
real* FortranData::beff;
real* FortranData::b2eff;
real* FortranData::emomM;
real* FortranData::emom;
real* FortranData::emom2;
real* FortranData::external_field;
real* FortranData::mmom;
real* FortranData::btorque;
real* FortranData::temperature;
real* FortranData::ipTemp;
unsigned int* FortranData::ipmcnstep;
real* FortranData::ipTemp_array;
unsigned int* FortranData::ipnstep;
real* FortranData::mmom0;
real* FortranData::mmom2;
real* FortranData::mmomi;

real* FortranData::j_tensor;
real* FortranData::kaniso;
real* FortranData::eaniso;
unsigned int* FortranData::taniso;
real* FortranData::sb;

real* FortranData::dxyz_vec;
int* FortranData::dxyz_atom;
int* FortranData::dxyz_list;

// GPU stuff
int* FortranData::gpu_mode;
int* FortranData::gpu_rng;
int* FortranData::gpu_rng_seed;

real* FortranData::mavg_buff;
real* FortranData::mavg2_buff;
real* FortranData::mavg4_buff;
real* FortranData::eavg_buff;
real* FortranData::eavg2_buff;


unsigned int* FortranData::avrg_step;
unsigned int* FortranData::avrg_buff;
unsigned int* FortranData::cumu_step;
unsigned int* FortranData::cumu_buff;
unsigned int* FortranData::eavrg_step;
unsigned int* FortranData::eavrg_buff;
unsigned int* FortranData::skyno_step;
unsigned int* FortranData::skyno_buff;

void FortranData::setFlagPointers(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy,
                                  char* p_do_avrg, char* p_do_proj_avrg, char* p_do_cumu,
                                  unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
                                  unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc){


   do_dm = p_do_dm;
   do_jtensor = p_do_jtensor;
   do_aniso = p_do_anisotropy;
   do_cuda_measurements = p_do_cuda_measurements;
   do_avrg = p_do_avrg;
   do_cumu = p_do_cumu;
   do_autocorr = p_do_autocorr;
   plotenergy = p_plotenergy;
   do_skyno = p_do_skyno;
   do_sc = p_do_sc;
}

void FortranData::setConstantPointers(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
                                      unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
                                      real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
                                      real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
                                      unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
                                      unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
                                      unsigned int* p_eavrg_step, unsigned int* p_eavrg_buff,  unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
                                      unsigned int* p_skyno_step, unsigned int* p_skyno_buff, unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
                                      unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep){

   stt = p_stt;
   SDEalgh = p_SDEalgh;
                                    
   rstep = p_rstep;
   nstep = p_nstep;
   Natom = p_Natom;
   nHam = p_nHam;
   Mensemble = p_Mensemble;
   max_no_neigh = p_max_no_neigh;
                                    
   delta_t = p_delta_t;
   gamma = p_gamma;
   k_bolt = p_k_bolt;
   mub = p_mub;
   damping = p_mplambda1;
                                    
   binderc = p_binderc;
   mavg = p_mavg;
                                    
   mompar = p_mompar;
   initexc = p_initexc;
                                    
                                      
   max_no_dmneigh = p_max_no_dmneigh;

   Temp = p_Temp;
   ipnphase = p_ipnphase;
   ipmcnphase = p_ipmcnphase;
   mcnstep = p_mcnstep;

   avrg_step = p_avrg_step;
   avrg_buff = p_avrg_buff;
   cumu_step = p_cumu_step;
   cumu_buff = p_cumu_buff;
   eavrg_step = p_eavrg_step;
   eavrg_buff = p_eavrg_buff;

   skyno_step = p_skyno_step;
   skyno_buff = p_skyno_buff;
   nq = p_nq;
   sc_window_fun = p_sc_window_fun;
   nw = p_nw;
   sc_sep = p_sc_sep;
   sc_step = p_sc_step; 
   sc_max_nstep = p_sc_max_nstep;
}

void FortranData::setHamiltonianPointers(real* p_ncoup, unsigned int* p_nlist, unsigned int* p_nlistsize,
                                         real* p_dmvect, unsigned int* p_dmlist, unsigned int* p_dmlistsize,
                                         real* p_kaniso, real* p_eaniso, unsigned int* p_taniso, real* p_sb,
                                         real* p_j_tensor, unsigned int* p_aHam, 
                                         real* p_external_field, real* p_btorque, real* p_Temp_array, 
                                         real * p_ipTemp, unsigned int * p_ipmcnstep,
                                         real * p_ipTemp_array, unsigned int* p_ipnstep){

   ncoup = p_ncoup;
   nlist = p_nlist;
   nlistsize = p_nlistsize;
   external_field = p_external_field;
   btorque = p_btorque;
   temperature = p_Temp_array;
   dmvect = p_dmvect;
   dmlist = p_dmlist;
   dmlistsize = p_dmlistsize;
   j_tensor = p_j_tensor;
   kaniso = p_kaniso;
   eaniso = p_eaniso;
   taniso = p_taniso;
   sb = p_sb;
   aHam = p_aHam;
   ipTemp = p_ipTemp;
   ipmcnstep = p_ipmcnstep;
   ipTemp_array = p_ipTemp_array;
   ipnstep = p_ipnstep;
}


void FortranData::setLatticePointers(real* p_beff, real* p_b2eff, real* p_emomM, real* p_emom, real* p_emom2, 
                                     real* p_mmom, real* p_mmom0, real* p_mmom2, real* p_mmomi,
                                     real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list){


   beff = p_beff;
   b2eff = p_b2eff;
   emomM = p_emomM;
   emom = p_emom;
   emom2 = p_emom2;
   mmom = p_mmom;
   mmom0 = p_mmom0;
   mmom2 = p_mmom2;
   mmomi = p_mmomi;

   dxyz_vec = p_dxyz_vec;
   dxyz_atom = p_dxyz_atom;
   dxyz_list = p_dxyz_list;
}

//TODO:binderc, autocorr_buff, spinwait
void FortranData::setMeasurablePointers(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
                                         real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
                                         real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
                                         real* p_eavg_buff, real* p_eavg2_buff, 
                                         real* p_spinwait, real* p_autocorr_buff, real* p_indxb_ac, 
                                         real* p_traj_step, real* p_traj_buff, real* p_traj_atom,
                                         real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj){

   mavg_buff = p_mavg_buff;
   mavg2_buff = p_mavg2_buff;
   mavg4_buff = p_mavg4_buff;
   eavg_buff = p_eavg_buff;
   eavg2_buff = p_eavg2_buff;


}

void FortranData::setCorrelationPointers(real* p_q, real* p_r_mid, real* p_coord, real* p_w,  std::complex<real>* p_m_k, 
                                        std::complex<real>* p_m_kw, std::complex<real>* p_m_kt){


   q = p_q;
   r_mid = p_r_mid;
   coord = p_coord;
   w = p_w;
   m_k = p_m_k;
   m_kw = p_m_kw;
   m_kt = p_m_kt;

}
/*void FortranData::setConstantPointers(char* p1, int* p2, unsigned int* p3, unsigned int* p4, unsigned int* p5,
                                      unsigned int* p6, unsigned int* p7, real* p8, real* p9, real* p10,
                                      real* p11, real* p12, real* p13, real* p14, int* p15, char* p16,
                                      unsigned int* p17, unsigned int* p18, unsigned int* p19,
                                      unsigned int* p20, unsigned int* p21, real * p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep,
                                      unsigned int * p_ipnphase) {
   stt = p1;
   SDEalgh = p2;

   rstep = p3;
   nstep = p4;
   Natom = p5;
   nHam = p21;
   Mensemble = p6;
   max_no_neigh = p7;

   delta_t = p8;
   gamma = p9;
   k_bolt = p10;
   mub = p11;
   damping = p12;

   binderc = p13;
   mavg = p14;

   mompar = p15;
   initexc = p16;

   do_dm = p17;
   max_no_dmneigh = p18;
   do_jtensor = p19;
   do_aniso = p20;
   Temp = p_Temp;
   ipnphase = p_ipnphase;
   ipmcnphase = p_ipmcnphase;
   mcnstep = p_mcnstep;
}*/

/*void FortranData::setMatrixPointers(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                    real* p6, real* p7, real* p8, real* p9, real* p10, real* p11, real* p12,
                                    real* p13, real* p14, real* p15, real* p16, unsigned int* p17,
                                    unsigned int* p18, real* p19, real* p20, real* p21, unsigned int* p22,
                                    real* p23, unsigned int* p24, real * p_ipTemp, unsigned int * p_ipmcnstep,
                                    real * p_ipTemp_array, unsigned int* p_ipnstep) {
   ncoup = p1;
   nlist = p2;
   nlistsize = p3;
   beff = p4;
   b2eff = p5;
   emomM = p6;
   emom = p7;
   emom2 = p8;
   external_field = p9;
   mmom = p10;
   btorque = p11;
   temperature = p12;
   mmom0 = p13;
   mmom2 = p14;
   mmomi = p15;
   dmvect = p16;
   dmlist = p17;
   dmlistsize = p18;
   j_tensor = p19;
   kaniso = p20;
   eaniso = p21;
   taniso = p22;
   sb = p23;
   aHam = p24;
   ipTemp = p_ipTemp;
   ipmcnstep = p_ipmcnstep;
   ipTemp_array = p_ipTemp_array;
   ipnstep = p_ipnstep;
}*/

void FortranData::setInputDataPointers(int* p1, int* p2, int* p3) {
   gpu_mode = p1;
   gpu_rng = p2;
   gpu_rng_seed = p3;
}

// Fortran helpers
extern "C" void fortrandata_setflags_(unsigned int* p_do_dm, unsigned int* p_do_jtensor, unsigned int* p_do_anisotropy, 
   char* p_do_avrg, char* p_do_proj_avrg, char* p_do_cumu,
   unsigned int* p_plotenergy, char* p_do_autocorr, char* p_do_tottraj,
   unsigned int* p_ntraj, char* p_do_cuda_measurements, char* p_do_skyno, char* p_do_sc) {
FortranData::setFlagPointers(
   p_do_dm, p_do_jtensor, p_do_anisotropy, p_do_avrg, p_do_proj_avrg, p_do_cumu,  p_plotenergy, 
   p_do_autocorr, p_do_tottraj, p_ntraj, p_do_cuda_measurements, p_do_skyno, p_do_sc);
}

extern "C" void fortrandata_setconstants_(char* p_stt, int* p_SDEalgh, unsigned int* p_rstep, unsigned int* p_nstep,
   unsigned int* p_Natom, unsigned int* p_Mensemble, unsigned int* p_max_no_neigh, 
   real* p_delta_t, real* p_gamma, real* p_k_bolt, real* p_mub, real* p_mplambda1,
   real* p_binderc, real* p_mavg, int* p_mompar, char* p_initexc, unsigned int* p_max_no_dmneigh,
   unsigned int* p_nHam, real* p_Temp, unsigned int* p_ipmcnphase, unsigned int* p_mcnstep, unsigned int* p_ipnphase,
   unsigned int* p_avrg_step, unsigned int* p_avrg_buff, unsigned int* p_cumu_step, unsigned int* p_cumu_buff,
   unsigned int* p_eavrg_step, unsigned int* p_eavrg_buff, unsigned int*p_tottraj_step, unsigned int*p_tottraj_buff,
   unsigned int* p_skyno_step, unsigned int* p_skyno_buff,  unsigned int* p_nq, unsigned int* p_sc_window_fun, unsigned int* p_nw,
   unsigned int* p_sc_sep, unsigned int* p_sc_step, unsigned int* p_sc_max_nstep) {
FortranData::setConstantPointers(
   p_stt, p_SDEalgh, p_rstep, p_nstep, p_Natom, p_Mensemble, p_max_no_neigh, p_delta_t, p_gamma, 
   p_k_bolt, p_mub, p_mplambda1, p_binderc, p_mavg, p_mompar, p_initexc, p_max_no_dmneigh, p_nHam, 
   p_Temp, p_ipmcnphase, p_mcnstep, p_ipnphase,
   p_avrg_step, p_avrg_buff, p_cumu_step, p_cumu_buff, p_eavrg_step, p_eavrg_buff, p_tottraj_step, p_tottraj_buff,
   p_skyno_step, p_skyno_buff, p_nq, p_sc_window_fun, p_nw, p_sc_sep, p_sc_step, p_sc_max_nstep);
}

extern "C" void fortrandata_sethamiltonian_(real* p_ncoup, unsigned int* p_nlist, unsigned int* p_nlistsize,
   real* p_dmvect, unsigned int* p_dmlist, unsigned int* p_dmlistsize,
   real* p_kaniso, real* p_eaniso, unsigned int* p_taniso, real* p_sb,
   real* p_j_tensor, unsigned int* p_aHam, 
   real* p_external_field, real* p_btorque, real* p_Temp_array, 
   real * p_ipTemp, unsigned int * p_ipmcnstep,
   real * p_ipTemp_array, unsigned int* p_ipnstep) {
FortranData::setHamiltonianPointers(
   p_ncoup, p_nlist, p_nlistsize, p_dmvect, p_dmlist, p_dmlistsize,  p_kaniso, p_eaniso, p_taniso, p_sb, 
   p_j_tensor, p_aHam, p_external_field, p_btorque, p_Temp_array, 
   p_ipTemp, p_ipmcnstep, p_ipTemp_array, p_ipnstep);
}

extern "C" void fortrandata_setlattice_(real* p_beff, real* p_b2eff, real* p_emomM, real* p_emom, real* p_emom2, 
   real* p_mmom, real* p_mmom0, real* p_mmom2, real* p_mmomi, real* p_dxyz_vec, int* p_dxyz_atom, int* p_dxyz_list) {
FortranData::setLatticePointers(
   p_beff, p_b2eff, p_emomM, p_emom, p_emom2, p_mmom, p_mmom0, p_mmom2, p_mmomi, p_dxyz_vec, p_dxyz_atom, p_dxyz_list);
}

extern "C" void fortrandata_setmeasurables_(real* p_mavg_buff, real* p_mavg2_buff, real* p_mavg4_buff,
   real* p_mavg_buff_proj, real* p_mavg2_buff_proj, real* p_mavg4_buff_proj, 
   real* p_binderc, real* p_avrgmcum, real* p_avrgm2cum, real* p_avrgm4cum, 
   real* p_eavg_buff, real* p_eavg2_buff, 
   real* p_spinwait, real* p_autocorr_buff, real* p_indxb_ac, 
   real* p_traj_step, real* p_traj_buff, real* p_traj_atom,
   real* p_mmomb, real* p_mmomb_traj, real* p_emomb, real* p_emomb_traj) {
FortranData::setMeasurablePointers(
   p_mavg_buff, p_mavg2_buff, p_mavg4_buff, p_mavg_buff_proj, p_mavg2_buff_proj, p_mavg4_buff_proj, 
   p_binderc, p_avrgmcum, p_avrgm2cum, p_avrgm4cum, p_eavg_buff, p_eavg2_buff, 
   p_spinwait, p_autocorr_buff, p_indxb_ac, p_traj_step, p_traj_buff, p_traj_atom,p_mmomb, p_mmomb_traj, p_emomb, p_emomb_traj);
}

extern "C" void fortrandata_setcorrelations_(real* p_q, real* p_r_mid, real* p_coord, real* p_w, std::complex<real>* p_m_k, 
                                             std::complex<real>* p_m_kw, std::complex<real>* p_m_kt) {
FortranData::setCorrelationPointers(
   p_q, p_r_mid, p_coord, p_w,  p_m_k, p_m_kw, p_m_kt);
}
/*extern "C" void fortrandata_setconstants_(char* p1, int* p2, unsigned int* p3, unsigned int* p4,
                                          unsigned int* p5, unsigned int* p6, unsigned int* p7, real* p8,
                                          real* p9, real* p10, real* p11, real* p12, real* p13, real* p14,
                                          int* p15, char* p16, unsigned int* p17, unsigned int* p18,
                                          unsigned int* p19, unsigned int* p20, unsigned int* p21, real* p_Temp, unsigned int* p_ipmcnphase,
                                          unsigned int* p_mcnstep, unsigned int* p_ipnphase) {
   FortranData::setConstantPointers(
       p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p_Temp, p_ipmcnphase, p_mcnstep, p_ipnphase);
}*/

/*extern "C" void fortrandata_setmatrices_(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                         real* p6, real* p7, real* p8, real* p9, real* p10, real* p11,
                                         real* p12, real* p13, real* p14, real* p15, real* p16,
                                         unsigned int* p17, unsigned int* p18, real* p19, real* p20,
                                         real* p21, unsigned int* p22, real* p23, unsigned int* p24, 
                                         real* p_ipTemp, unsigned int* p_ipmcnstep, real* p_ipTemp_array, unsigned int* p_ipnstep) {
   FortranData::setMatrixPointers(p1,
                                  p2,
                                  p3,
                                  p4,
                                  p5,
                                  p6,
                                  p7,
                                  p8,
                                  p9,
                                  p10,
                                  p11,
                                  p12,
                                  p13,
                                  p14,
                                  p15,
                                  p16,
                                  p17,
                                  p18,
                                  p19,
                                  p20,
                                  p21,
                                  p22,
                                  p23,
                                  p24,
                                  p_ipTemp,
                                  p_ipmcnstep,
                                  p_ipTemp_array,
                                  p_ipnstep);
}*/

extern "C" void fortrandata_setinputdata_(int* p1, int* p2, int* p3) {
   FortranData::setInputDataPointers(p1, p2, p3);
}

