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

real* FortranData::delta_t;
real* FortranData::gamma;
real* FortranData::k_bolt;
real* FortranData::mub;
real* FortranData::damping;

real* FortranData::binderc;
real* FortranData::mavg;

int* FortranData::mompar;
char* FortranData::initexc;

unsigned int* FortranData::do_dm;
unsigned int* FortranData::do_jtensor;
unsigned int* FortranData::do_aniso;
unsigned int* FortranData::max_no_dmneigh;

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
real* FortranData::mmom0;
real* FortranData::mmom2;
real* FortranData::mmomi;

real* FortranData::j_tensor;
real* FortranData::kaniso;
real* FortranData::eaniso;
unsigned int* FortranData::taniso;
real* FortranData::sb;

// GPU stuff
int* FortranData::gpu_mode;
int* FortranData::gpu_rng;
int* FortranData::gpu_rng_seed;

void FortranData::setConstantPointers(char* p1, int* p2, unsigned int* p3, unsigned int* p4, unsigned int* p5,
                                      unsigned int* p6, unsigned int* p7, real* p8, real* p9, real* p10,
                                      real* p11, real* p12, real* p13, real* p14, int* p15, char* p16,
                                      unsigned int* p17, unsigned int* p18, unsigned int* p19,
                                      unsigned int* p20, unsigned int* p21) {
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
}

void FortranData::setMatrixPointers(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                    real* p6, real* p7, real* p8, real* p9, real* p10, real* p11, real* p12,
                                    real* p13, real* p14, real* p15, real* p16, unsigned int* p17,
                                    unsigned int* p18, real* p19, real* p20, real* p21, unsigned int* p22,
                                    real* p23, unsigned int* p24) {
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
}

void FortranData::setInputDataPointers(int* p1, int* p2, int* p3) {
   gpu_mode = p1;
   gpu_rng = p2;
   gpu_rng_seed = p3;
}

// Fortran helpers
extern "C" void fortrandata_setconstants_(char* p1, int* p2, unsigned int* p3, unsigned int* p4,
                                          unsigned int* p5, unsigned int* p6, unsigned int* p7, real* p8,
                                          real* p9, real* p10, real* p11, real* p12, real* p13, real* p14,
                                          int* p15, char* p16, unsigned int* p17, unsigned int* p18,
                                          unsigned int* p19, unsigned int* p20, unsigned int* p21) {
   FortranData::setConstantPointers(
       p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21);
}

extern "C" void fortrandata_setmatrices_(real* p1, unsigned int* p2, unsigned int* p3, real* p4, real* p5,
                                         real* p6, real* p7, real* p8, real* p9, real* p10, real* p11,
                                         real* p12, real* p13, real* p14, real* p15, real* p16,
                                         unsigned int* p17, unsigned int* p18, real* p19, real* p20,
                                         real* p21, unsigned int* p22, real* p23, unsigned int* p24) {
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
                                  p24);
}

extern "C" void fortrandata_setinputdata_(int* p1, int* p2, int* p3) {
   FortranData::setInputDataPointers(p1, p2, p3);
}

