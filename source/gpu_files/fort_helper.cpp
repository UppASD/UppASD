//////////////////////////////////////////////////
// Helper functions that calls the C++ routines //
//////////////////////////////////////////////////

#include "cudaMdSimulation.hpp"
#include "mdSimulation.hpp"


#ifdef __cplusplus
extern "C" {
#endif


// C++
static MdSimulation cMdSim;


void cmdsim_initiateconstants_() {
   cMdSim.initiateConstants();
}


void cmdsim_initiatefortran_() {
   cMdSim.initiateFortran();
}


void cmdsim_initiateown_() {
   cMdSim.initiateOwn();
}


void cmdsim_measurementphase_() {
   cMdSim.measurementPhase();
}


void cmdsim_readmatrices_() {
   cMdSim.copyFromFortran();
}


void cmdsim_writematrices_() {
   cMdSim.copyToFortran();
}


// Cuda
static CudaMdSimulation cudaMdSim;


void cudamdsim_initiateconstants_() {
   cudaMdSim.initiateConstants();
}


void cudamdsim_initiatematrices_() {
   cudaMdSim.initiateMatrices();
}


void cudamdsim_measurementphase_() {
   cudaMdSim.measurementPhase();
}


#ifdef __cplusplus
}
#endif

