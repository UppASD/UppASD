//////////////////////////////////////////////////
// Helper functions that calls the C++ routines //
//////////////////////////////////////////////////

#include "cudaSimulation.hpp"
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
static CudaSimulation cudaSim;

void cudasim_initiateconstants_() {
   cudaSim.initiateConstants();
}

void cudasim_initiatematrices_() {
   cudaSim.initiateMatrices();
}

void cudasim_cudarunsimulation_(int *whichsim, int *whichphase, char* bf){
   cudaSim.cudaRunSimulation(*whichsim, *whichphase, *bf);
};

void cudasim_release_(){
   cudaSim.release();
};


#ifdef __cplusplus
}
#endif

