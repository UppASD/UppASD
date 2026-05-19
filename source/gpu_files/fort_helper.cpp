//////////////////////////////////////////////////
// Helper functions that calls the C++ routines //
//////////////////////////////////////////////////

#include "gpuSimulation.hpp"
//#include "mdSimulation.hpp"

#ifdef __cplusplus
extern "C" {
#endif

// C++
/*static MdSimulation cMdSim;

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
*/
// Cuda
static GpuSimulation gpuSim;

void gpusim_initiateconstants_() {
   gpuSim.initiateConstants();
}

void gpusim_initiatematrices_() {
   gpuSim.initiateMatrices();
}

void gpusim_gpurunsimulation_(int *whichsim, int *whichphase, char* bf){
   gpuSim.gpuRunSimulation(*whichsim, *whichphase, *bf);
};

void gpusim_release_(){
   gpuSim.release();
};


#ifdef __cplusplus
}
#endif

