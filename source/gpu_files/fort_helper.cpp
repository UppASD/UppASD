//////////////////////////////////////////////////
// Helper functions that calls the C++ routines //
//////////////////////////////////////////////////

#include "gpuSimulation.hpp"
//#include "mdSimulation.hpp"

#include "gpu_wrappers.h"
#include "measurement/memoryMeasurement.h"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

// Report a GPU failure that reached the Fortran boundary, then terminate.
//
// ASSERT_GPU (base.hpp) throws std::runtime_error carrying the driver error
// string for every failed GPU call, so an allocation failure otherwise
// propagates straight through these extern "C" entry points and aborts with a
// bare "terminate called after throwing...". This turns that into a readable
// banner: what failed, the tensor tracker breakdown, and the driver's
// free/total. The exception only carries the error string, so we distinguish
// out-of-memory from a generic GPU error by matching that string.
namespace {

[[noreturn]] void reportGpuFailureAndExit(const char* where, const std::exception& e) {
   std::cout.flush();
   std::fflush(stdout);

   const std::string msg = e.what();
   const bool oom = msg.find("out of memory") != std::string::npos;

   std::fprintf(stderr,
                "\n========================================================================\n");
   if(oom) {
      std::fprintf(stderr, " GPU OUT OF MEMORY in %s\n", where);
   } else {
      std::fprintf(stderr, " GPU ERROR in %s\n", where);
   }
   std::fprintf(stderr, "   %s\n", msg.c_str());
   std::fprintf(stderr,
                "------------------------------------------------------------------------\n");
   std::fflush(stderr);

   // Peak host/device consumption recorded so far by the tensor tracker.
   TensorMemoryTracker::printResults();
   std::cout.flush();

   // Driver's view of what is actually free right now.
   size_t freeBytes = 0, totalBytes = 0;
   const GPU_ERROR_T qerr = GPU_MEM_GET_INFO(&freeBytes, &totalBytes);
   if(qerr == GPU_SUCCESS) {
      std::fprintf(stderr, " GPU memory: %.2f GB free / %.2f GB total\n",
                   static_cast<double>(freeBytes) / 1e9,
                   static_cast<double>(totalBytes) / 1e9);
   } else {
      std::fprintf(stderr, " GPU memory: query failed (%s)\n", GPU_GET_ERROR_STRING(qerr));
   }
   std::fprintf(stderr,
                "========================================================================\n");
   std::fflush(stderr);

   // _Exit, not exit: output is already flushed above, and running the static
   // GpuSimulation destructor here would cudaFree on a context that a fatal (e.g.
   // sticky) GPU error may have torn down, aborting with "driver shutting down".
   std::_Exit(EXIT_FAILURE);
}

}  // namespace

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
   try {
      gpuSim.initiateConstants();
   } catch(const std::exception& e) {
      reportGpuFailureAndExit("gpusim_initiateconstants", e);
   }
}

// is_mc: 1 when the calling Fortran driver is a Monte Carlo phase (mc_driver),
// 0 for spin-dynamics (sd_driver). Lets initiateMatrices skip MC-only buffers
// (eneff, mmom0, mmom2) for SD runs. Each phase re-initiates, so this is per-phase.
void gpusim_initiatematrices_(int *is_mc) {
   try {
      gpuSim.initiateMatrices(*is_mc);
   } catch(const std::exception& e) {
      reportGpuFailureAndExit("gpusim_initiatematrices", e);
   }
}

void gpusim_gpurunsimulation_(int *whichsim, int *whichphase, char* bf){
   try {
      gpuSim.gpuRunSimulation(*whichsim, *whichphase, *bf);
   } catch(const std::exception& e) {
      reportGpuFailureAndExit("gpusim_gpurunsimulation", e);
   }
};

void gpusim_release_(){
   try {
      gpuSim.release();
   } catch(const std::exception& e) {
      reportGpuFailureAndExit("gpusim_release", e);
   }
};


#ifdef __cplusplus
}
#endif
