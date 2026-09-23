
#include "c_headers.hpp"
#include "c_helper.h"
#include "gpuDepondtIntegrator.hpp"

#include "gpuHamiltonianCalculations.hpp"
#include "gpuMeasurement.hpp"
#include "gpuMomentUpdater.hpp"
#include "gpuSimulation.hpp"
#include "gpuStructures.hpp"
#include "gpuStructures.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"
#include "gpuParallelizationHelper.hpp"
#include "measurementFactory.hpp"
#include "correlationFactory.hpp"
#include "measurementQueue.hpp"
#include "cpuRestMeasurement.hpp"

#include "gpu_wrappers.h"
#include "gpuCorrelations.hpp"

using ParallelizationHelper = GpuParallelizationHelper;

GpuSimulation::GpuMSSimulation::GpuMSSimulation() {
   // isInitiatedSD = false;
}

GpuSimulation::GpuMSSimulation::~GpuMSSimulation() {
}

// Printing simulation status
// Added copy to fortran line so that simulation status is printed correctly > Thomas Nystrand 14/09/09
void GpuSimulation::GpuMSSimulation::printMdStatus(std::size_t mstep, GpuSimulation& gpuSim) {

}

// Multiscale initial phase
void GpuSimulation::GpuMSSimulation::MSiphase(GpuSimulation& gpuSim) {
   
}

// Multiscale measurement phase
void GpuSimulation::GpuMSSimulation::MSmphase(GpuSimulation& gpuSim) {
  
}

