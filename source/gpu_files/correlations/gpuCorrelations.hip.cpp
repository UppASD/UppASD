#pragma once

#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpuCorrelations.hpp"
#include <numeric>
//#include <cooperative_groups.h>
//#include <cooperative_groups/reduce.h>
//#include <thrust/complex.h>
#include <hiprand/hiprand.h>

//namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif



// Constructor
GpuCorrelations::GpuCorrelations() {
    isallocated = 0; 
    printf("\n\n\nSorry, HIP version of correlations does not exist yet, please set do_gpu_correlations to 0 on AMD cards\n\n\n");

}
// Destructor
GpuCorrelations::~GpuCorrelations() {
    release();
}
// Initiator
bool GpuCorrelations::initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations) {
    // Assert that we're not already initialized
    release();
    return true;
}
void GpuCorrelations::release() {
    if (isallocated) {

}

void GpuCorrelations::runCorrelation_spin(const deviceLattice& gpuLattice, const int curstep) {

}

void GpuCorrelations::finalCorrelation(hostCorrelations& cpuCorrelations) {

}





