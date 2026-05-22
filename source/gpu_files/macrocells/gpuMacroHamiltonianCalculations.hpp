#pragma once
#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include <numeric>
#include "gpu_wrappers.h"

class GpuMacroHamiltonianCalculations{
private:
    bool isallocated;

    const bool do_j_tensor;
    const bool do_dm;
    const int do_aniso;
    const int do_ene;

    // System size
    const std::size_t N;
    const std::size_t M;
    const std::size_t NH;
    const std::size_t mnn;
    const std::size_t mnndm;

    //Jij
    const GpuTensor<unsigned int, 1>& nlistsize;
    const GpuTensor<unsigned int, 2>& nlist;
    const GpuTensor<real, 2>& ncoup;
    //aniso
    const GpuTensor<unsigned int, 1> taniso;
    const GpuTensor<real, 2> eaniso;
    const GpuTensor<real, 2> kaniso;
    const GpuTensor<real, 1> sb;
    //Zeeman field
    const GpuTensor<real, 3> extfield;
    //macrocells 
    const deviceMacrocell& macro;


    unsigned int maxThreads;
    unsigned int numThreads;
    unsigned int numBlocks;
    dim3 threads;
    dim3 blocks;

public:
    // Constructor
    GpuMacroHamiltonianCalculations(const Flag Flags, const SimulationParameters SimParam, const deviceHamiltonian& gpuHamiltonian, const deviceMacrocell& gpuMacro);
    // Destructor
    ~GpuMacroHamiltonianCalculations();
    void calculate(deviceLattice& gpuLattice, deviceEnergies& gpuEnergies, bool measure);
    // Initiator
   // bool initiate(const Flag Flags, const SimulationParameters SimParam, const hostCorrelations& cpuCorrelations);
    // Releaser
   // void release();
    // Measurements
   // void measure(std::size_t mstep);

private:
    //void measure_SC(std::size_t mstep);

};
