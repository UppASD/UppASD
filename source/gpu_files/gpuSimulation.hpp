#pragma once

#if defined (CUDA_V)
#include <curand.h>
#elif defined (HIP_V)
 #include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "tensor.hpp"
#include "gpuStructures.hpp"
#include "gpuAdaptiveRuntime.hpp"
#include "real_type.h"

class GpuHamiltonianCalculations;

class GpuSimulation {
private:
    Flag Flags; //those are constants
    SimulationParameters SimParam; //those are constants

    hostLattice cpuLattice; //those are host matrices
    hostMeasurables cpuMeasurebles;//those are host matrices
    hostHamiltonian cpuHamiltonian;//those are host matrices
    hostCorrelations cpuCorrelations;//those are host matrices
    // Precision-converted geometry retained for the lifetime of the GPU run.
    Tensor<real, 1> geometryC1, geometryC2, geometryC3, geometryBas;

    deviceLattice gpuLattice; //those are device matrices
    deviceMeasurables gpuMeasurebles;//those are device matrices
    deviceHamiltonian gpuHamiltonian;//those are device matrices
    deviceEnergies gpuEnergies;
    GpuAdaptiveRuntime gpuAdaptiveRuntime;
    bool adaptiveMaskEnabled = false;
    unsigned int adaptiveUpdateInterval = 1;
    GpuAdaptiveSelectorPolicy adaptiveSelectorPolicy{};
    GpuAdaptiveReconstructionPolicy adaptiveReconstructionPolicy{};
    real adaptivePolarizationThreshold = real(0.9);
    int adaptiveDiagnostics = 0;
    double adaptiveEnergyJumpLimitJ = 0.0;

    const unsigned int maxThreads;
    const unsigned int maxBlocks;

class GpuSDSimulation {
private:
   bool isInitiatedSD;   GpuTensor<real, 1> exchangeM;

   void printMdStatus(std::size_t mstep, GpuSimulation& gpuSim);

public:
   GpuSDSimulation();
   ~GpuSDSimulation();

   void SDmphase(GpuSimulation& gpuSim);
   void SDiphase(GpuSimulation& gpuSim);
};

   // class GpuMCSimulation;  // runs mc simulation

class GpuMCSimulation {
private:
   bool isInitiatedSD;
   void printMdStatus(std::size_t mstep, GpuSimulation& gpuSim);
   void printMdStatus_iphase(std::size_t mstep, GpuSimulation& gpuSim, int step);

public:
   GpuMCSimulation();
   ~GpuMCSimulation();

   void MCmphase(GpuSimulation& gpuSim);
   void MCiphase(GpuSimulation& gpuSim);
   void MCiphase_bf(GpuSimulation& gpuSim);
   void MCmphase_bf(GpuSimulation& gpuSim);
};

    bool isInitiated;
    bool isFreed;
    std::size_t estimatedDeviceBytes = 0;  // M1: projected device footprint, set in initiateMatrices
    bool runIsMC = false;       // M2: this phase is Monte Carlo (needs eneff/mmom0/mmom2)
    bool eneffAliased = false;  // M2: eneff points at beff (SD, non-aniso); do not Free it
    //void printConstants();

    void initiate_fortran_cpu_matrices(); //initiates cpu matrices
    bool gpuHasNoData();

public:
    GpuSimulation();
    ~GpuSimulation();

    void initiateConstants();  // initiates cpuFlags and cpuParameters
    bool initiateMatrices(int is_mc);   // allocates and initiates gpu matrices from cpu matrices using copyFromFortran, first calling initiate_fortran_cpu_matrices();
    void copyFromFortran();    // device to host 
    void copyToFortran();      // host to device
    void release();            // frees gpu matrices
    void updateAdaptiveBlockState(const int* blockState, std::size_t count);
    bool adaptiveEnabled() const { return gpuAdaptiveRuntime.ready(); }
    void advanceAdaptiveStep(std::size_t step,
                             GpuHamiltonianCalculations* hamiltonian = nullptr);

    void gpuRunSimulation(const int whichsim, const int whichphase, const char bf);

};

