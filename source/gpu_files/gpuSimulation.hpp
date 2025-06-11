#pragma once

#include <curand.h>
#include "c_headers.hpp"
#include "tensor.hpp"
#include "gpuStructures.hpp"
#include "real_type.h"

class GpuSimulation {
private:
    Flag Flags; //those are constants
    SimulationParameters SimParam; //those are constants

    hostLattice cpuLattice; //those are host matrices
    hostMeasurables cpuMeasurebles;//those are host matrices
    hostHamiltonian cpuHamiltonian;//those are host matrices

    deviceLattice gpuLattice; //those are device matrices
    deviceMeasurables gpuMeasurebles;//those are device matrices
    deviceHamiltonian gpuHamiltonian;//those are device matrices

class GpuSDSimulation {
private:
   bool isInitiatedSD;
   void printMdStatus(std::size_t mstep, GpuSimulation& gpuSim);

public:
   GpuSDSimulation();
   ~GpuSDSimulation();

   void SDmphase(GpuSimulation& gpuSim);
   void SDiphase(GpuSimulation& gpuSim);
};

   // class CudaMCSimulation;  // runs mc simulation

/*class CudaMCSimulation {
private:
   bool isInitiatedSD;
   void printMdStatus(std::size_t mstep, GpuSimulation& gpuSim);
   void printMdStatus_iphase(std::size_t mstep, GpuSimulation& gpuSim, int step);

public:
   CudaMCSimulation();
   ~CudaMCSimulation();

   void MCmphase(GpuSimulation& gpuSim);
   void MCiphase(GpuSimulation& gpuSim);
   void MCiphase_bf(GpuSimulation& gpuSim);
   void MCmphase_bf(GpuSimulation& gpuSim);
};
*/
    bool isInitiated;
    bool isFreed;
    //void printConstants();

    void initiate_fortran_cpu_matrices(); //initiates cpu matrices
    bool gpuHasNoData();

public:
    GpuSimulation();
    ~GpuSimulation();

    void initiateConstants();  // initiates cpuFlags and cpuParameters
    bool initiateMatrices();   // allocates and initiates gpu matrices from cpu matrices using copyFromFortran, first calling initiate_fortran_cpu_matrices();
    void copyFromFortran();    // device to host 
    void copyToFortran();      // host to device
    void release();            // frees gpu matrices

    void gpuRunSimulation(const int whichsim, const int whichphase, const char bf);

};

