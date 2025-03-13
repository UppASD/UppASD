#pragma once

#include <curand.h>
#include "c_headers.hpp"
#include "tensor.cuh"
#include "cudaStructures.hpp"
#include "real_type.h"

class CudaSimulation {
private:
    Flag Flags; //those are constants
    SimulationParameters SimParam; //those are constants

    hostLattice cpuLattice; //those are host matrices
    hostMeasurables cpuMeasurebles;//those are host matrices
    hostHamiltonian cpuHamiltonian;//those are host matrices

    cudaLattice gpuLattice; //those are device matrices
    cudaMeasurables gpuMeasurebles;//those are device matrices
    cudaHamiltonian gpuHamiltonian;//those are device matrices

class CudaSDSimulation {
private:
   bool isInitiatedSD;
   void printMdStatus(std::size_t mstep, CudaSimulation& cudaSim);

public:
   CudaSDSimulation();
   ~CudaSDSimulation();

   void SDmphase(CudaSimulation& cudaSim);
   void SDiphase(CudaSimulation& cudaSim);
};

   // class CudaMCSimulation;  // runs mc simulation

class CudaMCSimulation {
private:
   bool isInitiatedSD;
   void printMdStatus(std::size_t mstep, CudaSimulation& cudaSim);
   void printMdStatus_iphase(std::size_t mstep, CudaSimulation& cudaSim, int step);

public:
   CudaMCSimulation();
   ~CudaMCSimulation();

   void MCmphase(CudaSimulation& cudaSim);
   void MCiphase(CudaSimulation& cudaSim);
   void MCiphase_bf(CudaSimulation& cudaSim);
   void MCmphase_bf(CudaSimulation& cudaSim);
};

    bool isInitiated;
    bool isFreed;
    //void printConstants();

    void initiate_fortran_cpu_matrices(); //initiates cpu matrices
    bool gpuHasNoData();

public:
    CudaSimulation();
    ~CudaSimulation();

    void initiateConstants();  // initiates cpuFlags and cpuParameters
    bool initiateMatrices();   // allocates and initiates gpu matrices from cpu matrices using copyFromFortran, first calling initiate_fortran_cpu_matrices();
    void copyFromFortran();    // device to host 
    void copyToFortran();      // host to device
    void release();            // frees gpu matrices

    void cudaRunSimulation(const int whichsim, const int whichphase, const char bf);

};

