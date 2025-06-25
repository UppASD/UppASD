#pragma once
#include <cuda.h>
#include <curand.h>

#include "c_headers.hpp"
#include "c_helper.h"
#include "cudaGPUErrchk.hpp"
#include "cudaHamiltonianCalculations.hpp"
#include "cudaMetropolis.cuh"
#include "cudaMetropolis_bruteforce.cuh"
#include "cudaParallelizationHelper.hpp"
#include "cudaSimulation.hpp"
#include "cudaStructures.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.cuh"
#include "measurementFactory.cuh"

CudaSimulation::CudaMCSimulation::CudaMCSimulation() {
   // isInitiatedSD = false;
}

CudaSimulation::CudaMCSimulation::~CudaMCSimulation() {
}

void CudaSimulation::CudaMCSimulation::printMdStatus(std::size_t mstep, CudaSimulation& cudaSim) {
   if(cudaSim.SimParam.mcnstep > 20) {
      if(mstep % ((cudaSim.SimParam.rstep + cudaSim.SimParam.mcnstep) / 20) == 0) {
         cudaSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
         std::printf("CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (cudaSim.SimParam.rstep + cudaSim.SimParam.mcnstep),
                     *cudaSim.SimParam.mavg,
                     *cudaSim.SimParam.binderc);
      }
   } else {
      cudaSim.copyToFortran();
      fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
      std::printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *cudaSim.SimParam.mavg);
   }
}

void CudaSimulation::CudaMCSimulation::printMdStatus_iphase(std::size_t mstep, CudaSimulation& cudaSim, int step) {
   if(step > 20) {
      if(mstep % ((step) / 20) == 0) {
         //cudaMC.mom_update(cudaSim.gpuLattice);

         cudaSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
         std::printf("CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (step),
                     *cudaSim.SimParam.mavg,
                     *cudaSim.SimParam.binderc);
      }
   } else {
      //cudaMC.mom_update(CudaSim.gpuLattice);
      cudaSim.copyToFortran();
      fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
      std::printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *cudaSim.SimParam.mavg);
   }
}

// Monte Carlo initial phase (strict)
void CudaSimulation::CudaMCSimulation::MCiphase(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC initial phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // Metropolis
   CudaMetropolis cudaMC;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;
   //printf("HERE - 0\n");

   // Initiate MC and Hamiltonian
   num_subL = cudaMC.initiate(cudaSim.SimParam, cudaSim.cpuHamiltonian, cudaSim.cpuLattice);
   //printf("HERE - 1\n");
   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = cudaSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);
   //printf("HERE - 2\n");

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = cudaSim.cpuLattice.ipmcnstep(it);
   beta = 1/(cudaSim.cpuLattice.ipTemp(it) * cudaSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(cudaSim.gpuLattice);
   stopwatch.add("hamiltonian");
   //printf("HERE - 3\n");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, cudaSim, mcs);
      //printMdStatus(mstep, cudaSim);
         // Perform Metropolis sweep
         for(std::size_t sub = 0; sub < num_subL; sub++){
            cudaMC.MCrun(cudaSim.gpuLattice, beta, sub);
            stopwatch.add("montecarlo");
   //printf("HERE - 4\n");

             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(cudaSim.gpuLattice);
            stopwatch.add("hamiltonian");
         }

         // Check for error
         cudaError_t e = cudaGetLastError();
         if(e != cudaSuccess) {
            std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
            cudaDeviceReset();
            std::exit(EXIT_FAILURE);
         }
        // printf("mcs = %i\n", mstep);

   }  
   // End loop over simulation steps
   // Synchronize with device
   //cudaDeviceSynchronize();
   //printf("HERE - 3\n");
   //printf("HERE - 5\n");

   cudaMC.mom_update(cudaSim.gpuLattice);
   //printf("HERE - 6\n");

   cudaDeviceSynchronize();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (strict)
void CudaSimulation::CudaMCSimulation::MCmphase(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC measurement phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // MC
   CudaMetropolis cudaMC;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   num_subL = cudaMC.initiate(cudaSim.SimParam, cudaSim.cpuHamiltonian, cudaSim.cpuLattice);

   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   // Measurement
   const auto measurement = MeasurementFactory::create(cudaSim.gpuLattice, cudaSim.cpuLattice);
 
   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (cudaSim.SimParam.Temp * cudaSim.SimParam.k_bolt);



   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = cudaSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement->measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, cudaSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(cudaSim.gpuLattice);
      //stopwatch.add("hamiltonian");
      for(unsigned int sub; sub < num_subL; sub++){
         // Perform Metropolis sweep
         cudaMC.MCrun(cudaSim.gpuLattice, beta, sub);
         stopwatch.add("montecarlo");
         hamCalc.heisge(cudaSim.gpuLattice);
         stopwatch.add("hamiltonian");


      }



      // Check for error
      cudaError_t e = cudaGetLastError();
      if(e != cudaSuccess) {
         std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
         cudaDeviceReset();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
   measurement->measure(mcnstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement->flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   cudaDeviceSynchronize();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo initial phase (brute force)
void CudaSimulation::CudaMCSimulation::MCiphase_bf(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC BF initial phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // Metropolis
   CudaMetropolis_bruteforce cudaMC_bf;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;
   

   // Initiate MC and Hamiltonian
   //num_subL = cudaMC_bf.initiate(cudaSim.SimParam, cudaSim.cpuHamiltonian, cudaSim.cpuLattice);

   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!cudaMC_bf.initiate(cudaSim.SimParam)) {
      std::fprintf(stderr, "CudaMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = cudaSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = cudaSim.cpuLattice.ipmcnstep(it);
   beta = 1/(cudaSim.cpuLattice.ipTemp(it) * cudaSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(cudaSim.gpuLattice);
   stopwatch.add("hamiltonian");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, cudaSim, mcs);
      //printMdStatus(mstep, cudaSim);
         // Perform Metropolis sweep

            cudaMC_bf.MCrun(cudaSim.gpuLattice, beta);
            stopwatch.add("montecarlo");
             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(cudaSim.gpuLattice);
            stopwatch.add("hamiltonian");
         

         // Check for error
         cudaError_t e = cudaGetLastError();
         if(e != cudaSuccess) {
            std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
            cudaDeviceReset();
            std::exit(EXIT_FAILURE);
         }
        // printf("mcs = %i\n", mstep);

   }  
   // End loop over simulation steps
   // Synchronize with device
   //cudaDeviceSynchronize();
   //printf("HERE - 3\n");
   cudaMC_bf.mom_update(cudaSim.gpuLattice);
   cudaDeviceSynchronize();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (brute force)
void CudaSimulation::CudaMCSimulation::MCmphase_bf(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC BF measurement phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // MC
   CudaMetropolis_bruteforce cudaMC_bf;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!cudaMC_bf.initiate(cudaSim.SimParam)) {
      std::fprintf(stderr, "CudaMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

 // Measurement
   const auto measurement = MeasurementFactory::create(cudaSim.gpuLattice, cudaSim.cpuLattice);
 
   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (cudaSim.SimParam.Temp * cudaSim.SimParam.k_bolt);



   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = cudaSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement->measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, cudaSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(cudaSim.gpuLattice);
      //stopwatch.add("hamiltonian");

         // Perform Metropolis sweep
         cudaMC_bf.MCrun(cudaSim.gpuLattice, beta);
         stopwatch.add("montecarlo");
         hamCalc.heisge(cudaSim.gpuLattice);
         stopwatch.add("hamiltonian");


      



      // Check for error
      cudaError_t e = cudaGetLastError();
      if(e != cudaSuccess) {
         std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
         cudaDeviceReset();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
   measurement->measure(mcnstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement->flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   cudaDeviceSynchronize();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}
