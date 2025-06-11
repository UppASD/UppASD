#pragma once
#include <cuda.h>
#include <curand.h>

#include "c_headers.hpp"
#include "c_helper.h"
//#include "cudaGPUErrchk.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "cudaMetropolis.cuh"
#include "cudaMetropolis_bruteforce.cuh"
#include "gpuMeasurement.hpp"
#include "cudaParallelizationHelper.hpp"
#include "gpuSimulation.hpp"
#include "gpuStructures.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"
#include "tensor.hpp"

GpuSimulation::CudaMCSimulation::CudaMCSimulation() {
   // isInitiatedSD = false;
}

GpuSimulation::CudaMCSimulation::~CudaMCSimulation() {
}

void GpuSimulation::CudaMCSimulation::printMdStatus(std::size_t mstep, GpuSimulation& gpuSim) {
   if(gpuSim.SimParam.mcnstep > 20) {
      if(mstep % ((gpuSim.SimParam.rstep + gpuSim.SimParam.mcnstep) / 20) == 0) {
         gpuSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
         std::printf("CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (gpuSim.SimParam.rstep + gpuSim.SimParam.mcnstep),
                     *gpuSim.SimParam.mavg,
                     *gpuSim.SimParam.binderc);
      }
   } else {
      gpuSim.copyToFortran();
      fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
      std::printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *gpuSim.SimParam.mavg);
   }
}

void GpuSimulation::CudaMCSimulation::printMdStatus_iphase(std::size_t mstep, GpuSimulation& gpuSim, int step) {
   if(step > 20) {
      if(mstep % ((step) / 20) == 0) {
         //cudaMC.mom_update(gpuSim.gpuLattice);

         gpuSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
         std::printf("CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n",
                     mstep * 100 / (step),
                     *gpuSim.SimParam.mavg,
                     *gpuSim.SimParam.binderc);
      }
   } else {
      //cudaMC.mom_update(GpuSim.gpuLattice);
      gpuSim.copyToFortran();
      fortran_calc_simulation_status_variables(gpuSim.SimParam.mavg);
      std::printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *gpuSim.SimParam.mavg);
   }
}

// Monte Carlo initial phase (strict)
void GpuSimulation::CudaMCSimulation::MCiphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC initial phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // Metropolis
   CudaMetropolis cudaMC;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;
   //printf("HERE - 0\n");

   // Initiate MC and Hamiltonian
   num_subL = cudaMC.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);
   //printf("HERE - 1\n");
   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = gpuSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);
   //printf("HERE - 2\n");

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = gpuSim.cpuLattice.ipmcnstep(it);
   beta = 1/(gpuSim.cpuLattice.ipTemp(it) * gpuSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(gpuSim.gpuLattice);
   stopwatch.add("hamiltonian");
   //printf("HERE - 3\n");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, gpuSim, mcs);
      //printMdStatus(mstep, gpuSim);
         // Perform Metropolis sweep
         for(std::size_t sub = 0; sub < num_subL; sub++){
            cudaMC.MCrun(gpuSim.gpuLattice, beta, sub);
            stopwatch.add("montecarlo");
   //printf("HERE - 4\n");

             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(gpuSim.gpuLattice);
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

   cudaMC.mom_update(gpuSim.gpuLattice);
   //printf("HERE - 6\n");

   cudaDeviceSynchronize();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (strict)
void GpuSimulation::CudaMCSimulation::MCmphase(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // MC
   CudaMetropolis cudaMC;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   num_subL = cudaMC.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);

   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

 // Measurement
   GpuMeasurement measurement(gpuSim.gpuLattice.emomM,
                               gpuSim.gpuLattice.emom,
                               gpuSim.gpuLattice.mmom,
                               gpuSim.cpuLattice.emomM,
                               gpuSim.cpuLattice.emom,
                               gpuSim.cpuLattice.mmom);
 
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (gpuSim.SimParam.Temp * gpuSim.SimParam.k_bolt);



   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = gpuSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement.measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(gpuSim.gpuLattice);
      //stopwatch.add("hamiltonian");
      for(unsigned int sub; sub < num_subL; sub++){
         // Perform Metropolis sweep
         cudaMC.MCrun(gpuSim.gpuLattice, beta, sub);
         stopwatch.add("montecarlo");
         hamCalc.heisge(gpuSim.gpuLattice);
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
   measurement.measure(mcnstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement.flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   cudaDeviceSynchronize();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo initial phase (brute force)
void GpuSimulation::CudaMCSimulation::MCiphase_bf(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC BF initial phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }

   // Initiate default parallelization helper
    ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // Metropolis
   CudaMetropolis_bruteforce cudaMC_bf;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   

   // Initiate MC and Hamiltonian
   //num_subL = cudaMC_bf.initiate(gpuSim.SimParam, gpuSim.cpuHamiltonian, gpuSim.cpuLattice);

   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!cudaMC_bf.initiate(gpuSim.SimParam)) {
      std::fprintf(stderr, "CudaMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int ipmcnphase = gpuSim.SimParam.ipmcnphase;
   real beta; 
   unsigned int mcs;
   // Timing
   stopwatch.add("initiate");
   printf("\n\\ninphase = %i\n\n", ipmcnphase);

   for(unsigned int it = 0; it < ipmcnphase; it++){
   mcs = gpuSim.cpuLattice.ipmcnstep(it);
   beta = 1/(gpuSim.cpuLattice.ipTemp(it) * gpuSim.SimParam.k_bolt);
   // Apply Hamiltonian to obtain effective field
   hamCalc.heisge(gpuSim.gpuLattice);
   stopwatch.add("hamiltonian");

   //printf("mcs = %i\n", mcs);   
   // Time step loop
   for(std::size_t mstep = 0; mstep <= mcs; mstep++) {
      printMdStatus_iphase(mstep, gpuSim, mcs);
      //printMdStatus(mstep, gpuSim);
         // Perform Metropolis sweep

            cudaMC_bf.MCrun(gpuSim.gpuLattice, beta);
            stopwatch.add("montecarlo");
             // Apply Hamiltonian to obtain effective field
            hamCalc.heisge(gpuSim.gpuLattice);
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
   cudaMC_bf.mom_update(gpuSim.gpuLattice);
   cudaDeviceSynchronize();   
   }



   //cudaMC.release();
   stopwatch.add("final synchronize");
}

// Monte Carlo measurement phase (brute force)
void GpuSimulation::CudaMCSimulation::MCmphase_bf(GpuSimulation& gpuSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaMCSimulation: MC BF measurement phase starting\n");

   // Initiated?
   if(!gpuSim.isInitiated) {
      std::fprintf(stderr, "GpuSimulation: not initiated!\n");
      return;
   }
   // Initiate default parallelization helper
   ParallelizationHelperInstance.initiate(gpuSim.SimParam.N, gpuSim.SimParam.M, gpuSim.SimParam.NH);
   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

    // MC
   CudaMetropolis_bruteforce cudaMC_bf;

   // Hamiltonian calculations
   GpuHamiltonianCalculations hamCalc;
   unsigned int num_subL = 0;

 
   // Initiate MC and Hamiltonian
   if(!hamCalc.initiate(gpuSim.Flags, gpuSim.SimParam, gpuSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaMCSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   if(!cudaMC_bf.initiate(gpuSim.SimParam)) {
      std::fprintf(stderr, "CudaMCSimulation_bf: Hamiltonian failed to initiate!\n");
      return;
   }

 // Measurement
   GpuMeasurement measurement(gpuSim.gpuLattice.emomM,
                               gpuSim.gpuLattice.emom,
                               gpuSim.gpuLattice.mmom,
                               gpuSim.cpuLattice.emomM,
                               gpuSim.cpuLattice.emom,
                               gpuSim.cpuLattice.mmom);
 
   int mnn = gpuSim.cpuHamiltonian.j_tensor.extent(2);
   int l = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = gpuSim.cpuHamiltonian.j_tensor.extent(3);
   real beta = 1 / (gpuSim.SimParam.Temp * gpuSim.SimParam.k_bolt);



   // Timing
   stopwatch.add("initiate");

   size_t mcnstep = gpuSim.SimParam.mcnstep;

   // Time step loop
   for(std::size_t mstep = 1; mstep <= mcnstep; mstep++) {
      // Measure
      //printf("STEP = %i\n", mstep);
      measurement.measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, gpuSim);

      // Apply Hamiltonian to obtain effective field
      //hamCalc.heisge(gpuSim.gpuLattice);
      //stopwatch.add("hamiltonian");

         // Perform Metropolis sweep
         cudaMC_bf.MCrun(gpuSim.gpuLattice, beta);
         stopwatch.add("montecarlo");
         hamCalc.heisge(gpuSim.gpuLattice);
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
   measurement.measure(mcnstep + 1);  // TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement.flushMeasurements(mcnstep + 1);  // TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   cudaDeviceSynchronize();
   //cudaMC.release();
   stopwatch.add("final synchronize");
}
