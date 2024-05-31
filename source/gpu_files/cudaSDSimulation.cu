#include <cuda.h>
#include <curand.h>

#include "c_headers.hpp"
#include "c_helper.h"
#include "cudaDepondtIntegrator.hpp"
#include "cudaGPUErrchk.hpp"
#include "tensor.cuh"
#include "cudaHamiltonianCalculations.hpp"

#include "cudaStructures.hpp"
#include "cudaSimulation.hpp"
#include "cudaMeasurement.hpp"
#include "cudaMomentUpdater.hpp"
#include "cudaParallelizationHelper.hpp"
#include "fortranData.hpp"
#include "real_type.h"
#include "stopwatch.hpp"
#include "stopwatchDeviceSync.hpp"
#include "stopwatchPool.hpp"

CudaSimulation::CudaSDSimulation::CudaSDSimulation() {
   //isInitiatedSD = false;
}

CudaSimulation::CudaSDSimulation::~CudaSDSimulation() {
}

// Printing simulation status
// Added copy to fortran line so that simulation status is printed correctly > Thomas Nystrand 14/09/09
void CudaSimulation::CudaSDSimulation::printMdStatus(std::size_t mstep, CudaSimulation& cudaSim)  {
   if(cudaSim.SimParam.nstep > 20) {
      if(mstep % ((cudaSim.SimParam.rstep + cudaSim.SimParam.nstep) / 20) == 0) {
         cudaSim.copyToFortran();  // This is run so seldomly it has not impact on overall performance
         fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
         std::printf(
             "CUDA: %3ld%% done. Mbar: %10.6f. U: %8.5f.\n", mstep * 100 / (cudaSim.SimParam.rstep + cudaSim.SimParam.nstep), *cudaSim.SimParam.mavg, *cudaSim.SimParam.binderc);
      }
   } else {
      cudaSim.copyToFortran();
      fortran_calc_simulation_status_variables(cudaSim.SimParam.mavg);
      std::printf("CUDA: Iteration %ld Mbar %13.6f\n", mstep, *cudaSim.SimParam.mavg);
   }
}

// Spin Dynamics measurement phase
void CudaSimulation::CudaSDSimulation::SDiphase(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaSDSimulation: SD initial phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);
   // Depontd integrator
   CudaDepondtIntegrator integrator;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;

   // Moment updater
   CudaMomentUpdater momUpdater(cudaSim.gpuLattice, cudaSim.SimParam.mompar, cudaSim.SimParam.initexc); 


   // Initiate integrator and Hamiltonian
   if(!integrator.initiate(cudaSim.SimParam)) {
   std::fprintf(stderr, "CudaMdSimulation: integrator failed to initiate!\n");
   return;
   }
   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {
      std::fprintf(stderr, "CudaSDSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   // TEMPORARY PRINTING
   std::printf("\n");
   std::printf("________DEBUG System Information:___________ \n");
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(0));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(1));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(2));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(3));
   std::printf("______________________________________\n");
   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   std::printf("_______________________________________________\n");

   // Initiate constants for integrator
   integrator.initiateConstants(cudaSim.SimParam, cudaSim.cpuLattice);
   // Timing
   stopwatch.add("initiate");
   size_t nstep = cudaSim.SimParam.nstep;
   size_t rstep = cudaSim.SimParam.rstep;
   // Time step loop
   for(std::size_t mstep = rstep + 1; mstep <= rstep + nstep; mstep++) {
      // Print simulation status for each 5% of the simulation length
      //printMdStatus(mstep); -- Do we need it in initial phase?

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(cudaSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
     integrator.evolveFirst(cudaSim.gpuLattice);
     stopwatch.add("evolution");

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(cudaSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(cudaSim.gpuLattice);
      stopwatch.add("evolution");

      // Update magnetic moments after time evolution step
      momUpdater.update();
      stopwatch.add("moments");

      // Check for error
      cudaError_t e = cudaGetLastError();
      if(e != cudaSuccess) {
         std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
         cudaDeviceReset();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps
   // Synchronize with device
   cudaDeviceSynchronize();
   stopwatch.add("final synchronize");
}

// Spin Dynamics measurement phase
void CudaSimulation::CudaSDSimulation::SDmphase(CudaSimulation& cudaSim) {
   // Unbuffered printf
   std::setbuf(stdout, nullptr);
   std::setbuf(stderr, nullptr);
   std::printf("CudaSDSimulation: SD measurement phase starting\n");

   // Initiated?
   if(!cudaSim.isInitiated) {
      std::fprintf(stderr, "CudaSimulation: not initiated!\n");
      return;
   }

   // Timer
   StopwatchDeviceSync stopwatch = StopwatchDeviceSync(GlobalStopwatchPool::get("Cuda measurement phase"));

   // Initiate default parallelization helper
   CudaParallelizationHelper::def.initiate(cudaSim.SimParam.N, cudaSim.SimParam.M, cudaSim.SimParam.NH);
   // Depontd integrator
   CudaDepondtIntegrator integrator;

   // Hamiltonian calculations
   CudaHamiltonianCalculations hamCalc;

   // Moment updater
   CudaMomentUpdater momUpdater(cudaSim.gpuLattice, cudaSim.SimParam.mompar, cudaSim.SimParam.initexc); 
   // Measurement
  CudaMeasurement measurement(cudaSim.gpuLattice.emomM, cudaSim.gpuLattice.emom, cudaSim.gpuLattice.mmom, cudaSim.cpuLattice.emomM, cudaSim.cpuLattice.emom, cudaSim.cpuLattice.mmom); 

   // Initiate integrator and Hamiltonian
   if(!integrator.initiate(cudaSim.SimParam)) { //TODO
      std::fprintf(stderr, "CudaMdSimulation: integrator failed to initiate!\n");
      return;
   }
  
   if(!hamCalc.initiate(cudaSim.Flags, cudaSim.SimParam, cudaSim.gpuHamiltonian)) {//TODO
      std::fprintf(stderr, "CudaSDSimulation: Hamiltonian failed to initiate!\n");
      return;
   }

   // TEMPORARY PRINTING
   std::printf("\n");
   std::printf("________DEBUG System Information:___________ \n");
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(0));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(1));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(2));
   std::printf("%zu\n", cudaSim.cpuHamiltonian.j_tensor.extent(3));
   std::printf("______________________________________\n");
   int mnn = cudaSim.cpuHamiltonian.j_tensor.extent(2);
   int l = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   int NH = cudaSim.cpuHamiltonian.j_tensor.extent(3);
   std::printf("_______________________________________________\n");

   // Initiate constants for integrator
   integrator.initiateConstants(cudaSim.SimParam, cudaSim.cpuLattice);


   // Timing
   stopwatch.add("initiate");

   size_t nstep = cudaSim.SimParam.nstep;
   size_t rstep = cudaSim.SimParam.rstep;

   // Time step loop
   for(std::size_t mstep = rstep + 1; mstep <= rstep + nstep; mstep++) {
       // Measure
      measurement.measure(mstep);
      stopwatch.add("measurement");

      // Print simulation status for each 5% of the simulation length
      printMdStatus(mstep, cudaSim);

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(cudaSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform first step of SDE solver
      integrator.evolveFirst(cudaSim.gpuLattice);
      stopwatch.add("evolution");

      // Apply Hamiltonian to obtain effective field
      hamCalc.heisge(cudaSim.gpuLattice);
      stopwatch.add("hamiltonian");

      // Perform second (corrector) step of SDE solver
      integrator.evolveSecond(cudaSim.gpuLattice);
      stopwatch.add("evolution");
      // Update magnetic moments after time evolution step
      momUpdater.update();
      stopwatch.add("moments");

      // Check for error
      cudaError_t e = cudaGetLastError();
      if(e != cudaSuccess) {
         std::printf("Uncaught CUDA error %d: %s\n", e, cudaGetErrorString(e));
         cudaDeviceReset();
         std::exit(EXIT_FAILURE);
      }

   }  // End loop over simulation steps

   // Final measure
    measurement.measure(rstep + nstep + 1);//TODO
   stopwatch.add("measurement");

   // Print remaining measurements
   measurement.flushMeasurements(rstep + nstep + 1);//TODO
   stopwatch.add("flush measurement");

   // Synchronize with device
   cudaDeviceSynchronize();
   stopwatch.add("final synchronize");
}

