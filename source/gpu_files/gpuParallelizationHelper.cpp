#pragma once
#include "gpuParallelizationHelper.hpp"
#include "gpu_wrappers.h"


GpuParallelizationHelper ParallelizationHelperInstance;
//CudaParallelizationHelper CudaParallelizationHelper::def;

  // Constructors
  GpuParallelizationHelper::GpuParallelizationHelper() {
      // System parameters
      N = 0;
      M = 0;
      NH = 0;

      // Streams
      workStream = 0;
      copyStream = 0;
   }

   GpuParallelizationHelper::GpuParallelizationHelper(unsigned int Natom, unsigned int Mensemble, unsigned int nHam) {
      // Zero streams
      workStream = 0;
      copyStream = 0;

      // Initiate
      initiate(Natom, Mensemble, nHam);
   }

   // Initiate
   void GpuParallelizationHelper::initiate(unsigned int Natom, unsigned int Mensemble, unsigned int nHam) {
      // System size
      N = Natom;
      M = Mensemble;
      NH = nHam;

      // Free previous streams
      free();

      // Create streams
      GPU_STREAM_CREATE(&workStream);
      GPU_STREAM_CREATE(&copyStream);
   }

   // Free
   void GpuParallelizationHelper::free() {
      if(workStream != 0) {
         GPU_STREAM_DESTROY(workStream);
      }
      if(copyStream != 0) {
         GPU_STREAM_DESTROY(copyStream);
      }
      workStream = copyStream = 0;
   }

   // Stream access
   GPU_STREAM_T GpuParallelizationHelper::getWorkStream() {
      return workStream;
   }

   GPU_STREAM_T GpuParallelizationHelper::getCopyStream() {
      return copyStream;
   }

