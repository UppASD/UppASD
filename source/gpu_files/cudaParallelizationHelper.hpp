// TODO Fix that THREAD_COUNT so its not hard coded

#pragma once

#include "cudaGPUErrchk.hpp"
#include "gridHelper.hpp"
#include "real_type.h"

#if !defined(THREAD_COUNT)
#define THREAD_COUNT 64
#endif
#if !defined(USE_BIG_GRID)
#define USE_BIG_GRID false
#endif

////////////////////////////////////////////////////////////////////////////////
// Kernel templates
////////////////////////////////////////////////////////////////////////////////
template <usd_int threads, bool big, typename O>
__global__ void atom_kernel(O op) {
   usd_int atom;
   if(GridHelper<threads, big>::index1d(&atom, op.NM)) {
      op.each(atom);
   }
}

template <usd_int threads, bool big, typename O>
__global__ void site_kernel(O op) {
   usd_int site;
   if(GridHelper<threads, big>::index1d(&site, op.N)) {
      op.each(site);
   }
}

template <usd_int threads, bool big, typename O>
__global__ void atom_site_kernel(O op) {
   usd_int site, ensemble;
   if(GridHelper<threads, big>::index2d(&site, &ensemble, op.N, op.M)) {
      usd_int atom = ensemble * op.N + site;
      op.each(atom, site);
   }
}

template <usd_int threads, bool big, typename O>
__global__ void atom_site_ensemble_kernel(O op) {
   usd_int site, ensemble;
   if(GridHelper<threads, big>::index2d(&site, &ensemble, op.N, op.M)) {
      usd_int atom = ensemble * op.N + site;
      op.each(atom, site, ensemble);
   }
}

template <usd_int threads, bool big, typename O>
__global__ void element_axis_site_ensemble_kernel(O op) {
   usd_int axis, site, ensemble;
   if(GridHelper<threads, big>::index3d(&axis, &site, &ensemble, 3, op.N, op.M)) {
      usd_int element = axis + 3 * site + 3 * op.N * ensemble;
      op.each(element, axis, site, ensemble);
   }
}

template <usd_int threads, bool big, typename O>
__global__ void element_kernel(O op) {
   usd_int element;
   if(GridHelper<threads, big>::index1d(&element, op.NM3)) {
      op.each(element);
   }
}

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper class
////////////////////////////////////////////////////////////////////////////////

class CudaParallelizationHelper {
private:
   // System size
   usd_int N;
   usd_int M;

   // Streams
   cudaStream_t workStream;
   cudaStream_t copyStream;

   // Grid helper
   GridHelper<THREAD_COUNT, USE_BIG_GRID> gridHelper;

public:
   // Default helper
   static CudaParallelizationHelper def;

   // Parallelization base classes
   struct Atom {
      usd_int NM;
   };

   struct Site {
      usd_int N;
   };

   struct AtomSite {
      usd_int N;
      usd_int M;
   };

   struct AtomSiteEnsemble {
      usd_int N;
      usd_int M;
   };

   struct ElementAxisSiteEnsemble {
      usd_int N;
      usd_int M;
   };

   struct Element {
      usd_int NM3;
   };

   // Constructors
   CudaParallelizationHelper() {
      // System parameters
      N = 0;
      M = 0;

      // Streams
      workStream = 0;
      copyStream = 0;
   }

   CudaParallelizationHelper(usd_int Natom, usd_int Mensemble) {
      // Zero streams
      workStream = 0;
      copyStream = 0;

      // Initiate
      initiate(Natom, Mensemble);
   }

   // Initiate
   void initiate(usd_int Natom, usd_int Mensemble) {
      // System size
      N = Natom;
      M = Mensemble;

      // Free previous streams
      free();

      // Create streams
      cudaStreamCreate(&workStream);
      cudaStreamCreate(&copyStream);
   }

   // Free
   void free() {
      if(workStream != 0) {
         cudaStreamDestroy(workStream);
      }
      if(copyStream != 0) {
         cudaStreamDestroy(copyStream);
      }
      workStream = copyStream = 0;
   }

   // Stream access
   cudaStream_t getWorkStream() {
      return workStream;
   }

   cudaStream_t getCopyStream() {
      return copyStream;
   }

   // Call helpers
   template <typename O>
   void cudaAtomCall(O op) {
      // Assert that O is derived from Atom
      (void)static_cast<Atom*>((O*)0);

      // Setup size
      op.NM = N * M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N * M);
      atom_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void cudaSiteCall(O op) {
      // Assert that O is derived from Site
      (void)static_cast<Site*>((O*)0);

      // Setup size
      op.N = N;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N);
      site_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void cudaAtomSiteCall(O op) {
      // Assert that O is derived from AtomSite
      (void)static_cast<AtomSite*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim2d(&block, &grid, N, M);
      atom_site_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void cudaAtomSiteEnsembleCall(O op) {
      // Assert that O is derived from AtomSiteEnsemble
      (void)static_cast<AtomSiteEnsemble*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim2d(&block, &grid, N, M);
      atom_site_ensemble_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   template <typename O>
   void cudaElementAxisSiteEnsembleCall(O op) {
      // Assert that O is derived from ElementAxisSiteEnsemble
      (void)static_cast<ElementAxisSiteEnsemble*>((O*)0);

      // Setup size
      op.N = N;
      op.M = M;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim3d(&block, &grid, 3, N, M);
      element_axis_site_ensemble_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
   }

   // Call helpers
   template <typename O>
   void cudaElementCall(O op) {
      // Assert that O is derived from Element
      (void)static_cast<Element*>((O*)0);

      // Setup size
      op.NM3 = N * M * 3;

      // Call kernel
      dim3 block, grid;
      gridHelper.dim1d(&block, &grid, N * M * 3);
      element_kernel<THREAD_COUNT, USE_BIG_GRID><<<grid, block, 0, workStream>>>(op);
      // gpuErrchk(cudaPeekAtLastError());
      // gpuErrchk(cudaDeviceSynchronize());
   }
};

