// TODO Fix that THREAD_COUNT so its not hard coded

#pragma once
#include "gpu_wrappers.h"
#include "c_headers.hpp"
//#include "cudaGPUErrchk.hpp"
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
template <std::size_t threads, bool big, typename O>
__global__ void atom_kernel(O op) {
   unsigned int atom;
   if(GridHelper<threads, big>::index1d(&atom, op.NM)) {
      op.each(atom);
   }
}

template <std::size_t threads, bool big, typename O>
__global__ void site_kernel(O op) {
   unsigned int site;
   if(GridHelper<threads, big>::index1d(&site, op.N)) {
      op.each(site);
   }
}

template <std::size_t threads, bool big, typename O>
__global__ void atom_site_kernel(O op) {
   unsigned int site, ensemble;
   if(GridHelper<threads, big>::index2d(&site, &ensemble, op.N, op.M)) {
      unsigned int atom = ensemble * op.N + site;
      op.each(atom, site);
   }
}

template <std::size_t threads, bool big, typename O>
__global__ void atom_site_ensemble_kernel(O op) {
   unsigned int site, ensemble;
   if(GridHelper<threads, big>::index2d(&site, &ensemble, op.N, op.M)) {
      unsigned int atom = ensemble * op.N + site;
      op.each(atom, site, ensemble);
   }
}

template <std::size_t threads, bool big, typename O>
__global__ void element_axis_site_ensemble_kernel(O op) {
   unsigned int axis, site, ensemble;
   if(GridHelper<threads, big>::index3d(&axis, &site, &ensemble, 3, op.N, op.M)) {
      unsigned int element = axis + 3 * site + 3 * op.N * ensemble;
      op.each(element, axis, site, ensemble);
   }
}

template <std::size_t threads, bool big, typename O>
__global__ void element_kernel(O op) {
   unsigned int element;
   if(GridHelper<threads, big>::index1d(&element, op.NM3)) {
      op.each(element);
   }
}

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper class
////////////////////////////////////////////////////////////////////////////////

class GpuParallelizationHelper {
private:
   // System size
   unsigned int N;
   unsigned int NH;
   unsigned int M;

   // Streams
   GPU_STREAM_T workStream;
   GPU_STREAM_T copyStream;

   // Grid helper
   GridHelper<THREAD_COUNT, USE_BIG_GRID> gridHelper;

public:
   // Default helper
   static GpuParallelizationHelper def;

   // Parallelization base classes
   struct Atom {unsigned int NM;};
   struct Site {unsigned int N; unsigned int NH;};
   struct AtomSite {unsigned int N; unsigned int M;};
   struct AtomSiteEnsemble {unsigned int N; unsigned int NH; unsigned int M;};
   struct ElementAxisSiteEnsemble {unsigned int N; unsigned int M; unsigned int NH;};
   struct Element {unsigned int NM3;};

   // Constructors
   GpuParallelizationHelper();
   GpuParallelizationHelper(unsigned int Natom, unsigned int Mensemble, unsigned int nHam);

   // Initiate
   void initiate(unsigned int Natom, unsigned int Mensemble, unsigned int nHam);
   // Free
   void free();

   // Stream access
   GPU_STREAM_T getWorkStream();
   GPU_STREAM_T getCopyStream();

   // Call helpers
   template <typename O> void gpuAtomCall(O op);
   template <typename O> void gpuSiteCall(O op);
   template <typename O> void gpuAtomSiteCall(O op);
   template <typename O> void gpuAtomSiteEnsembleCall(O op);
   template <typename O> void gpuElementAxisSiteEnsembleCall(O op); 
   template <typename O> void gpuElementCall(O op);
};

extern GpuParallelizationHelper ParallelizationHelperInstance;



#if defined(HIP_V)
#include "hipParallelizationHelper.tpp"
#elif defined(CUDA_V)
#include "cudaParallelizationHelper.tpp"
#endif