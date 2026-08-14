// TODO Fix that THREAD_COUNT so its not hard coded

#pragma once
#include "gpu_wrappers.h"
#include "c_headers.hpp"
//#include "cudaGPUErrchk.hpp"
#include "gridHelper.hpp"
#include "real_type.h"
#if !defined(THREAD_COUNT)
#define THREAD_COUNT 256
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

// The list contains one-based physical site identifiers.  A compact-list slot
// is deliberately not an atom identity: the same site is visited once for
// every ensemble and then translated back to the ordinary flattened (site +
// ensemble*N) storage index before the operation sees it.
template <std::size_t threads, bool big, typename O>
__global__ void active_atom_kernel(O op, const int* oneBasedAtoms,
                                   unsigned int activeCount, unsigned int N,
                                   unsigned int M) {
   unsigned int work;
   if(GridHelper<threads, big>::index1d(&work, activeCount * M)) {
      const unsigned int slot = work % activeCount;
      const unsigned int ensemble = work / activeCount;
      const int oneBasedAtom = oneBasedAtoms[slot];
      op.each(static_cast<unsigned int>(oneBasedAtom - 1) + ensemble * N);
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
   // Execute an Atom operation for a compact list of one-based physical sites
   // in every ensemble.  The list must contain valid, unique site ids in
   // [1,N]; callers retain ownership because adaptive compaction already owns
   // its device-resident lists.
   template <typename O> void gpuActiveAtomCall(O op, const int* oneBasedAtoms,
                                                unsigned int activeCount);
   template <typename O> void gpuSiteCall(O op);
   template <typename O> void gpuAtomSiteCall(O op);
   template <typename O> void gpuAtomSiteEnsembleCall(O op);
   template <typename O> void gpuElementAxisSiteEnsembleCall(O op); 
   template <typename O> void gpuElementCall(O op);
};

extern GpuParallelizationHelper ParallelizationHelperInstance;



#include "gpuParallelizationHelper.tpp"
