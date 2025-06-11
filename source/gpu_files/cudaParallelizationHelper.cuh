#pragma once
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


