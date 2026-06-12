#pragma once

#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "tensor.hpp"

struct GpuMultiscaleDescriptor {
   unsigned int atoms = 0;
   unsigned int ensembles = 0;
   unsigned int cells = 0;
   unsigned int max_atoms_per_cell = 0;
};

struct GpuMultiscaleData {
   GpuTensor<unsigned int, 1> atom_cell_index;
   GpuTensor<unsigned int, 1> cell_atom_count;
   GpuTensor<unsigned int, 2> cell_atom_list;
   GpuTensor<real, 3> cell_moment;
   GpuTensor<real, 3> cell_field;
};

class GpuMultiscaleField {
public:
   bool supports(const GpuMultiscaleDescriptor& descriptor) const {
      return descriptor.atoms > 0 && descriptor.ensembles > 0 &&
             descriptor.cells > 0 && descriptor.max_atoms_per_cell > 0;
   }

   bool initiate(const GpuMultiscaleDescriptor& descriptor, GpuMultiscaleData data) {
      if(!supports(descriptor)) return false;
      desc = descriptor;
      buffers = data;
      initiated = true;
      return true;
   }

   void release() {
      initiated = false;
   }

   bool isInitiated() const {
      return initiated;
   }

   void addField(deviceLattice& gpuLattice, GPU_STREAM_T stream);

private:
   GpuMultiscaleDescriptor desc;
   GpuMultiscaleData buffers;
   bool initiated = false;
};
