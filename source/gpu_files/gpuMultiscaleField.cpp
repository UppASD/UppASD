#include "gpuMultiscaleField.hpp"

void GpuMultiscaleField::addField(deviceLattice& gpuLattice, GPU_STREAM_T stream) {
   (void)gpuLattice;
   (void)stream;
   // Planned pipeline:
   // 1. accumulate atomistic emomM into cell_moment
   // 2. solve long-range/coarse dipole or micromagnetic field on cell_field
   // 3. prolong cell_field back into atomistic beff/eneff as an additive term
}
