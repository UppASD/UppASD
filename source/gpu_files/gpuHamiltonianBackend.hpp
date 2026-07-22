#pragma once

enum class GpuHamiltonianBackend {
   DirectSparse,
   LatticeConvolution,
   FftDipole,
   GridDipole,
   MultiscaleDipole
};

struct GpuHamiltonianBackendSelection {
   GpuHamiltonianBackend exchange = GpuHamiltonianBackend::DirectSparse;
   GpuHamiltonianBackend dmi = GpuHamiltonianBackend::DirectSparse;
   GpuHamiltonianBackend dipole = GpuHamiltonianBackend::DirectSparse;

   bool convolution_ready = false;
   bool multiscale_ready = false;
   bool convolution_kernel_ready = false;
};
