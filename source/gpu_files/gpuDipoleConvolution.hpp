#pragma once

#include "gpu_wrappers.h"
#include "gpuFftWrapper.hpp"
#include "dipoleEwaldKernel.hpp"
#include "real_type.h"
#include "stopwatchDeviceSync.hpp"
#include "tensor.hpp"

#include <cstddef>
#include <array>
#include <memory>
#include <vector>

// CV6 deliberately separates the grid used by the dipole solver from the
// exchange-convolution grid.  Dipole boundary conditions define the physical
// long-range sum and must not be inferred from exchange's current PBC-only
// convolution eligibility.
enum class GpuDipoleBoundaryMode {
   Open,
   Periodic3D,
   Periodic2D
};

enum class GpuDipoleDiscretization {
   Atomistic,
   MacrospinGrid
};

struct GpuDipoleGridShape {
   std::size_t n1 = 0;
   std::size_t n2 = 0;
   std::size_t n3 = 0;

   std::size_t cells() const;
   bool valid() const;
};

// Buffer dimensions for the regular-grid R2C/C2R path.  Keeping this as a
// value type makes the FFT-library integer limits and the memory contract
// testable before a plan or a device allocation exists.
struct GpuDipoleFftLayout {
   // The source macro tensor has active_cells cells.  The real R2C allocation
   // has fft_cells cells and may include zero padding for a finite test
   // convolution.  They coincide for the accepted periodic operator.
   GpuDipoleGridShape active_grid{};
   GpuDipoleGridShape fft_grid{};
   GpuDipoleGridShape spectral_grid{};
   std::size_t active_cells = 0;
   std::size_t active_macros = 0;
   std::size_t fft_cells = 0;
   std::size_t spectral_cells = 0;
   std::size_t field_batches = 0;   // 3 * basis * ensembles
   std::size_t kernel_batches = 0;  // 9 * basis * basis

   bool valid() const;
   bool fitsFftLibrary() const;
   std::size_t persistentBytes() const;
   std::size_t constructionBytes() const;
};

// Initialization-only diagnostics for the complete uploaded spectrum.  They
// are intentionally separate from evaluate(): normal runtime evaluation must
// not read back or synchronize device data.
struct GpuDipoleSpectrumDiagnostics {
   double max_reciprocity_error = 0.0;
   double max_conjugacy_error = 0.0;
   double max_hermitian_error = 0.0;
};

// Process-local Builder-B cache diagnostics.  The cache retains immutable host
// tensors only; every convolution instance still owns and uploads its own GPU
// spectrum and FFT plans.  This makes repeated equivalent initializations
// cheaper without coupling lifetimes across devices or streams.
struct GpuDipoleKernelCacheStats {
   std::size_t hits = 0;
   std::size_t misses = 0;
   std::size_t entries = 0;
};

struct GpuDipoleConvolutionDescriptor {
   // Atomistic grid dimensions.  For MacrospinGrid these describe the source
   // lattice; macro_grid describes the active macro source grid.
   GpuDipoleGridShape atomistic_grid{};
   GpuDipoleGridShape macro_grid{};
   // Empty selects the canonical grid for the boundary mode.  A non-empty
   // value is used by diagnostic OPEN_FFT plumbing only and is validated
   // against the selected active grid.
   GpuDipoleGridShape fft_grid{};
   unsigned int basis = 0;
   unsigned int ensembles = 0;
   GpuDipoleBoundaryMode boundary = GpuDipoleBoundaryMode::Open;
   GpuDipoleDiscretization discretization = GpuDipoleDiscretization::Atomistic;

   // Primitive-cell vectors.  fullCellMatrix() scales these by atomistic_grid
   // to form the periodic supercell H; Cartesian grid extents alone are never
   // a valid substitute for that (possibly skew) matrix.
   const real* c1 = nullptr;
   const real* c2 = nullptr;
   const real* c3 = nullptr;
   // Cartesian primitive-cell basis offsets.  The descriptor owns an fp64
   // copy so Ewald construction does not depend on Fortran buffer lifetime.
   std::vector<std::array<double, 3>> basis_offsets;
   // Device-resident [3, macrocell] centres from the CPU-owned PME map.
   const real* macro_centers = nullptr;
   std::size_t macro_count = 0;
   // Physical conversion is deliberately fp64 and applied only by the future
   // production field boundary.  The dormant kernel remains dimensionless.
   double alat = 0.0;
   double tolerance = 1.0e-10;
   double field_prefactor = 0.0;

   // Only used for Periodic2D.  The selected axis is the non-periodic slab
   // normal and is linearly padded; the remaining axes are periodic.
   unsigned int open_axis = 2;

   GpuDipoleGridShape activeGrid() const;
   GpuDipoleGridShape fftGrid() const;
   std::array<real, 9> fullCellMatrix() const;
   real cellVolume() const;
   std::array<real, 9> reciprocalCellMatrix() const;
   GpuDipoleFftLayout fftLayout() const;
   bool valid() const;
};

// One bridge value is consumed by both the memory preflight and runtime plan
// allocation.  Keeping this conversion here prevents those paths from quietly
// diverging as the descriptor gains physics or storage fields.
struct GpuDipoleDescriptorInput {
   GpuDipoleGridShape atomistic_grid{};
   GpuDipoleGridShape macro_grid{};
   unsigned int basis = 0;
   unsigned int ensembles = 0;
   std::array<char, 3> boundaries{{'0', '0', '0'}};
   const real* c1 = nullptr;
   const real* c2 = nullptr;
   const real* c3 = nullptr;
   const real* basis_offsets = nullptr; // [3, basis], Fortran column-major
   const real* macro_centers = nullptr;
   std::size_t macro_count = 0;
   double alat = 0.0;
   double tolerance = 1.0e-10;
};

bool makeEwald3dFftDipoleDescriptor(const GpuDipoleDescriptorInput& input,
                                    GpuDipoleConvolutionDescriptor& descriptor);
// OPEN_FFT is a separate finite Hamiltonian, not a periodic descriptor with a
// different padding choice.  Its factory intentionally accepts only all-open
// boundary flags and selects the active/padded open layout.
bool makeOpenFftDipoleDescriptor(const GpuDipoleDescriptorInput& input,
                                 GpuDipoleConvolutionDescriptor& descriptor);

// Lifecycle and geometry shell for the regular-grid operator.  Plan readiness
// and kernel readiness are deliberately separate: a valid allocation must
// never be mistaken for a complete periodic dipole operator.
class GpuDipoleConvolution {
public:
   GpuDipoleConvolution() = default;
   ~GpuDipoleConvolution();
   GpuDipoleConvolution(const GpuDipoleConvolution&) = delete;
   GpuDipoleConvolution& operator=(const GpuDipoleConvolution&) = delete;

   bool initiate(const GpuDipoleConvolutionDescriptor& descriptor, GPU_STREAM_T work_stream);
   void release();

   bool isInitiated() const;
   bool kernelReady() const;
   const GpuDipoleConvolutionDescriptor& descriptor() const;
   const GpuDipoleFftLayout& fftLayout() const;

   // Upload the complete Builder A displacement tensor, transform every
   // kernel batch, and store FFT(K)/Ngrid.  The input layout is
   // [cell + Ngrid*kernel_batch], with n1 fastest and
   // kernel_batch = row + 3*(column + 3*(target_basis + basis*source_basis)).
   // This construction API remains outside production dispatch.
   // ``validate_physics`` is false only for deliberately non-reciprocal delta
   // tensors in the FFT-plumbing suite.  Complete Ewald kernels always use
   // the default validation before becoming ready.
   void uploadCompleteKernelForTesting(const std::vector<double>& complete_kernel,
                                       bool validate_physics = true);
   void uploadCompleteKernelForTesting(const DipoleKernelBuildResult& complete_kernel);
   // Production OPEN_FFT uploads the already-built finite displacement tensor
   // once.  The implementation shares the proven transform/normalization
   // path with the diagnostic upload API, but is deliberately named as a
   // production operation so callers cannot mistake it for periodic Builder B.
   void uploadOpenKernel(const std::vector<double>& complete_kernel);

   // Production construction uses the fp64 reciprocal-alias Builder B.
   // Builder A remains available above as the validation/reference path.
   void buildPeriodicKernel();

   // The only runtime operator primitive in this slice.  It performs packed
   // R2C, block contraction, and raw C2R; the stored spectrum owns the sole
   // 1/Ngrid normalization.
   void evaluate(const GpuTensor<real, 3>& macro_moments);

   // The accepted OPEN_FFT block-one slice is basis-resolved.  The field
   // produced by evaluate is scaled from its dimensionless kernel convention exactly once here and
   // added to the already assembled Hamiltonian fields.  Energy is reduced
   // from those exact packed macro fields, in the pre-mRy (Tesla * mu_B)
   // convention used by deviceEnergies::energyM.
   void addFieldsToAtoms(GpuTensor<real, 3>& beff, GpuTensor<real, 3>& eneff,
                         const unsigned int* macro_cell_index,
                         std::size_t atom_count);
   void accumulateEnergy(GpuTensor<real, 2>& energyM, std::size_t atom_count);

   // Explicitly diagnostic readback APIs.  They synchronize and are intended
   // only for the standalone GPU tests while simulation dispatch remains off.
   std::vector<real> diagnosticFieldsForTesting() const;
   // Exact padded FFT-grid buffers, exposed only for layout validation.  They
   // make stale-padding errors observable without changing production APIs.
   std::vector<real> diagnosticPackedMomentsForTesting() const;
   std::vector<real> diagnosticPaddedFieldsForTesting() const;
   // Accumulate the downloaded device-real field and moment values in fp64.
   // This diagnoses the operator error without adding a second, unrelated
   // fp32 host reduction error to the acceptance comparison.
   std::vector<double> diagnosticEnergiesForTesting() const;
   bool diagnosticConstructionStorageAllocatedForTesting() const;
   GpuDipoleSpectrumDiagnostics diagnosticSpectrumForTesting() const;
   static GpuDipoleKernelCacheStats kernelCacheStatsForTesting();
   static void clearKernelCacheForTesting();

   // Persistent field and spectral buffers required by the eventual regular
   // grid solver.  Tensor construction staging is deliberately separate
   // because it is released after the kernel FFT.
   static std::size_t estimatePersistentBytes(const GpuDipoleConvolutionDescriptor& descriptor);
   static std::size_t estimateWorkspaceBytes(const GpuDipoleConvolutionDescriptor& descriptor);
   // Construction staging differs by physical mode.  OPEN_FFT uploads only
   // the finite real-space tensor; the periodic builder also stages its
   // reciprocal alias spectrum.
   static std::size_t estimateConstructionBytes(const GpuDipoleConvolutionDescriptor& descriptor);
   static std::size_t estimateBytes(const GpuDipoleConvolutionDescriptor& descriptor);

private:
   GpuDipoleConvolutionDescriptor desc{};
   GpuDipoleFftLayout layout{};
   GPU_STREAM_T stream{};
   bool initiated = false;
   bool kernel_ready = false;
   bool kernel_real_allocated = false;
   GpuDipoleSpectrumDiagnostics spectrum_diagnostics{};
   bool allocated = false;
   bool forward_plan_created = false;
   bool backward_plan_created = false;
   bool kernel_plan_created = false;
   GpuFftHandle forward_plan{};
   GpuFftHandle backward_plan{};
   GpuFftHandle kernel_plan{};
   GpuTensor<real, 2> moments_real;
   GpuTensor<real, 2> fields_real;
   GpuTensor<GpuFftComplex, 2> moments_fft;
   GpuTensor<GpuFftComplex, 2> fields_fft;
   GpuTensor<GpuFftComplex, 2> kernel_fft;
   // Construction-only [real_cell, kernel_batch] storage.  It is allocated
   // only while the complete host tensor is uploaded and FFT'd, then freed.
   GpuTensor<real, 2> kernel_real;
   GpuTensor<real, 2> cell_vectors;
   GpuTensor<real, 2> reciprocal_vectors;
   GpuTensor<real, 1> cell_volume;
   GpuTensor<unsigned char, 1> fft_workspace;
   std::unique_ptr<StopwatchDeviceSync> stopwatch;

   void packMacroMoments(const GpuTensor<real, 3>& macro_moments);
   void forwardTransformMoments();
   void applySpectralKernel();
   void inverseTransformFields();
   void uploadRealKernelAndAliasSpectrum(const std::vector<double>& real_kernel,
                                         const std::vector<std::complex<double>>& reciprocal_alias_spectrum,
                                         bool validate_physics);
};
