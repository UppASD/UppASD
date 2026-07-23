
#include "c_headers.hpp"
#include "gpuHamiltonianDeviceOps.hpp"
#include "gpuHamiltonianCalculations.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "gpuParallelizationHelper.hpp"
#include "measurementData.h"

using ParallelizationHelper = GpuParallelizationHelper;
namespace hamdev = gpu::hamiltonian;

__device__ void sum_warp_energy(real& val){
   #pragma unroll
   for (int offset = WARPSIZE / 2; offset > 0; offset >>= 1)
   {
      val += SHFL_DOWN(val, offset);
   }

}
__global__ void null_energy(GpuVector<EnergyData> energy){
   energy(0) = {};
}

// Possible improvements

////////////////////////////////////////////////////////////////////////////////
// Parallelization helper classes
////////////////////////////////////////////////////////////////////////////////

// The neighbour list setup helper
//
// Note (Thomas):
// For Heisenberg
// Class sets everything between neighbours and maxneighbours
// to zero since hamiltonian implementation always runs to max neighbours
class GpuHamiltonianCalculations::SetupNeighbourList : public ParallelizationHelper::Site {
private:
   real* coup;
   unsigned int* pos;
   const unsigned int* size;
   const unsigned int* aham;
   unsigned int mnn;

public:
   SetupNeighbourList(Exchange& ex, HamRed& redHam) {
      coup = ex.coupling.data();
      size = ex.neighbourCount.data();
      pos = ex.neighbourPos.data();
      aham = redHam.redNeibourCount.data();
      mnn = ex.mnn;
   }

   __device__ void each(unsigned int site) {
      unsigned int rsite = aham[site] - 1;  // to fix the indexing dismatch between Fortran and C++
      real* myCoup = &coup[rsite];
      unsigned int* myPos = &pos[site];
      unsigned int mySize = size[rsite];
      if(site < NH) {  // to avoid calculating the same ncoup in case of reduced Hamiltonian too many times
         for(unsigned int i = 0; i < mnn; i++) {
            if(i < mySize) {
               myPos[i * N]--;
            } else {
               myCoup[i * NH] = (real)0.0;
               myPos[i * N] = 0;
            }
         }
      } else {
         for(unsigned int i = 0; i < mnn; i++) {
            if(i < mySize) {
               myPos[i * N]--;
            } else {
               myPos[i * N] = 0;
            }
         }
      }
   }
};

// CV6.1: macro moments must be derived from the current device moments every
// Hamiltonian evaluation.  The host macro arrays are only initialization data;
// using them after GPU integration would leave the dipole field one or more
// time steps stale. cell_index keeps Fortran's one-based convention.
class GpuHamiltonianCalculations::UpdateMacroMoments : public ParallelizationHelper::AtomSite {
private:
   const unsigned int* cellIndex;
   const real* moments;
   real* macroMoments;
   unsigned int macroCount;

public:
   UpdateMacroMoments(const unsigned int* indices, const GpuTensor<real, 3>& source,
                      GpuTensor<real, 3>& destination, unsigned int count)
      : cellIndex(indices), moments(source.data()), macroMoments(destination.data()), macroCount(count) {}

   __device__ void each(unsigned int atom, unsigned int site) {
      const unsigned int oneBasedCell = cellIndex[site];
      if(oneBasedCell == 0 || oneBasedCell > macroCount) return;
      const unsigned int cell = oneBasedCell - 1;
      const unsigned int ensemble = atom / N;
      atomicAdd(&macroMoments[0 + 3 * (cell + macroCount * ensemble)], moments[0 + 3 * (site + N * ensemble)]);
      atomicAdd(&macroMoments[1 + 3 * (cell + macroCount * ensemble)], moments[1 + 3 * (site + N * ensemble)]);
      atomicAdd(&macroMoments[2 + 3 * (cell + macroCount * ensemble)], moments[2 + 3 * (site + N * ensemble)]);
   }
};

// The neighbour list setup helper
//
// For Tensorial Exchange
class GpuHamiltonianCalculations::SetupNeighbourListExchangeTensor : public ParallelizationHelper::Site {
private:
   real* tensor;
   unsigned int* pos;
   const unsigned int* size;
   unsigned int mnn;
   unsigned int* aham;

public:
   SetupNeighbourListExchangeTensor(TensorialExchange& tenEx, HamRed& redHam) {
      tensor = tenEx.tensor.data();
      size = tenEx.neighbourCount.data();
      pos = tenEx.neighbourPos.data();
      mnn = tenEx.mnn;
      aham = redHam.redNeibourCount.data();
   }

   __device__ real& tensor_ind(unsigned int i, unsigned int j, unsigned int k, unsigned int l) {
      return tensor[l + NH * (i + 3 * (j + 3 * k))];
   }

   __device__ void each(unsigned int site) {
      unsigned int rsite = aham[site] - 1;
      if(site < NH) {
         for(unsigned int k = 0; k < mnn; k++) {
            if(pos[site * mnn + k] != 0) {
               pos[site * mnn + k]--;
            } else {
               pos[site * mnn + k] = 0;

               tensor_ind(0, 0, k, rsite) = {};
               tensor_ind(0, 1, k, rsite) = {};
               tensor_ind(0, 2, k, rsite) = {};
               tensor_ind(1, 0, k, rsite) = {};
               tensor_ind(1, 1, k, rsite) = {};
               tensor_ind(1, 2, k, rsite) = {};
               tensor_ind(2, 0, k, rsite) = {};
               tensor_ind(2, 1, k, rsite) = {};
               tensor_ind(2, 2, k, rsite) = {};
            }
         }
      } else {
         for(unsigned int k = 0; k < mnn; k++) {
            if(pos[site * mnn + k] != 0) {
               pos[site * mnn + k]--;
            } else {
               pos[site * mnn + k] = 0;
            }
         }
      }
   }
};

// unnecessary for anisotropy probably
class GpuHamiltonianCalculations::SetupAnisotropy : public ParallelizationHelper::Site {
private:
   real* kaniso;
   real* eaniso;
   unsigned int* taniso;

public:
   SetupAnisotropy(Anisotropy& aniso) {
      kaniso = aniso.kaniso.data();
      eaniso = aniso.eaniso.data();
      taniso = aniso.taniso.data();
   }

   __device__ void each(unsigned int site) {
   }
};

// Note (Thomas):
// For DM interaction
// Class sets everything between neighbours and maxneighbours
// to zero since hamiltonian implementation always runs to max neighbours
class GpuHamiltonianCalculations::SetupNeighbourListDM : public ParallelizationHelper::Site {
private:
   real* coup;
   unsigned int* pos;
   const unsigned int* size;
   unsigned int* aham;
   unsigned int mnn;

public:
   SetupNeighbourListDM(DMinteraction& dm, HamRed& redHam) {
      coup = dm.interaction.data();
      size = dm.neighbourCount.data();
      pos = dm.neighbourPos.data();
      mnn = dm.mnn;
      aham = redHam.redNeibourCount.data();
   }

   __device__ void each(unsigned int site) {
      // Phil's
      unsigned int rsite = aham[site] - 1;
      if(site < NH) {
         for(unsigned int i = 0; i < mnn; i++) {
            if(pos[site * mnn + i] != 0) {
               pos[site * mnn + i]--;
            } else {
               pos[site * mnn + i] = 0;

               // Dimension of the DM vector: (dim1,dim2,dim3)  <--> (3,mnn,N)
               coup[0 + 3 * i + rsite * mnn * 3] = (real)0.0;
               coup[1 + 3 * i + rsite * mnn * 3] = (real)0.0;
               coup[2 + 3 * i + rsite * mnn * 3] = (real)0.0;
            }
         }
      } else {
         for(unsigned int i = 0; i < mnn; i++) {
            if(pos[site * mnn + i] != 0) {
               pos[site * mnn + i]--;
            } else {
               pos[site * mnn + i] = 0;
            }
         }
      }
   }
};

// Note: (Thomas)
// Calculating the magnetic field from various effects
// such as the heisenberg field and DM interactions
// Added DM effect 2014/09/23
template<bool HasDM, bool HasAniso, bool HasTensor, bool Measure>
class GpuHamiltonianCalculations::Heisge : public ParallelizationHelper::AtomSiteEnsemble {
private:
   real* __restrict__ beff;
   real* __restrict__ eneff;
   GpuTensor<real, 2> energyM;
   const real* __restrict__ coup;
   const unsigned int* __restrict__ pos;
   unsigned int mnn;
   const real* __restrict__ dmcoup;
   const unsigned int* __restrict__ dmpos;
   unsigned int dmmnn;
   const real* __restrict__ tensor;
   const unsigned int* __restrict__ tensorpos;
   unsigned int tensormnn;
   const real* __restrict__ emomM;
   const real* __restrict__ ext_f;
   const real* __restrict__ kaniso;
   const real* __restrict__ eaniso;
   const unsigned int* __restrict__ taniso;
   const real* __restrict__ sb;
   const unsigned int* __restrict__ aham;

public:
   Heisge(GpuTensor<real, 3>& p_beff, GpuTensor<real, 3>& p_eneff, GpuTensor<real, 2> p_energyM,
          const GpuTensor<real, 3>& p_emomM, const GpuTensor<real, 3>& p_ext_f,
          const Exchange& p_ex, const DMinteraction& p_dm, const TensorialExchange& p_tenEx,
          const Anisotropy& p_aniso, const HamRed& p_redHam)
       : energyM(p_energyM), beff(p_beff.data()), eneff(p_eneff.data()), emomM(p_emomM.data()),
         ext_f(p_ext_f.data()), coup(p_ex.coupling.data()), pos(p_ex.neighbourPos.data()), mnn(p_ex.mnn),
         dmcoup(p_dm.interaction.data()), dmpos(p_dm.neighbourPos.data()), dmmnn(p_dm.mnn),
         tensor(p_tenEx.tensor.data()), tensorpos(p_tenEx.neighbourPos.data()), tensormnn(p_tenEx.mnn),
         kaniso(p_aniso.kaniso.data()), eaniso(p_aniso.eaniso.data()), taniso(p_aniso.taniso.data()),
         sb(p_aniso.sb.data()), aham(p_redHam.redNeibourCount.data()) {}

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;
      real dm_x = (real)0.0;
      real dm_y = (real)0.0;
      real dm_z = (real)0.0;
      const unsigned int rsite = aham[site] - 1;
      const real* __restrict__ my_emomM = &emomM[ensemble * N * 3];

      if constexpr (HasTensor) {
         for(unsigned int i = 0; i < tensormnn; i++) {
            const unsigned int neighbor = tensorpos[site * tensormnn + i];
            const unsigned int offset = neighbor * 3;
            const unsigned int base = rsite + NH * (9 * i);
            const real Sx = my_emomM[offset + 0];
            const real Sy = my_emomM[offset + 1];
            const real Sz = my_emomM[offset + 2];
            x += tensor[base + NH * 0] * Sx + tensor[base + NH * 3] * Sy + tensor[base + NH * 6] * Sz;
            y += tensor[base + NH * 1] * Sx + tensor[base + NH * 4] * Sy + tensor[base + NH * 7] * Sz;
            z += tensor[base + NH * 2] * Sx + tensor[base + NH * 5] * Sy + tensor[base + NH * 8] * Sz;
         }
      } else {
         const real* __restrict__ site_coup = &coup[rsite];
         const unsigned int* __restrict__ site_pos = &pos[site];
         for(unsigned int i = 0; i < mnn; i++) {
            const unsigned int offset = site_pos[i * N] * 3;
            const real c = site_coup[i * NH];
            x += c * my_emomM[offset + 0];
            y += c * my_emomM[offset + 1];
            z += c * my_emomM[offset + 2];
         }
      }

      if constexpr (HasDM) {
         for(unsigned int i = 0; i < dmmnn; i++) {
            const unsigned int offset = dmpos[site * dmmnn + i] * 3;
            const auto dm_term = hamdev::dm_field(
                hamdev::make_vec3(dmcoup[3 * i + rsite * dmmnn * 3],
                                  dmcoup[1 + 3 * i + rsite * dmmnn * 3],
                                  dmcoup[2 + 3 * i + rsite * dmmnn * 3]),
                hamdev::make_vec3(my_emomM[offset + 0], my_emomM[offset + 1], my_emomM[offset + 2]));
            dm_x += dm_term.x;
            dm_y += dm_term.y;
            dm_z += dm_term.z;
         }
         x += dm_x;
         y += dm_y;
         z += dm_z;
      }

      real ax = (real)0.0;
      real ay = (real)0.0;
      real az = (real)0.0;
      real ax_en = (real)0.0;
      real ay_en = (real)0.0;
      real az_en = (real)0.0;
      real Sx = (real)0.0;
      real Sy = (real)0.0;
      real Sz = (real)0.0;
      if constexpr (HasAniso || Measure) {
         Sx = emomM[atom * 3 + 0];
         Sy = emomM[atom * 3 + 1];
         Sz = emomM[atom * 3 + 2];
      }
      if constexpr (HasAniso) {
         const unsigned int type = taniso[site];
         const real ex = eaniso[site * 3 + 0];
         const real ey = eaniso[site * 3 + 1];
         const real ez = eaniso[site * 3 + 2];
         const real k1 = kaniso[site * 2 + 0];
         const real k2 = kaniso[site * 2 + 1];
         if(type == 1 || type == 7) {
            const real tt1 = Sx * ex + Sy * ey + Sz * ez;
            const real tt3 = (real)2.0 * tt1 * (k1 + (real)2.0 * k2 * ((real)1.0 - tt1 * tt1));
            const real tt3_en = tt1 * (k1 + k2 * ((real)2.0 - tt1 * tt1));
            ax -= tt3 * ex; ay -= tt3 * ey; az -= tt3 * ez;
            ax_en -= tt3_en * ex; ay_en -= tt3_en * ey; az_en -= tt3_en * ez;
         }
         if(type == 2 || type == 7) {
            real k1_cubic = k1;
            real k2_cubic = k2;
            if(type == 7) { k1_cubic *= sb[site]; k2_cubic *= sb[site]; }
            ax += (real)2.0 * k1_cubic * Sx * (Sy * Sy + Sz * Sz) + (real)2.0 * k2_cubic * Sx * Sy * Sy * Sz * Sz;
            ay += (real)2.0 * k1_cubic * Sy * (Sz * Sz + Sx * Sx) + (real)2.0 * k2_cubic * Sy * Sz * Sz * Sx * Sx;
            az += (real)2.0 * k1_cubic * Sz * (Sx * Sx + Sy * Sy) + (real)2.0 * k2_cubic * Sz * Sx * Sx * Sy * Sy;
            ax_en += k1_cubic * Sx * (Sy * Sy + Sz * Sz) / (real)2.0 + k2_cubic * Sx * Sy * Sy * Sz * Sz / (real)3.0;
            ay_en += k1_cubic * Sy * (Sz * Sz + Sx * Sx) / (real)2.0 + k2_cubic * Sy * Sz * Sz * Sx * Sx / (real)3.0;
            az_en += k1_cubic * Sz * (Sx * Sx + Sy * Sy) / (real)2.0 + k2_cubic * Sz * Sx * Sx * Sy * Sy / (real)3.0;
         }
      }

      const real ext_x = ext_f[atom * 3 + 0];
      const real ext_y = ext_f[atom * 3 + 1];
      const real ext_z = ext_f[atom * 3 + 2];
      beff[atom * 3 + 0] = x + ax + ext_x;
      beff[atom * 3 + 1] = y + ay + ext_y;
      beff[atom * 3 + 2] = z + az + ext_z;
      eneff[atom * 3 + 0] = x + ax_en + ext_x;
      eneff[atom * 3 + 1] = y + ay_en + ext_y;
      eneff[atom * 3 + 2] = z + az_en + ext_z;

      if constexpr (Measure) {
         real interaction = (x * Sx + y * Sy + z * Sz) * (real)-0.5;
         if constexpr (HasDM) interaction = ((x - dm_x) * Sx + (y - dm_y) * Sy + (z - dm_z) * Sz) * (real)-0.5;
         real dm_energy = (dm_x * Sx + dm_y * Sy + dm_z * Sz) * (real)-0.5;
         // On-site anisotropy is not a pairwise term, so it must NOT carry the 1/2
         // double-counting factor used for the bilinear exchange/DM sums; the ax_en
         // prefactors are built so E_ani = -(a_en . S) (verified for uniaxial and
         // cubic against the Fortran calc_energy convention). Using -0.5 here yielded
         // exactly half the reference anisotropy energy.
         real anisotropy = (ax_en * Sx + ay_en * Sy + az_en * Sz) * (real)-1.0;
         real external = ext_x * Sx + ext_y * Sy + ext_z * Sz;
         sum_warp_energy(interaction);
         if constexpr (HasDM) sum_warp_energy(dm_energy);
         if constexpr (HasAniso) sum_warp_energy(anisotropy);
         sum_warp_energy(external);
         if((threadIdx.x & (WARPSIZE - 1)) == 0) {
            interaction /= static_cast<real>(N);
            dm_energy /= static_cast<real>(N);
            anisotropy /= static_cast<real>(N);
            external /= static_cast<real>(N);
            const real total = interaction + (HasDM ? dm_energy : (real)0.0) +
                               (HasAniso ? anisotropy : (real)0.0) + external;
            atomicAdd(&energyM(ensemble, HasTensor ? 3 : 0), interaction);
            if constexpr (HasDM) atomicAdd(&energyM(ensemble, 2), dm_energy);
            if constexpr (HasAniso) atomicAdd(&energyM(ensemble, 1), anisotropy);
            atomicAdd(&energyM(ensemble, 4), external);
            atomicAdd(&energyM(ensemble, 5), total);
         }
      }
   }
};

class GpuHamiltonianCalculations::HeisgeJijElement
    : public ParallelizationHelper::ElementAxisSiteEnsemble {
private:
   real* beff;
   real* eneff;
   const real* coup;
   const unsigned int* pos;
   const unsigned int* size;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const unsigned int* aham;

public:
   HeisgeJijElement(GpuTensor<real, 3>& p_beff,GpuTensor<real, 3>& p_eneff, const GpuTensor<real, 3>& p_emomM, const GpuTensor<real, 3>& p_ext_f, const Exchange& ex, const HamRed& redHam) {
      beff = p_beff.data();
      eneff = p_eneff.data();
      emomM = p_emomM.data();
      ext_f = p_ext_f.data();
      coup = ex.coupling.data();
      pos = ex.neighbourPos.data();
      size = ex.neighbourCount.data();

      mnn = ex.mnn;
      aham = redHam.redNeibourCount.data();
   }

   __device__ void each(unsigned int element, unsigned int axis, unsigned int site, unsigned int ensemble) {
      // Field
      real f = (real)0.0;

      // Pointers with fixed indices
      const unsigned int rsite = aham[site] - 1;
      const real* site_coup = &coup[rsite];
      const unsigned int* site_pos = &pos[site];
      const real* ensemble_emomM = &emomM[ensemble * N * 3];

      // Exchange term loop
      //		const unsigned int s = size[i];
      //		for (int j = 0; j < s; j++) {
      for(unsigned int i = 0; i < mnn; i++) {
         unsigned int offset = site_pos[i * N] * 3;
         f += site_coup[i * NH] * ensemble_emomM[offset + axis];
      }

      // Save field
      beff[element] = f + ext_f[element];
      eneff[element] = f + ext_f[element];
   }
};


///////////////////////////////////ext/////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

GpuHamiltonianCalculations::GpuHamiltonianCalculations() : parallel(ParallelizationHelperInstance) {
   initiated = false;
}

GpuHamiltonianCalculations::~GpuHamiltonianCalculations() {
   convolution.release();
   macroMoments.Free();
}

bool GpuHamiltonianCalculations::canUseLatticeConvolution(const Flag Flags, const SimulationParameters SimParam,
                                                          const deviceHamiltonian& gpuHamiltonian) const {
   if(!SimParam.do_gpu_convolution) return false;
   const auto reject = [](const char* reason) {
      std::printf("Gpu: device lattice convolution requested but disabled: %s; using sparse Hamiltonian.\n", reason);
      return false;
   };
   if(SimParam.BC1 != 'P' || SimParam.BC2 != 'P' || SimParam.BC3 != 'P')
      return reject("the FFT backend currently requires periodic boundary conditions in all directions");
   if(SimParam.N1 == 0 || SimParam.N2 == 0 || SimParam.N3 == 0 || SimParam.NA == 0)
      return reject("the cell grid or basis size is zero");
   if(SimParam.N != SimParam.N1 * SimParam.N2 * SimParam.N3 * SimParam.NA)
      return reject("Natom does not match N1*N2*N3*NA");
   if(SimParam.NH != SimParam.NA)
      return reject("the reduced Hamiltonian is required (set do_reduced Y so nHam equals NA)");
   if(gpuHamiltonian.nlist.empty() || gpuHamiltonian.nlistsize.empty())
      return reject("the exchange neighbour list is unavailable");
   if(Flags.do_aniso != 0 &&
      (gpuHamiltonian.kaniso.empty() || gpuHamiltonian.eaniso.empty() || gpuHamiltonian.taniso.empty()))
      return reject("the requested anisotropy data is unavailable");

   if(Flags.do_jtensor) {
      if(gpuHamiltonian.j_tensor.empty())
         return reject("the requested tensor-exchange data is unavailable");
      return true;
   }

   if(gpuHamiltonian.ncoup.empty())
      return reject("the exchange coupling data is unavailable");
   return true;
}

bool GpuHamiltonianCalculations::canUseMultiscaleDipole(const Flag Flags, const SimulationParameters SimParam,
                                                        const deviceHamiltonian& gpuHamiltonian) const {
   (void)Flags;
   (void)SimParam;
   (void)gpuHamiltonian;

   // Future backend selection:
   // - atom-to-cell and cell-to-atom interpolation data available
   // - coarse demag/dipole grid dimensions known
   // - additive field contribution can be accumulated before integration
   return false;
}

bool GpuHamiltonianCalculations::initiate(const Flag Flags, const SimulationParameters SimParam, deviceHamiltonian& gpuHamiltonian) {
   N = SimParam.N;   // Number of atoms
   NH = SimParam.NH;    // Number of reduced atoms
   mnn = SimParam.mnn;
   mnndm = SimParam.mnndm;
   redHam.redNeibourCount = gpuHamiltonian.aHam;
   external_field=gpuHamiltonian.extfield;
   do_ene = Flags.do_ene;
   do_j_tensor = false;
   do_dm = false;
   do_aniso = 0;
   convolution.release();
   dipoleConvolution.release();
   macroCellIndex = nullptr;
   macroMoments.Free();
   numMacro = 0;
   refreshMacroMoments = false;
   const bool gpuDipoleRequested = FortranData::gpu_dipole_mode && *FortranData::gpu_dipole_mode != 0;
   if(gpuDipoleRequested && FortranData::pme_num_macro && *FortranData::pme_num_macro > 0 &&
      FortranData::pme_macro_grid && FortranData::NA && *FortranData::NA > 0 &&
      !gpuHamiltonian.macro_cell_index.empty() && !gpuHamiltonian.macro_nlistsize.empty()) {
      numMacro = *FortranData::pme_num_macro;
      macroCellIndex = gpuHamiltonian.macro_cell_index.data();
      macroMoments.Allocate(3, numMacro, SimParam.M);
      refreshMacroMoments = true;
      GpuDipoleDescriptorInput input{};
      input.atomistic_grid = {SimParam.N1, SimParam.N2, SimParam.N3};
      input.macro_grid = {static_cast<std::size_t>(FortranData::pme_macro_grid[0]),
                          static_cast<std::size_t>(FortranData::pme_macro_grid[1]),
                          static_cast<std::size_t>(FortranData::pme_macro_grid[2])};
      input.basis = *FortranData::NA;
      input.ensembles = static_cast<unsigned int>(SimParam.M);
      input.boundaries = {{SimParam.BC1, SimParam.BC2, SimParam.BC3}};
      input.c1 = SimParam.C1;
      input.c2 = SimParam.C2;
      input.c3 = SimParam.C3;
      input.macro_centers = gpuHamiltonian.macro_center.data();
      input.macro_count = numMacro;
      input.alat = SimParam.alat;
      input.tolerance = SimParam.gpu_dipole_tol;
      GpuDipoleConvolutionDescriptor dipoleDescriptor{};
      if(!makeEwald3dFftDipoleDescriptor(input, dipoleDescriptor) ||
         !dipoleConvolution.initiate(dipoleDescriptor, parallel.getWorkStream())) {
         throw std::runtime_error("GPU EWALD3D_FFT descriptor is invalid");
      }
      const auto grid = dipoleConvolution.descriptor().activeGrid();
      const auto padded = dipoleConvolution.fftLayout().real_grid;
      const auto& layout = dipoleConvolution.fftLayout();
      std::printf("Gpu: EWALD3D_FFT geometry staged (%zu x %zu x %zu coarse cells, %u basis channel%s; FFT grid %zu x %zu x %zu, %zu half-spectrum points; tol %.3e, field prefactor %.16e T).\n",
                  grid.n1, grid.n2, grid.n3, dipoleDescriptor.basis,
                  dipoleDescriptor.basis == 1 ? "" : "s", padded.n1, padded.n2, padded.n3,
                  layout.spectral_cells, dipoleDescriptor.tolerance, dipoleDescriptor.field_prefactor);
      std::printf("Gpu: EWALD3D_FFT memory: %.3f MiB persistent, %.3f MiB peak including workspace and construction.\n",
                  static_cast<double>(GpuDipoleConvolution::estimatePersistentBytes(dipoleDescriptor)) / (1024.0 * 1024.0),
                  static_cast<double>(GpuDipoleConvolution::estimateBytes(dipoleDescriptor)) / (1024.0 * 1024.0));
      std::printf("Gpu: CV6 macrocell moment aggregation prepared (%u cells, %zu ensemble%s); dipole FFT pending.\n",
                  numMacro, SimParam.M, SimParam.M == 1 ? "" : "s");
   } else if(gpuDipoleRequested) {
      throw std::runtime_error("GPU EWALD3D_FFT requested without staged macrocell data");
   }
   // This is deliberately after descriptor creation: a requested valid mode
   // exercises the ownership, geometry, and memory contract, but never enters
   // the unfinished physical field/energy path.
   if(gpuDipoleRequested) {
      throw std::runtime_error(
         "GPU dipole mode was requested, but its field/energy operator is not yet available; "
         "use gpu_dipole_mode OFF until CV6 field validation is complete");
   }
   backend = {};
   backend.convolution_ready = canUseLatticeConvolution(Flags, SimParam, gpuHamiltonian);
   backend.multiscale_ready = canUseMultiscaleDipole(Flags, SimParam, gpuHamiltonian);
   if(backend.multiscale_ready) {
      backend.dipole = GpuHamiltonianBackend::MultiscaleDipole;
   }

   if(redHam.redNeibourCount.empty()) {
      initiated = false;
       //std::printf("HERE - 1\n");  
      return false;
   }
   if(N % 32 != 0) {
      std::printf("Note: Performance is better if the number of atoms is a multiple of 32.\n");
   }

   //------- Anisotropy -------//
   if(Flags.do_aniso != 0) {
      aniso.kaniso = gpuHamiltonian.kaniso;
      aniso.eaniso = gpuHamiltonian.eaniso;
      aniso.taniso = gpuHamiltonian.taniso;
      aniso.sb = gpuHamiltonian.sb;
      GpuHamiltonianCalculations::do_aniso = Flags.do_aniso;  

      if(aniso.kaniso.empty() || aniso.eaniso.empty() || aniso.taniso.empty()|| aniso.sb.empty()) {
         initiated = false;
         return false;
      }
   }

   //------- Tensorial Exchange -------//
   if(Flags.do_jtensor == 1) {
      GpuHamiltonianCalculations::do_j_tensor = true;
      tenEx.mnn = mnn;
      tenEx.neighbourCount = gpuHamiltonian.nlistsize;
      tenEx.neighbourPos = gpuHamiltonian.nlist;
      tenEx.tensor= gpuHamiltonian.j_tensor;


      // for(unsigned int site = 0; site < N; site++) {
      //	const unsigned int * myPos  = &(nlist.get_data())[site];
      //	const unsigned int   mySize = nlistsize.get_data()[site];
      //	std::printf(" %d ", myPos[0]);
      //	std::printf("| ");
      //	for (unsigned int i = 0; i < tenEx.mnn; i++)
      //	{
      //		std::printf(" %d ", myPos[i * N]);
      //	}
      //	std::printf("\n");
      // }

      if(!gpuHamiltonian.neighbourListsPrepared) {
         parallel.gpuSiteCall(SetupNeighbourListExchangeTensor(tenEx, redHam));
      }

      // Did we get the memory?
      if(tenEx.tensor.empty() || tenEx.neighbourCount.empty() || tenEx.neighbourPos.empty()) {
         initiated = false;
         return false;
     

      }
      if(backend.convolution_ready) {
         GpuLatticeConvolutionDescriptor conv_desc{};
         conv_desc.n1 = static_cast<unsigned int>(SimParam.N1);
         conv_desc.n2 = static_cast<unsigned int>(SimParam.N2);
         conv_desc.n3 = static_cast<unsigned int>(SimParam.N3);
         conv_desc.basis = static_cast<unsigned int>(SimParam.NA);
         conv_desc.ensembles = static_cast<unsigned int>(SimParam.M);
         conv_desc.bc1 = SimParam.BC1;
         conv_desc.bc2 = SimParam.BC2;
         conv_desc.bc3 = SimParam.BC3;
         conv_desc.c1 = SimParam.C1;
         conv_desc.c2 = SimParam.C2;
         conv_desc.c3 = SimParam.C3;
         conv_desc.basis_positions = SimParam.Bas;
         backend.convolution_ready = convolution.initiate(conv_desc, parallel.getWorkStream());
         if(backend.convolution_ready) {
            backend.convolution_kernel_ready =
               convolution.buildTensorKernel(tenEx.tensor, tenEx.neighbourPos, tenEx.mnn,
                                             parallel.getWorkStream());
            if(backend.convolution_kernel_ready) {
               backend.exchange = GpuHamiltonianBackend::LatticeConvolution;
               std::printf("Gpu: device lattice convolution active (do_gpu_convolution Y; %u x %u x %u cells, %u basis, %u ensemble%s).\n",
                           conv_desc.n1, conv_desc.n2, conv_desc.n3, conv_desc.basis, conv_desc.ensembles,
                           conv_desc.ensembles == 1 ? "" : "s");
            } else {
               std::printf("Gpu: device lattice convolution requested but disabled: tensor FFT kernel construction failed; using sparse Hamiltonian.\n");
            }
         } else {
            std::printf("Gpu: device lattice convolution requested but disabled: FFT plan or buffer initialization failed; using sparse Hamiltonian.\n");
         }
      }
      gpuHamiltonian.neighbourListsPrepared = true;

      // Flag
      initiated = true;
      return true;
   }
else{
   //------- Heisenberg Exchange -------//
   ex.mnn = mnn;  // Max number of neighbours
   ex.coupling = gpuHamiltonian.ncoup;
   ex.neighbourCount = gpuHamiltonian.nlistsize;
   ex.neighbourPos = gpuHamiltonian.nlist;
   
   // Did we get the memory?
   if(ex.coupling.empty() || ex.neighbourCount.empty() || ex.neighbourPos.empty()) {
      initiated = false;
      return false;
   }

   // List setup kernel call
   if(!gpuHamiltonian.neighbourListsPrepared) {
      parallel.gpuSiteCall(SetupNeighbourList(ex, redHam));
   }

   //------- DM Interaction -------//
   dm.mnn = 0;
   if(Flags.do_dm) {
      GpuHamiltonianCalculations::do_dm = true;

      dm.mnn = mnndm;  // Max number of DM neighbours  // I CHANGED THE INDEX FROM 0 TO 1!!!

      dm.interaction = gpuHamiltonian.dmvect;
      dm.neighbourCount = gpuHamiltonian.dmlistsize;
      dm.neighbourPos = gpuHamiltonian.dmlist;

      if(dm.interaction.empty() || dm.neighbourCount.empty() || dm.neighbourPos.empty()) {
         initiated = false;
     

         return false;
      }
      if(!gpuHamiltonian.neighbourListsPrepared) {
         parallel.gpuSiteCall(SetupNeighbourListDM(dm, redHam));
      }
   }

   if(backend.convolution_ready) {
      GpuLatticeConvolutionDescriptor conv_desc{};
      conv_desc.n1 = static_cast<unsigned int>(SimParam.N1);
      conv_desc.n2 = static_cast<unsigned int>(SimParam.N2);
      conv_desc.n3 = static_cast<unsigned int>(SimParam.N3);
      conv_desc.basis = static_cast<unsigned int>(SimParam.NA);
      conv_desc.ensembles = static_cast<unsigned int>(SimParam.M);
      conv_desc.bc1 = SimParam.BC1;
      conv_desc.bc2 = SimParam.BC2;
      conv_desc.bc3 = SimParam.BC3;
      conv_desc.c1 = SimParam.C1;
      conv_desc.c2 = SimParam.C2;
      conv_desc.c3 = SimParam.C3;
      conv_desc.basis_positions = SimParam.Bas;
      backend.convolution_ready = convolution.initiate(conv_desc, parallel.getWorkStream());
      if(backend.convolution_ready) {
         backend.convolution_kernel_ready =
            convolution.buildIsotropicDmKernel(ex.coupling, ex.neighbourPos, ex.mnn,
                                               dm.interaction, dm.neighbourPos, dm.mnn,
                                               do_dm, parallel.getWorkStream());
         if(backend.convolution_kernel_ready) {
            backend.exchange = GpuHamiltonianBackend::LatticeConvolution;
            backend.dmi = do_dm ? GpuHamiltonianBackend::LatticeConvolution : GpuHamiltonianBackend::DirectSparse;
            std::printf("Gpu: device lattice convolution active (do_gpu_convolution Y; %u x %u x %u cells, %u basis, %u ensemble%s).\n",
                        conv_desc.n1, conv_desc.n2, conv_desc.n3, conv_desc.basis, conv_desc.ensembles,
                        conv_desc.ensembles == 1 ? "" : "s");
         } else {
            std::printf("Gpu: device lattice convolution requested but disabled: isotropic/DM FFT kernel construction failed; using sparse Hamiltonian.\n");
         }
      } else {
         std::printf("Gpu: device lattice convolution requested but disabled: FFT plan or buffer initialization failed; using sparse Hamiltonian.\n");
      }
   }

   gpuHamiltonian.neighbourListsPrepared = true;

   // Flag
   initiated = true;
   return true;
}
}

void GpuHamiltonianCalculations::heisge(deviceLattice& gpuLattice, deviceEnergies& gpuEnergies,
                                        bool measure, bool includeAnisotropy) {
   // CV6.1: keep macro blocks coherent with the current device moments.  This
   // runs before any future macro-grid dipole field evaluation and is harmless
   // while the FFT backend remains disabled.
   if(refreshMacroMoments) {
      macroMoments.zeros_async(parallel.getWorkStream());
      parallel.gpuAtomSiteCall(UpdateMacroMoments(macroCellIndex, gpuLattice.emomM,
                                                  macroMoments, numMacro));
   }

   // Kernel call
   //null_energy<<<1,1>>>(gpuLattice.energy);
   // The FFT convolution backend serves both field-only and measuring steps
   // (CV2). On measure steps it reduces energyM directly from the convolved
   // field: the bilinear (exchange+DM) energy lands in the exchange column for
   // isotropic exchange or the tensor column for do_jtensor, matching the CPU
   // convention. On-site anisotropy and external are added in the unpack kernel.
   // Note: for isotropic exchange WITH DM the exchange/DM projection is not
   // separated (col 0 carries the combined bilinear energy, col 2 stays 0); the
   // total (col 5) is still correct. Run such systems as do_jtensor=1 for a
   // matching per-column split.
   if(backend.convolution_ready && backend.convolution_kernel_ready) {
      GpuLatticeConvolutionAnisotropy anis{};
      if(do_aniso != 0) {
         anis.kaniso = aniso.kaniso.data();
         anis.eaniso = aniso.eaniso.data();
         anis.taniso = aniso.taniso.data();
         anis.sb = aniso.sb.data();
         anis.present = true;
      }
      const unsigned int energy_col = do_j_tensor ? 3u : 0u;
      if(measure) {
         gpuEnergies.energyM.zeros_async(parallel.getWorkStream());
         convolution.apply(gpuLattice, external_field, anis, includeAnisotropy, &gpuEnergies,
                           energy_col, parallel.getWorkStream());
      } else {
         convolution.apply(gpuLattice, external_field, anis, includeAnisotropy, nullptr,
                           energy_col, parallel.getWorkStream());
      }
      return;
   }

   if(measure){
      //gpuEnergies.totalM.zeros();
      //gpuEnergies.extM.zeros();
      gpuEnergies.energyM.zeros_async(parallel.getWorkStream());
   }


    
   if(do_j_tensor) {
     //if(measure) gpuEnergies.tensorM.zeros();

      if(do_aniso != 0 && includeAnisotropy) {
         //if(measure) gpuEnergies.aniM.zeros();
         if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<false, true, true, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                          gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
         else parallel.gpuAtomSiteEnsembleCall(Heisge<false, true, true, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                          gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
      }
      else{
         if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<false, false, true, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                          gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
         else parallel.gpuAtomSiteEnsembleCall(Heisge<false, false, true, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                          gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
      }

   } 
   else {
     // if(measure) gpuEnergies.exchM.zeros();
      if(do_dm){
        //if(measure) gpuEnergies.dmM.zeros();
         if(do_aniso != 0 && includeAnisotropy){
            //if(measure) gpuEnergies.aniM.zeros();
            if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<true, true, false, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
            else parallel.gpuAtomSiteEnsembleCall(Heisge<true, true, false, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
            //printf("dmaniso\n");
         }
         else{
            if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<true, false, false, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                            gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
            else parallel.gpuAtomSiteEnsembleCall(Heisge<true, false, false, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                            gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
         
            //printf("dm\n");
                                          }
      }
      else{
         if(do_aniso != 0 && includeAnisotropy){
            //if(measure) gpuEnergies.aniM.zeros();
            if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<false, true, false, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
            else parallel.gpuAtomSiteEnsembleCall(Heisge<false, true, false, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
         
            //printf("aniso\n");
                                          }
         else{
            if(measure) parallel.gpuAtomSiteEnsembleCall(Heisge<false, false, false, true>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
            else parallel.gpuAtomSiteEnsembleCall(Heisge<false, false, false, false>(gpuLattice.beff, gpuLattice.eneff, gpuEnergies.energyM,
                                             gpuLattice.emomM, external_field, ex, dm, tenEx, aniso, redHam));
         
            //printf("plain\n");
                                          }
      }     
   }
   return;
}
