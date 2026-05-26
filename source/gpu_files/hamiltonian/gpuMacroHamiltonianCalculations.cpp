#pragma once
#include "c_headers.hpp"
#include "tensor.hpp"
#include "real_type.h"
#include "gpuStructures.hpp"
#include <numeric>
#include "gpu_wrappers.h"
#include "gpuMacroHamiltonianCalculations.hpp"

/*class GpuHamiltonianCalculations::HeisgeJijAniso : public ParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   real* eneff;
   //GpuTensor<real, 1> exchM;
   //GpuTensor<real, 1> aniM;
   //GpuTensor<real, 1> extM;
   //GpuTensor<real, 1> totalM;
   GpuTensor<real, 2> energyM;
   const real* coup;
   const unsigned int* pos;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const real* kaniso;
   const real* eaniso;
   const unsigned int* taniso;
   const real* sb;
   const unsigned int* aham;
   //int do_ene;
   bool measure;

public:
   HeisgeJijAniso(GpuTensor<real, 3>& p_beff, GpuTensor<real, 3>& p_eneff, GpuTensor<real, 2> p_energyM, const GpuTensor<real, 3>& p_emomM, const GpuTensor<real, 3>& p_ext_f, 
                  const Exchange& ex, const DMinteraction& dm, const Anisotropy& aniso, const HamRed& redHam, const int p_do_ene, const bool p_measure)
             : energyM(p_energyM)
    {
      beff = p_beff.data();
      eneff = p_eneff.data();
      emomM = p_emomM.data();
      ext_f = p_ext_f.data();

      coup = ex.coupling.data();
      pos = ex.neighbourPos.data();
      mnn = ex.mnn;

      kaniso = aniso.kaniso.data();
      eaniso = aniso.eaniso.data();
      taniso = aniso.taniso.data();
      sb = aniso.sb.data();

      aham = redHam.redNeibourCount.data();
      //do_ene = p_do_ene;
      measure = p_measure;

   }

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      // Field
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;
      real ax = (real)0.0;
      real ay = (real)0.0;
      real az = (real)0.0;
      real ax_en = (real)0.0;
      real ay_en = (real)0.0;
      real az_en = (real)0.0;
      real ex = (real)0.0;
      real ey = (real)0.0;
      real ez = (real)0.0;
      const unsigned int type = taniso[site];  // type of the anisotropy: 0 = none, 1 = uniaxial, 2 = cubic

      // Pointers with fixed indices
      const unsigned int rsite = aham[site] - 1;
      const real* site_coup = &coup[rsite];
      const unsigned int* site_pos = &pos[site];
      const real* my_emomM = &emomM[ensemble * N * 3];
      // Exchange term loop
      for(unsigned int i = 0; i < mnn; i++) {
         const auto x_offset = site_pos[i * N] * 3;
         const auto c = site_coup[i * NH];
         // printf("%f\n", c);
         x += c * my_emomM[x_offset + 0];
         y += c * my_emomM[x_offset + 1];
         z += c * my_emomM[x_offset + 2];
      }

      //Anisotropy
      const real Sx = emomM[atom * 3 + 0];
      const real Sy = emomM[atom * 3 + 1];
      const real Sz = emomM[atom * 3 + 2];

      // direction of uniaxial anisotropy
      ex = eaniso[0 + site * 3];
      ey = eaniso[1 + site * 3];
      ez = eaniso[2 + site * 3];

      // anisotropy constants
      const real k1 = kaniso[0 + site * 2];
      const real k2 = kaniso[1 + site * 2];

      if(type == 1 || type == 7)  // uniaxial anisotropy
      {
         const real tt1 = Sx * ex + Sy * ey + Sz * ez;
         const real tt2 = k1 + (real)2.0 * k2 * (1 - tt1 * tt1);
         const real tt3 = (real)2.0 * tt1 * tt2;

         const real tt2_en = k1 + k2 * ((real)2.0 - tt1 * tt1);
         const real tt3_en = tt1 * tt2_en;

         ax += -tt3 * ex;
         ay += -tt3 * ey;
         az += -tt3 * ez;

         ax_en += -tt3_en * ex; //To Anders: sign????
         ay_en += -tt3_en * ey;
         az_en += -tt3_en * ez;
      }
      if(type == 2 || type == 7) {  // cubic anisotropy

         real k1_cubic = k1;
         real k2_cubic = k2;

         if(type == 7) {  // Apply uniaxial and cubic anisotropy: The Cubic Anisotropy constant = Uniaxial
                          // constant x sb
            k1_cubic *= sb[site];
            k2_cubic *= sb[site];
         }

         ax += (real)2.0 * k1_cubic * Sx * (Sy * Sy + Sz * Sz) + (real)2.0 * k2_cubic * Sx * Sy * Sy * Sz * Sz;
         ay += (real)2.0 * k1_cubic * Sy * (Sz * Sz + Sx * Sx) + (real)2.0 * k2_cubic * Sy * Sz * Sz * Sx * Sx;
         az += (real)2.0 * k1_cubic * Sz * (Sx * Sx + Sy * Sy) + (real)2.0 * k2_cubic * Sz * Sx * Sx * Sy * Sy;

         ax_en += k1_cubic * Sx * (Sy * Sy + Sz * Sz)/((real)2.0) + k2_cubic * Sx * Sy * Sy * Sz * Sz/((real)3.0);//To Anders: sign???
         ay_en += k1_cubic * Sy * (Sz * Sz + Sx * Sx)/((real)2.0) + k2_cubic * Sy * Sz * Sz * Sx * Sx/((real)3.0);
         az_en += k1_cubic * Sz * (Sx * Sx + Sy * Sy)/((real)2.0) + k2_cubic * Sz * Sx * Sx * Sy * Sy/((real)3.0);
      }

      const real ext_x = ext_f[atom * 3 + 0];
      const real ext_y = ext_f[atom * 3 + 1];
      const real ext_z = ext_f[atom * 3 + 2];

      // Save field
      beff[atom * 3 + 0] = x + ax + ext_x;
      beff[atom * 3 + 1] = y + ay + ext_y;
      beff[atom * 3 + 2] = z + az + ext_z;

      eneff[atom * 3 + 0] = x + ax_en + ext_x;
      eneff[atom * 3 + 1] = y + ay_en + ext_y;
      eneff[atom * 3 + 2] = z + az_en + ext_z;
};*/
__global__ void Jij(GpuTensor<real, 3> beff, GpuTensor<real, 3> eneff, const GpuTensor<real, 3> emomM, const GpuTensor<real, 3>ext_f, 
                const GpuTensor<unsigned int, 1> nlistsize, const GpuTensor<unsigned int, 2> nlist, const GpuTensor<real, 2> ncoup, 
                const GpuTensor<unsigned int, 1> taniso, const GpuTensor<real, 2> eaniso, const GpuTensor<real, 2> kaniso, const GpuTensor<real, 1> sb,
                const GpuTensor<unsigned int, 1> cell_index, const GpuTensor<unsigned int, 1> macro_alistsize, const GpuTensor<unsigned int, 2> macro_atom_alist,
                const GpuTensor<unsigned int, 1> macro_halo_nlistsize, const GpuTensor<unsigned int, 2> macro_halo_to_global, const GpuTensor<unsigned int, 2> macro_atom_to_global,
                const GpuTensor<unsigned int, 2> macro_atom_local_nlistsize, const GpuTensor<unsigned int, 3> macro_atom_local_nlist, 
                const GpuTensor<unsigned int, 1> macro_cell_nlistsize, const GpuTensor<unsigned int, 2> macro_cell_nlist,
                const unsigned int max_macro_halo_size, const unsigned int max_num_atom_macro_cell)
{
    //TODO: ADD MACROCELLS <= blockIdx.y check
   // threads = {maxThreads, 1, 1};
   // blocks = {Num_macro, M, 1};
   unsigned int tid = threadIdx.x;
   unsigned int threadNum = blockDim.x;
   unsigned int mcInd = blockIdx.x; //macrocell index
   unsigned int mInd = blockIdx.y; //Mensemble index

   unsigned int unique_neigb_in_cell = macro_halo_nlistsize(mcInd);
   unsigned int atoms_in_cell = macro_alistsize(mcInd);
   unsigned int tasks = unique_neigb_in_cell * atoms_in_cell;

   //max_macro_halo_size     - max neighb
   //max_num_atom_macro_cell - max atoms
    //sh_neighb[3][max_macro_halo_size]; 
    //sh_beff[3][max_num_atom_macro_cell]; 

    extern __shared__ real shmem[];
    real* sh_neighb = shmem;
    real* sh_beff = sh_neighb + 3 * max_macro_halo_size;

    #define SH_N(j,i) sh_neighb[(i)*max_macro_halo_size + (j)]
    #define SH_B(j,i) sh_beff[(i)*max_num_atom_macro_cell + (j)]

    for(int i = 0; i < (3*(max_macro_halo_size + max_num_atom_macro_cell)); i+=threadNum){
        shmem[i] = 0;
    }

   for(unsigned int loc_nInd = 0; loc_nInd < unique_neigb_in_cell; loc_nInd+= threadNum){
        unsigned int global_nInd = macro_halo_to_global(loc_nInd, mcInd) - 1;

        SH_N(0, loc_nInd) = emomM(0, global_nInd, mInd);
        SH_N(1, loc_nInd) = emomM(1, global_nInd, mInd);
        SH_N(2, loc_nInd) = emomM(2, global_nInd, mInd);

    }

    __syncthreads();

    for(unsigned int i = 0; i < tasks; i+= threadNum){
        unsigned int loc_nInd = i % unique_neigb_in_cell; //neighbour index in halo
        unsigned int loc_aInd = i / unique_neigb_in_cell; //atom index in cell

        unsigned int global_nInd = macro_halo_to_global(loc_nInd, mcInd) - 1;
        unsigned int global_aInd = macro_atom_to_global(loc_nInd, mcInd) - 1;
        real c = ncoup(global_aInd, global_nInd);
        unsigned int nn = macro_atom_local_nlist(loc_nInd, loc_aInd, mcInd);
        
        atomicAdd(&SH_B(0, loc_aInd), c*SH_N(0, nn));
        atomicAdd(&SH_B(1, loc_aInd), c*SH_N(1, nn));
        atomicAdd(&SH_B(2, loc_aInd), c*SH_N(2, nn));

         /*for(unsigned int i = 0; i < mnn; i++) {
         const auto x_offset = site_pos[i * N] * 3;
         const auto c = site_coup[i * NH];
         // printf("%f\n", c);
         x += c * my_emomM[x_offset + 0];
         y += c * my_emomM[x_offset + 1];
         z += c * my_emomM[x_offset + 2];
        }*/

    }

       for(unsigned int loc_aInd = 0; loc_aInd < atoms_in_cell; loc_aInd+= threadNum){
        unsigned int global_aInd = macro_atom_to_global(loc_aInd, mcInd) - 1;

        beff(0, global_aInd, mInd) = SH_B(0, loc_aInd);
        beff(1, global_aInd, mInd) = SH_B(1, loc_aInd);
        beff(2, global_aInd, mInd) = SH_B(2, loc_aInd);

        //eneff(0, global_aInd, mInd) = SH_B(0, loc_aInd);
        //eneff(1, global_aInd, mInd) = SH_B(1, loc_aInd);
        //eneff(2, global_aInd, mInd) = SH_B(2, loc_aInd);
    
    }

}

__global__ void aniso(GpuTensor<real, 3> beff, GpuTensor<real, 3> eneff, const GpuTensor<real, 3> emomM, const GpuTensor<real, 3>ext_f, 
                const GpuTensor<unsigned int, 1> taniso, const GpuTensor<real, 2> eaniso, const GpuTensor<real, 2> kaniso, const GpuTensor<real, 1> sb,
                unsigned int N)
{
    unsigned int nInd = blockIdx.x * blockDim.x + threadIdx.x;
    if(nInd < N){
      unsigned int mInd = blockIdx.y; //Mensemble index
      real ax = (real)0.0;
      real ay = (real)0.0;
      real az = (real)0.0;
      real ax_en = (real)0.0;
      real ay_en = (real)0.0;
      real az_en = (real)0.0;
      real ex = (real)0.0;
      real ey = (real)0.0;
      real ez = (real)0.0;
      const unsigned int type = taniso(nInd);  // type of the anisotropy: 0 = none, 1 = uniaxial, 2 = cubic

      //Anisotropy
      const real Sx = emomM(0, nInd, mInd);
      const real Sy = emomM(1, nInd, mInd);
      const real Sz = emomM(2, nInd, mInd);

      // direction of uniaxial anisotropy
      ex = eaniso(0, nInd);
      ey = eaniso(1, nInd);
      ez = eaniso(2, nInd);

      // anisotropy constants
      const real k1 = kaniso(0, nInd);
      const real k2 = kaniso(1, nInd);

      if(type == 1 || type == 7)  // uniaxial anisotropy
      {
         const real tt1 = Sx * ex + Sy * ey + Sz * ez;
         const real tt2 = k1 + (real)2.0 * k2 * (1 - tt1 * tt1);
         const real tt3 = (real)2.0 * tt1 * tt2;

         const real tt2_en = k1 + k2 * ((real)2.0 - tt1 * tt1);
         const real tt3_en = tt1 * tt2_en;

         ax += -tt3 * ex;
         ay += -tt3 * ey;
         az += -tt3 * ez;

         ax_en += -tt3_en * ex; //To Anders: sign????
         ay_en += -tt3_en * ey;
         az_en += -tt3_en * ez;
      }
      if(type == 2 || type == 7) {  // cubic anisotropy

         real k1_cubic = k1;
         real k2_cubic = k2;

         if(type == 7) {  // Apply uniaxial and cubic anisotropy: The Cubic Anisotropy constant = Uniaxial
                          // constant x sb
            k1_cubic *= sb(nInd);
            k2_cubic *= sb(nInd);
         }

         ax += (real)2.0 * k1_cubic * Sx * (Sy * Sy + Sz * Sz) + (real)2.0 * k2_cubic * Sx * Sy * Sy * Sz * Sz;
         ay += (real)2.0 * k1_cubic * Sy * (Sz * Sz + Sx * Sx) + (real)2.0 * k2_cubic * Sy * Sz * Sz * Sx * Sx;
         az += (real)2.0 * k1_cubic * Sz * (Sx * Sx + Sy * Sy) + (real)2.0 * k2_cubic * Sz * Sx * Sx * Sy * Sy;

         ax_en += k1_cubic * Sx * (Sy * Sy + Sz * Sz)/((real)2.0) + k2_cubic * Sx * Sy * Sy * Sz * Sz/((real)3.0);//To Anders: sign???
         ay_en += k1_cubic * Sy * (Sz * Sz + Sx * Sx)/((real)2.0) + k2_cubic * Sy * Sz * Sz * Sx * Sx/((real)3.0);
         az_en += k1_cubic * Sz * (Sx * Sx + Sy * Sy)/((real)2.0) + k2_cubic * Sz * Sx * Sx * Sy * Sy/((real)3.0);
      }

      const real ext_x = ext_f(0, nInd, mInd);
      const real ext_y = ext_f(1, nInd, mInd);
      const real ext_z = ext_f(2, nInd, mInd);

      // Save field
      beff(0, nInd, mInd) += ax + ext_x;
      beff(1, nInd, mInd) += ay + ext_y;
      beff(2, nInd, mInd) += az + ext_z;

      eneff(0, nInd, mInd) += ax_en + ext_x;
      eneff(1, nInd, mInd) += ay_en + ext_y;
      eneff(2, nInd, mInd) += az_en + ext_z;
    }  
}

GpuMacroHamiltonianCalculations::GpuMacroHamiltonianCalculations(const Flag Flags, const SimulationParameters SimParam, const deviceHamiltonian& gpuHamiltonian, const deviceMacrocell& gpuMacro)
:   do_j_tensor(Flags.do_jtensor)
,   do_dm(Flags.do_dm)
,   do_aniso(Flags.do_aniso)
,   do_ene(Flags.do_ene)
,   N(SimParam.N)
,   M(SimParam.M)
,   NH(SimParam.NH)
,   mnn(SimParam.mnn)
,   mnndm(SimParam.mnndm)
,   nlistsize(gpuHamiltonian.nlistsize)
,   nlist(gpuHamiltonian.nlist)
,   ncoup(gpuHamiltonian.ncoup)
,   taniso(gpuHamiltonian.taniso)
,   eaniso(gpuHamiltonian.eaniso)
,   kaniso(gpuHamiltonian.kaniso)
,   sb(gpuHamiltonian.sb)
,   extfield(gpuHamiltonian.extfield)
,   macro(gpuMacro)
,   maxThreads_jij(256)
,   maxThreads_ani(512)
{
    threads_jij = {static_cast<uint32_t>(maxThreads_jij), 1, 1};
    blocks_jij = {static_cast<uint32_t>(macro.Num_macro), static_cast<uint32_t>(M), 1};

    threads_ani = {static_cast<uint32_t>(maxThreads_ani), 1, 1};
    blocks_ani = {static_cast<uint32_t>((N + maxThreads_ani - 1)/maxThreads_ani), static_cast<uint32_t>(M), 1};
}


// Destructor
GpuMacroHamiltonianCalculations::~GpuMacroHamiltonianCalculations(){}

void GpuMacroHamiltonianCalculations::calculate(deviceLattice& gpuLattice, deviceEnergies& gpuEnergies, bool measure) {
   if(do_j_tensor) {
      if(do_aniso != 0) {
        //tensor aniso
      }
      else {
        //Tensor aniso
      }
    } 

    else {
        if(do_dm) {
            if(do_aniso !=0) {
                // Jij DM aniso
            }
            else{
                //Jij DM
            }
        }
        else{
            if(do_aniso !=0) {
                //Jij aniso  
                size_t shmem_size = sizeof(float) * (3 * macro.max_halo_size + 3 * macro.max_num_atom);

                Jij<<<blocks_jij, threads_jij, shmem_size>>>(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, extfield, 
                                            nlistsize, nlist, ncoup, taniso, eaniso, kaniso, sb,
                                            macro.cell_index, macro.alistsize, macro.atom_alist,
                                            macro.halo_nlistsize, macro.halo_to_global, macro.atom_to_global,
                                            macro.atom_local_nlistsize, macro.atom_local_nlist, 
                                            macro.cell_nlistsize, macro.cell_nlist,
                                            macro.max_halo_size, macro.max_num_atom);
                aniso<<<blocks_ani, threads_ani, shmem_size>>>(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, extfield, 
                                            taniso, eaniso, kaniso, sb, N);
            }
            else {
                //Jji
                size_t shmem_size = sizeof(float) * (3 * macro.max_halo_size + 3 * macro.max_num_atom);

                Jij<<<blocks_jij, threads_jij, shmem_size>>>(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, extfield, 
                                            nlistsize, nlist, ncoup, taniso, eaniso, kaniso, sb,
                                            macro.cell_index, macro.alistsize, macro.atom_alist,
                                            macro.halo_nlistsize, macro.halo_to_global, macro.atom_to_global,
                                            macro.atom_local_nlistsize, macro.atom_local_nlist, 
                                            macro.cell_nlistsize, macro.cell_nlist,
                                            macro.max_halo_size, macro.max_num_atom);
            }     
        }
    }

}
 