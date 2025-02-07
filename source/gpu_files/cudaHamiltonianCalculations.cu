#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "cudaHamiltonianCalculations.hpp"
#include "tensor.cuh"
#include "real_type.h"
#include "cudaStructures.hpp"
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
class CudaHamiltonianCalculations::SetupNeighbourList : public CudaParallelizationHelper::Site {
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

// The neighbour list setup helper
//
// For Tensorial Exchange
class CudaHamiltonianCalculations::SetupNeighbourListExchangeTensor : public CudaParallelizationHelper::Site {
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
      return tensor[i + 3 * (j + 3 * (k + mnn * l))];
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
class CudaHamiltonianCalculations::SetupAnisotropy : public CudaParallelizationHelper::Site {
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
class CudaHamiltonianCalculations::SetupNeighbourListDM : public CudaParallelizationHelper::Site {
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
class CudaHamiltonianCalculations::HeisgeJij : public CudaParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   real* eneff;
   const real* coup;
   const unsigned int* pos;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const real* dmcoup;
   const unsigned int* dmpos;
   unsigned int dmmnn;
   const unsigned int* aham;

public:
   HeisgeJij(CudaTensor<real, 3>& p_beff, CudaTensor<real, 3>& p_eneff, const CudaTensor<real, 3>& p_emomM, const CudaTensor<real, 3>& p_ext_f, const Exchange& ex, const DMinteraction& dm,
             const HamRed& redHam) {
      beff = p_beff.data();
      eneff = p_eneff.data();
      emomM = p_emomM.data();
      ext_f = p_ext_f.data();

      coup = ex.coupling.data();
      pos = ex.neighbourPos.data();
      mnn = ex.mnn;

      dmcoup = dm.interaction.data();
      dmpos = dm.neighbourPos.data();
      dmmnn = dm.mnn;
      aham = redHam.redNeibourCount.data();
   }

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      // Field
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;

      // Pointers with fixed indices
      const unsigned int rsite = aham[site] - 1;
      const real* site_coup = &coup[rsite];
      const unsigned int* site_pos = &pos[site];
      const real* my_emomM = &emomM[ensemble * N * 3];
      // Exchange term loop
      for(unsigned int i = 0; i < mnn; i++) {
         unsigned int x_offset = site_pos[i * N] * 3;
         real c = site_coup[i * NH];
         // printf("%f\n", c);
         x += c * my_emomM[x_offset + 0];
         y += c * my_emomM[x_offset + 1];
         z += c * my_emomM[x_offset + 2];
      }

      // Phil's DM interaction implementation (still only incorporated into the isotropic Heisenberg exchange)
      for(unsigned int i = 0; i < dmmnn; i++) {
         unsigned int neighborPosIndex
             = dmpos[site * dmmnn + i];  // neighbor position in the site enemble given in 0,1,2,...,N-1

         unsigned int x_offset = neighborPosIndex * 3;

         real Sx = my_emomM[x_offset + 0];
         real Sy = my_emomM[x_offset + 1];
         real Sz = my_emomM[x_offset + 2];
         real Dx = dmcoup[0 + 3 * i + rsite * dmmnn * 3];
         real Dy = dmcoup[1 + 3 * i + rsite * dmmnn * 3];
         real Dz = dmcoup[2 + 3 * i + rsite * dmmnn * 3];

         x += -Dz * Sy + Dy * Sz;
         y += -Dx * Sz + Dz * Sx;
         z += -Dy * Sx + Dx * Sy;
      }

      // Save field
      beff[atom * 3 + 0] = x + ext_f[atom * 3 + 0];
      beff[atom * 3 + 1] = y + ext_f[atom * 3 + 1];
      beff[atom * 3 + 2] = z + ext_f[atom * 3 + 2];

      eneff[atom * 3 + 0] = x + ext_f[atom * 3 + 0];
      eneff[atom * 3 + 1] = y + ext_f[atom * 3 + 1];
      eneff[atom * 3 + 2] = z + ext_f[atom * 3 + 2];
   }
};

class CudaHamiltonianCalculations::HeisJijTensor : public CudaParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   real* eneff;
   const real* tensor;
   const unsigned int* pos;
   const unsigned int* size;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const unsigned int* aham;

public:
   HeisJijTensor(CudaTensor<real, 3>& p_beff, CudaTensor<real, 3>& p_eneff, const CudaTensor<real, 3>& p_emomM, const CudaTensor<real, 3>& p_ext_f, const TensorialExchange& tenEx,
                 const HamRed& redHam) {
      beff = p_beff.data();
      eneff = p_eneff.data();
      emomM = p_emomM.data();
      ext_f = p_ext_f.data();

      tensor = tenEx.tensor.data();
      pos = tenEx.neighbourPos.data();
      size = tenEx.neighbourCount.data();
      mnn = tenEx.mnn;
      aham = redHam.redNeibourCount.data();
   }

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      // Field
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;

      // Pointers with fixed indices
      const real* my_emomM = &emomM[ensemble * N * 3];
      const unsigned int rsite = aham[site] - 1;
      // emomM <--> (3,N,M)
      // tensor <---> (3,3,mnn,N)
      // pos   <--> (mnn,N)

      // Tensorial exchange coupling
      for(unsigned int i = 0; i < mnn; i++) {
         unsigned int neighborPosIndex
             = pos[site * mnn + i];  // neighbor position in the site enemble given in 0,1,2,...,N-1

         unsigned int x_offset = neighborPosIndex * 3;

         unsigned int k = i;
         unsigned int l = rsite;

         real J11 = tensor[0 + 3 * (0 + 3 * (k + mnn * l))];  // i=0,j=0
         real J12 = tensor[0 + 3 * (1 + 3 * (k + mnn * l))];  // i=0,j=1
         real J13 = tensor[0 + 3 * (2 + 3 * (k + mnn * l))];  // i=0,j=2
         real J21 = tensor[1 + 3 * (0 + 3 * (k + mnn * l))];  // i=1,j=0
         real J22 = tensor[1 + 3 * (1 + 3 * (k + mnn * l))];  // i=1,j=1
         real J23 = tensor[1 + 3 * (2 + 3 * (k + mnn * l))];  // i=1,j=2
         real J31 = tensor[2 + 3 * (0 + 3 * (k + mnn * l))];  // i=2,j=0
         real J32 = tensor[2 + 3 * (1 + 3 * (k + mnn * l))];  // i=2,j=1
         real J33 = tensor[2 + 3 * (2 + 3 * (k + mnn * l))];  // i=2,j=2

         // magnetic moment of current neighbor
         real Sx = my_emomM[x_offset + 0];
         real Sy = my_emomM[x_offset + 1];
         real Sz = my_emomM[x_offset + 2];

         x += J11 * Sx + J12 * Sy + J13 * Sz;
         y += J21 * Sx + J22 * Sy + J23 * Sz;
         z += J31 * Sx + J32 * Sy + J33 * Sz;
      }

      // Save field
      beff[atom * 3 + 0] = x + ext_f[atom * 3 + 0];
      beff[atom * 3 + 1] = y + ext_f[atom * 3 + 1];
      beff[atom * 3 + 2] = z + ext_f[atom * 3 + 2];

      eneff[atom * 3 + 0] = x + ext_f[atom * 3 + 0];
      eneff[atom * 3 + 1] = y + ext_f[atom * 3 + 1];
      eneff[atom * 3 + 2] = z + ext_f[atom * 3 + 2];
   }
};

class CudaHamiltonianCalculations::HeisgeJijAniso : public CudaParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   real* eneff;
   const real* emomM;
   const real* kaniso;
   const real* eaniso;
   const unsigned int* taniso;
   const real* sb;

public:
   HeisgeJijAniso(CudaTensor<real, 3>& p_beff, CudaTensor<real, 3>& p_eneff, const CudaTensor<real, 3>&  p_emomM, const Anisotropy& aniso) {
      beff = p_beff.data();
      eneff = p_eneff.data();
      emomM = p_emomM.data();
      kaniso = aniso.kaniso.data();
      eaniso = aniso.eaniso.data();
      taniso = aniso.taniso.data();
      sb = aniso.sb.data();
   }

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      // Field
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;

      real x_en = (real)0.0;
      real y_en = (real)0.0;
      real z_en = (real)0.0;
      // Magnetic moment at current site/atom
      real Sx = (real)0.0;
      real Sy = (real)0.0;
      real Sz = (real)0.0;
      // Uniaxial anisotropy unit vector
      real ex = (real)0.0;
      real ey = (real)0.0;
      real ez = (real)0.0;

      const unsigned int type = taniso[site];  // type of the anisotropy: 0 = none, 1 = uniaxial, 2 = cubic

      Sx = emomM[atom * 3 + 0];
      Sy = emomM[atom * 3 + 1];
      Sz = emomM[atom * 3 + 2];

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

         x += -tt3 * ex;
         y += -tt3 * ey;
         z += -tt3 * ez;

         x_en += -tt3_en * ex; //To Anders: sign????
         y_en += -tt3_en * ey;
         z_en += -tt3_en * ez;
      }
      if(type == 2 || type == 7) {  // cubic anisotropy

         real k1_cubic = k1;
         real k2_cubic = k2;

         if(type == 7) {  // Apply uniaxial and cubic anisotropy: The Cubic Anisotropy constant = Uniaxial
                          // constant x sb
            k1_cubic *= sb[site];
            k2_cubic *= sb[site];
         }

         x += (real)2.0 * k1_cubic * Sx * (Sy * Sy + Sz * Sz) + (real)2.0 * k2_cubic * Sx * Sy * Sy * Sz * Sz;
         y += (real)2.0 * k1_cubic * Sy * (Sz * Sz + Sx * Sx) + (real)2.0 * k2_cubic * Sy * Sz * Sz * Sx * Sx;
         z += (real)2.0 * k1_cubic * Sz * (Sx * Sx + Sy * Sy) + (real)2.0 * k2_cubic * Sz * Sx * Sx * Sy * Sy;

         x_en += k1_cubic * Sx * (Sy * Sy + Sz * Sz)/((real)2.0) + k2_cubic * Sx * Sy * Sy * Sz * Sz/((real)3.0);//To Anders: sign???
         y_en += k1_cubic * Sy * (Sz * Sz + Sx * Sx)/((real)2.0) + k2_cubic * Sy * Sz * Sz * Sx * Sx/((real)3.0);
         z_en += k1_cubic * Sz * (Sx * Sx + Sy * Sy)/((real)2.0) + k2_cubic * Sz * Sx * Sx * Sy * Sy/((real)3.0);
      }

      // Save field
      beff[atom * 3 + 0] += x;
      beff[atom * 3 + 1] += y;
      beff[atom * 3 + 2] += z;

      eneff[atom * 3 + 0] += x_en;
      eneff[atom * 3 + 1] += y_en;
      eneff[atom * 3 + 2] += z_en;
   }
};

class CudaHamiltonianCalculations::HeisgeJijElement
    : public CudaParallelizationHelper::ElementAxisSiteEnsemble {
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
   HeisgeJijElement(CudaTensor<real, 3>& p_beff,CudaTensor<real, 3>& p_eneff, const CudaTensor<real, 3>& p_emomM, const CudaTensor<real, 3>& p_ext_f, const Exchange& ex, const HamRed& redHam) {
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


////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

CudaHamiltonianCalculations::CudaHamiltonianCalculations() : parallel(CudaParallelizationHelper::def) {
   initiated = false;
}

bool CudaHamiltonianCalculations::initiate(const Flag Flags, const SimulationParameters SimParam, cudaHamiltonian& gpuHamiltonian) {
   N = SimParam.N;   // Number of atoms
   NH = SimParam.NH;    // Number of reduced atoms
   mnn = SimParam.mnn;
   mnndm = SimParam.mnndm;
   redHam.redNeibourCount = gpuHamiltonian.aHam;
   external_field=gpuHamiltonian.extfield;

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
      CudaHamiltonianCalculations::do_aniso = Flags.do_aniso;  

      if(aniso.kaniso.empty() || aniso.eaniso.empty() || aniso.taniso.empty()|| aniso.sb.empty()) {
         initiated = false;
         return false;
      }
   }

   //------- Tensorial Exchange -------//
   if(Flags.do_jtensor == 1) {
      CudaHamiltonianCalculations::do_j_tensor = true;
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

      parallel.cudaSiteCall(SetupNeighbourListExchangeTensor(tenEx, redHam));

      // Did we get the memory?
      if(tenEx.tensor.empty() || tenEx.neighbourCount.empty() || tenEx.neighbourPos.empty()) {
         initiated = false;
         return false;
     

      }
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
   parallel.cudaSiteCall(SetupNeighbourList(ex, redHam));

   //------- DM Interaction -------//
   dm.mnn = 0;
   if(Flags.do_dm) {
  

      dm.mnn = mnndm;  // Max number of DM neighbours  // I CHANGED THE INDEX FROM 0 TO 1!!!

      dm.interaction = gpuHamiltonian.dmvect;
      dm.neighbourCount = gpuHamiltonian.dmlistsize;
      dm.neighbourPos = gpuHamiltonian.dmlist;

      if(dm.interaction.empty() || dm.neighbourCount.empty() || dm.neighbourPos.empty()) {
         initiated = false;
     

         return false;
      }
      parallel.cudaSiteCall(SetupNeighbourListDM(dm, redHam));
   }

   // Flag
   initiated = true;
   return true;
}
}

void CudaHamiltonianCalculations::heisge(cudaLattice& gpuLattice) {
   // Kernel call

   if(do_j_tensor == 1) {
      parallel.cudaAtomSiteEnsembleCall(HeisJijTensor(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, external_field, tenEx, redHam));

   } else {
      parallel.cudaAtomSiteEnsembleCall(HeisgeJij(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, external_field, ex, dm, redHam));
   }

   if(do_aniso != 0) {
      parallel.cudaAtomSiteEnsembleCall(HeisgeJijAniso(gpuLattice.beff, gpuLattice.eneff, gpuLattice.emomM, aniso));
   }

   return;
}
