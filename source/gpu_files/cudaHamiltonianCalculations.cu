#include <cuda_runtime.h>

#include "c_headers.hpp"
#include "cudaHamiltonianCalculations.hpp"
#include "cudaMatrix.hpp"
#include "hostMatrix.hpp"
#include "real_type.h"


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
   SetupNeighbourList(const Exchange& ex, const HamRed& redHam) {
      coup = ex.coupling;
      size = ex.neighbourCount;
      pos = ex.neighbourPos;
      aham = redHam.redNeibourCount;
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
   SetupNeighbourListExchangeTensor(const TensorialExchange& tenEx, const HamRed& redHam) {
      tensor = tenEx.tensor;
      size = tenEx.neighbourCount;
      pos = tenEx.neighbourPos;
      mnn = tenEx.mnn;
      aham = redHam.redNeibourCount;
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
   SetupAnisotropy(const Anisotropy& aniso) {
      kaniso = aniso.kaniso;
      eaniso = aniso.eaniso;
      taniso = aniso.taniso;
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
   SetupNeighbourListDM(const DMinteraction& dm, const HamRed& redHam) {
      coup = dm.interaction;
      size = dm.neighbourCount;
      pos = dm.neighbourPos;
      mnn = dm.mnn;
      aham = redHam.redNeibourCount;
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
   HeisgeJij(real* p1, const real* p2, const real* p3, const Exchange& ex, const DMinteraction& dm,
             const HamRed& redHam) {
      beff = p1;
      emomM = p2;
      ext_f = p3;

      coup = ex.coupling;
      pos = ex.neighbourPos;
      mnn = ex.mnn;

      dmcoup = dm.interaction;
      dmpos = dm.neighbourPos;
      dmmnn = dm.mnn;
      aham = redHam.redNeibourCount;
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
   }
};


class CudaHamiltonianCalculations::HeisJijTensor : public CudaParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   const real* tensor;
   const unsigned int* pos;
   const unsigned int* size;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const unsigned int* aham;

public:
   HeisJijTensor(real* p1, const real* p2, const real* p3, const TensorialExchange& tenEx,
                 const HamRed& redHam) {
      beff = p1;
      emomM = p2;
      ext_f = p3;

      tensor = tenEx.tensor;
      pos = tenEx.neighbourPos;
      size = tenEx.neighbourCount;
      mnn = tenEx.mnn;
      aham = redHam.redNeibourCount;
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
   }
};


class CudaHamiltonianCalculations::HeisgeJijAniso : public CudaParallelizationHelper::AtomSiteEnsemble {
private:
   real* beff;
   const real* emomM;
   const real* kaniso;
   const real* eaniso;
   const unsigned int* taniso;
   const real* sb;

public:
   HeisgeJijAniso(real* p1, const real* p2, const Anisotropy& aniso) {
      beff = p1;
      emomM = p2;
      kaniso = aniso.kaniso;
      eaniso = aniso.eaniso;
      taniso = aniso.taniso;
      sb = aniso.sb;
   }

   __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
      // Field
      real x = (real)0.0;
      real y = (real)0.0;
      real z = (real)0.0;
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

         x += -tt3 * ex;
         y += -tt3 * ey;
         z += -tt3 * ez;
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
      }

      // Save field
      beff[atom * 3 + 0] += x;
      beff[atom * 3 + 1] += y;
      beff[atom * 3 + 2] += z;
   }
};


class CudaHamiltonianCalculations::HeisgeJijElement
    : public CudaParallelizationHelper::ElementAxisSiteEnsemble {
private:
   real* beff;
   const real* coup;
   const unsigned int* pos;
   const unsigned int* size;
   const real* emomM;
   const real* ext_f;
   unsigned int mnn;
   const unsigned int* aham;

public:
   HeisgeJijElement(real* p1, const real* p5, const real* p6, const Exchange& ex, const HamRed& redHam) {
      beff = p1;
      coup = ex.coupling;
      pos = ex.neighbourPos;
      size = ex.neighbourCount;
      emomM = p5;
      ext_f = p6;
      mnn = ex.mnn;
      aham = redHam.redNeibourCount;
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
   }
};


////////////////////////////////////////////////////////////////////////////////
// Helpers
////////////////////////////////////////////////////////////////////////////////
template <typename T>
static void transpose(T* A, const T* B, std::size_t M, std::size_t N) {
   for(std::size_t y = 0; y < M; ++y) {
      for(std::size_t x = 0; x < N; ++x) {
         A[(x * M) + y] = B[(y * N) + x];
      }
   }
}


template <typename T, std::size_t I, std::size_t J, std::size_t K>
static void transpose(hostMatrix<T, 2, I, J, K>& A, const hostMatrix<T, 2, I, J, K>& B) {
   // Sizes
   std::size_t M = A.dimension_size(0);
   std::size_t N = A.dimension_size(1);

   if(B.dimension_size(1) != M || B.dimension_size(0) != N) {
      std::fprintf(stderr, "Error: illegal matrix transpose\n");
      std::exit(EXIT_FAILURE);
   }

   transpose(A.get_data(), B.get_data(), M, N);
}


// Function for testing time impact of optimal neighbour alignment
// Will not produce correct results
void alignOptimal(hostMatrix<unsigned int, 2>& nlist, bool same) {
   // Sizes
   std::size_t N = nlist.dimension_size(0);
   std::size_t mnn = nlist.dimension_size(1);

   for(std::size_t m = 0; m < mnn; ++m) {
      for(std::size_t n = 0; n < N; ++n) {
         nlist(n, m) = same ? ((m % N) + 1) : (((n + 32 * m) % N) + 1);
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////
CudaHamiltonianCalculations::CudaHamiltonianCalculations() : parallel(CudaParallelizationHelper::def) {
   initiated = false;
}


bool CudaHamiltonianCalculations::initiate(
    const hostMatrix<real, 2>& ncoup, const hostMatrix<unsigned int, 2>& nlist,
    const hostMatrix<unsigned int, 1>& nlistsize, const hostMatrix<real, 3, 3>& dm_ncoup,
    const hostMatrix<unsigned int, 2>& dm_nlist, const hostMatrix<unsigned int, 1>& dm_nlistsize,
    const int do_dm, const int do_j_tensor, const hostMatrix<real, 4, 3, 3> j_tensor, const int do_aniso,
    const hostMatrix<real, 2, 2> kaniso, const hostMatrix<real, 2, 3> eaniso,
    const hostMatrix<unsigned int, 1> taniso, const hostMatrix<real, 1> sb,
    const hostMatrix<unsigned int, 1>& aHam) {
   N = nlist.dimension_size(1);   // Number of atoms
   NH = ncoup.dimension_size(1);  // Number of reduced atoms
   redHam.redNeibourCount.clone(aHam);
   if(!redHam.redNeibourCount.has_data()) {
      release();
      return false;
   }
   if(N % 32 != 0) {
      std::printf("Note: Performance is better if the number of atoms is a multiple of 32.\n");
   }

   //------- Anisotropy -------//
   if(do_aniso != 0) {
      aniso.kaniso.clone(kaniso);
      aniso.eaniso.clone(eaniso);
      aniso.taniso.clone(taniso);
      aniso.sb.clone(sb);
      CudaHamiltonianCalculations::do_aniso = do_aniso;
      if(!aniso.kaniso.has_data() || !aniso.eaniso.has_data() || !aniso.taniso.has_data()
         || !aniso.sb.has_data()) {
         release();
         return false;
      }
   }

   //------- Tensorial Exchange -------//
   if(do_j_tensor == 1) {
      CudaHamiltonianCalculations::do_j_tensor = true;

      NH = j_tensor.dimension_size(3);

      // Matrixes are not transposed when using tensorial exchange
      // hostMatrix<real,4,3,3>         j_tensor_t;
      // hostMatrix<unsigned int,2> nlist_t;
      // j_tensor_t.initiate(3,3,N,tenEx.mnn);
      // nlist_t.initiate(N,tenEx.mnn);
      // transpose(j_tensor_t, j_tensor);
      // transpose(nlist_t, nlist);

      tenEx.mnn = j_tensor.dimension_size(2);
      tenEx.neighbourCount.clone(nlistsize);
      tenEx.neighbourPos.clone(nlist);
      tenEx.tensor.clone(j_tensor);

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
      if(!tenEx.tensor.has_data() || !tenEx.neighbourCount.has_data() || !tenEx.neighbourPos.has_data()) {
         release();
         return false;
      }
      // Flag
      initiated = true;
      return true;
   }

   //------- Heisenberg Exchange -------//
   ex.mnn = ncoup.dimension_size(0);  // Max number of neighbours

   // Transposing the matrices will make CUDA calculations faster
   hostMatrix<real, 2> ncoup_t;
   hostMatrix<unsigned int, 2> nlist_t;

   ncoup_t.initiate(NH, ex.mnn);
   nlist_t.initiate(N, ex.mnn);

   transpose(ncoup_t, ncoup);
   transpose(nlist_t, nlist);

   // TEST
   // alignOptimal(nlist_t, true);
   // std::printf("blubb: %f",ex.coupling);

   ex.coupling.clone(ncoup_t);
   ex.neighbourCount.clone(nlistsize);
   ex.neighbourPos.clone(nlist_t);

   // Did we get the memory?
   if(!ex.coupling.has_data() || !ex.neighbourCount.has_data() || !ex.neighbourPos.has_data()) {
      release();
      return false;
   }

   // List setup kernel call
   parallel.cudaSiteCall(SetupNeighbourList(ex, redHam));

   //------- DM Interaction -------//
   dm.mnn = 0;
   if(do_dm) {
      dm.mnn
          = dm_ncoup.dimension_size(1);  // Max number of DM neighbours  // I CHANGED THE INDEX FROM 0 TO 1!!!

      dm.interaction.clone(dm_ncoup);
      dm.neighbourCount.clone(dm_nlistsize);
      dm.neighbourPos.clone(dm_nlist);

      if(!dm.interaction.has_data() || !dm.neighbourCount.has_data() || !dm.neighbourPos.has_data()) {
         release();
         return false;
      }
      parallel.cudaSiteCall(SetupNeighbourListDM(dm, redHam));
   }

   // Flag
   initiated = true;
   return true;
}


void CudaHamiltonianCalculations::release() {
   ex.coupling.free();
   ex.neighbourCount.free();
   ex.neighbourPos.free();
   dm.interaction.free();
   dm.neighbourCount.free();
   dm.neighbourPos.free();
   tenEx.neighbourCount.free();
   tenEx.neighbourPos.free();
   tenEx.tensor.free();
   aniso.kaniso.free();
   aniso.eaniso.free();
   aniso.taniso.free();
   aniso.sb.free();
   redHam.redNeibourCount.free();
   initiated = false;
}


void CudaHamiltonianCalculations::heisge(cudaMatrix<real, 3, 3>& beff, const cudaMatrix<real, 3, 3>& emomM,
                                         const cudaMatrix<real, 3, 3>& external_field) {
   // Kernel call

   if(do_j_tensor == 1) {
      parallel.cudaAtomSiteEnsembleCall(HeisJijTensor(beff, emomM, external_field, tenEx, redHam));

   } else {
      parallel.cudaAtomSiteEnsembleCall(HeisgeJij(beff, emomM, external_field, ex, dm, redHam));
   }

   if(do_aniso != 0) {
      parallel.cudaAtomSiteEnsembleCall(HeisgeJijAniso(beff, emomM, aniso));
   }

   return;
}
