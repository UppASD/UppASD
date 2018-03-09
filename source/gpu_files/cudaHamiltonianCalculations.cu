#include <stdio.h>

#include <cuda_runtime.h>

using namespace std;

#include "real_type.h"
#include <iostream>

#include "hostMatrix.hpp"
#include "cudaMatrix.hpp"

#include "cudaHamiltonianCalculations.hpp"

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
class CudaHamiltonianCalculations::SetupNeighbourList :
	public CudaParallelizationHelper::Site {
private:
	real *               coup;
	unsigned int *       pos;
	const unsigned int * size;
	unsigned int         mnn;

public:
	SetupNeighbourList(const Exchange &ex) {
		coup  = ex.coupling;
		size  = ex.neighbourCount;
		pos   = ex.neighbourPos;
		mnn   = ex.mnn;
	}

        __device__ void each(unsigned int site) {
		real *         myCoup = &coup[site];
		unsigned int * myPos  = &pos[site];
		unsigned int   mySize = size[site];
		for (unsigned int i = 0; i < mnn; i++) {
			if (i < mySize)
				myPos[i * N]--;
			else {
				myCoup[i * N] = (real)0.0;
				myPos[i * N]  = 0;
			}
		}
	}
};


// Note (Thomas):
// For DM interaction
// Class sets everything between neighbours and maxneighbours
// to zero since hamiltonian implementation always runs to max neighbours
class CudaHamiltonianCalculations::SetupNeighbourListDM :
	public CudaParallelizationHelper::Site {
private:
	real *               coup;
	unsigned int *       pos;
	const unsigned int * size;
	unsigned int         mnn;

public:
	SetupNeighbourListDM(const DMinteraction & dm) {
		coup  = dm.interaction;
		size  = dm.neighbourCount;
		pos   = dm.neighbourPos;
		mnn   = dm.mnn;
	}

        __device__ void each(unsigned int site) {
		real *         myCoup = &coup[site * 3];
		unsigned int * myPos  = &pos[site];
		unsigned int   mySize = size[site];
		for (unsigned int i = 0; i < mnn; i++) {
			if (i < mySize)
				myPos[i * N]--;
			else {
				myCoup[i * N + 0] = (real)0.0;
				myCoup[i * N + 1] = (real)0.0;
				myCoup[i * N + 2] = (real)0.0;
				myPos[i * N]      = 0;
			}
		}
	}
};

// Note: (Thomas)
// Calculating the magnetic field from various effects
// such as the heisenberg field and DM interactions
// Added DM effect 2014/09/23
class CudaHamiltonianCalculations::HeisgeJij :
	public CudaParallelizationHelper::AtomSiteEnsemble {
private:
	real *               beff;
	const real *         coup;
	const unsigned int * pos;
	const real *         emomM;
	const real *         ext_f;
	unsigned int         mnn;
	const real *         dmcoup;
	const unsigned int * dmpos;
	unsigned int         dmmnn;
public:
	HeisgeJij(real * p1, const real * p2, const real * p3,
			 const Exchange & ex, const DMinteraction & dm) {
		beff   = p1;
		emomM  = p2;
		ext_f  = p3;

		coup   = ex.coupling;
		pos    = ex.neighbourPos;
		mnn    = ex.mnn;

		dmcoup = dm.interaction;
		dmpos  = dm.neighbourPos;
		dmmnn   = dm.mnn; 
	}

        __device__ void each(unsigned int atom, unsigned int site, unsigned int ensemble) {
		// Field
		real x = (real)0.0;
		real y = (real)0.0;
		real z = (real)0.0;

		// Pointers with fixed indices
		const real *         site_coup = &coup[site];
		const unsigned int * site_pos  = &pos[site];
		const real *         my_emomM  = &emomM[ensemble * N * 3];

		const real *         site_dmcoup    = &dmcoup[site];
		const unsigned int * site_dmpos     = &dmpos[site];

		// Exchange term loop
		for (unsigned int i = 0; i < mnn; i++) {
			unsigned int x_offset = site_pos[i * N] * 3; 
			real c = site_coup[i * N];
			x += c * my_emomM[x_offset + 0];
			y += c * my_emomM[x_offset + 1];
			z += c * my_emomM[x_offset + 2];
		}

		// DM interaction, almost no performance impact if dmmnn is 0	
		for (unsigned int i = 0; i < dmmnn; i++) {
			unsigned int x_offset = site_dmpos[i * N] * 3; 
			x += -site_dmcoup[i*N+2]*my_emomM[x_offset+1] + site_dmcoup[i*N+1]*my_emomM[x_offset+2];
			y += -site_dmcoup[i*N+0]*my_emomM[x_offset+2] + site_dmcoup[i*N+2]*my_emomM[x_offset+0];
			z += -site_dmcoup[i*N+1]*my_emomM[x_offset+0] + site_dmcoup[i*N+0]*my_emomM[x_offset+1];
		}

		// Save field
		beff[atom * 3 + 0] = x + ext_f[atom * 3 + 0];
		beff[atom * 3 + 1] = y + ext_f[atom * 3 + 1];
		beff[atom * 3 + 2] = z + ext_f[atom * 3 + 2];
	}
};


class CudaHamiltonianCalculations::HeisgeJijElement :
	public CudaParallelizationHelper::ElementAxisSiteEnsemble {
private:
	real *               beff;
	const real *         coup;
	const unsigned int * pos;
	const unsigned int * size;
	const real *         emomM;
	const real *         ext_f;
	unsigned int         mnn;
public:
	HeisgeJijElement(real * p1, const real * p5, const real * p6, const Exchange & ex) {
		beff   = p1;
		coup   = ex.coupling;
		pos    = ex.neighbourPos;
		size   = ex.neighbourCount;
		emomM  = p5;
		ext_f  = p6;
		mnn    = ex.mnn;
	}

        __device__ void each(unsigned int element, unsigned int axis, unsigned int site, unsigned int ensemble) {
		// Field
		real f = (real)0.0;

		// Pointers with fixed indices
		const real *         site_coup      = &coup[site];
		const unsigned int * site_pos       = &pos[site];
		const real *         ensemble_emomM = &emomM[ensemble * N * 3];

		// Exchange term loop
//		const unsigned int s = size[i];
//		for (int j = 0; j < s; j++) {
		for (unsigned int i = 0; i < mnn; i++) {
			unsigned int offset = site_pos[i * N] * 3;
			f += site_coup[i * N] * ensemble_emomM[offset + axis];
		}

		// Save field
		beff[element] = f + ext_f[element];
	}
};


////////////////////////////////////////////////////////////////////////////////
// Helpers
////////////////////////////////////////////////////////////////////////////////
template<typename T>
static void transpose(T * A, const T * B, size_t M, size_t N) {
	for (size_t y = 0; y < M; ++y)
		for (size_t x = 0; x < N; ++x)
			A[(x * M) + y] = B[(y * N) + x];
}

template <typename T, size_t I, size_t J, size_t K>
static void transpose(hostMatrix<T,2,I,J,K> &A, const hostMatrix<T,2,I,J,K> &B) {
	// Sizes
	size_t M = A.dimension_size(0);
	size_t N = A.dimension_size(1);

	if (B.dimension_size(1) != M || B.dimension_size(0) != N) {
		fprintf(stderr, "Error: illegal matrix transpose\n");
		exit(EXIT_FAILURE);
	}

	transpose(A.get_data(), B.get_data(), M, N);
}

// Function for testing time impact of optimal neighbour alignment
// Will not produce correct results
void alignOptimal(hostMatrix<unsigned int,2> &nlist, bool same) {
	// Sizes
	size_t N   = nlist.dimension_size(0);
	size_t mnn = nlist.dimension_size(1);

	for (size_t m = 0; m < mnn; ++m)
		for (size_t n = 0; n < N; ++n)
			nlist(n,m) = same ? ((m % N) + 1) : (((n + 32 * m) % N) + 1);
}


////////////////////////////////////////////////////////////////////////////////
// Class members
////////////////////////////////////////////////////////////////////////////////

CudaHamiltonianCalculations::CudaHamiltonianCalculations() :
	parallel(CudaParallelizationHelper::def) {
	initiated = false;
}

bool CudaHamiltonianCalculations::initiate(
		const hostMatrix<real,2>         &ncoup, 
		const hostMatrix<unsigned int,2> &nlist,
		const hostMatrix<unsigned int,1> &nlistsize,
		const hostMatrix<real,3,3>       &dm_ncoup, 
		const hostMatrix<unsigned int,2> &dm_nlist,
		const hostMatrix<unsigned int,1> &dm_nlistsize,
		const int 			 do_dm) {

	// Memory access is better if N is multiple of 32
	// (alignment of 128 bytes, see Cuda Best Parctice Guide)
	N      = ncoup.dimension_size(1);    // Number of atoms
	if (N % 32 != 0) {
		printf("Note: Performance is better if the number of atoms is a multiple of 32.\n");
	}

	//------- Heisenberg Exchange -------//
	ex.mnn    = ncoup.dimension_size(0);    // Max number of neighbours

	// Transposing the matrices will make CUDA calculations faster
	hostMatrix<real,2>         ncoup_t;
	hostMatrix<unsigned int,2> nlist_t;

	ncoup_t.initiate(N,ex.mnn);
	nlist_t.initiate(N,ex.mnn);

	transpose(ncoup_t, ncoup);
	transpose(nlist_t, nlist);

// TEST
//alignOptimal(nlist_t, true);
	//printf("blubb: %f",ex.coupling);

	ex.coupling.clone(ncoup_t);
	ex.neighbourCount.clone(nlistsize);
	ex.neighbourPos.clone(nlist_t);

	// Did we get the memory?
	if (!ex.coupling.has_data()       ||
	    !ex.neighbourCount.has_data() ||
	    !ex.neighbourPos.has_data()) {
		release();
		return false;
	}

	// List setup kernel call
	parallel.cudaSiteCall(SetupNeighbourList(ex));

	//------- DM Interaction -------//
	dm.mnn = 0;
	if (do_dm) {
		dm.mnn = dm_ncoup.dimension_size(0); // Max number of DM neighbours

		dm.interaction.clone(dm_ncoup);
		dm.neighbourCount.clone(dm_nlistsize);
		dm.neighbourPos.clone(dm_nlist);
		
		if (!dm.interaction.has_data()       ||
			!dm.neighbourCount.has_data() ||
			!dm.neighbourPos.has_data()) {
			release();
			return false;
		}
		parallel.cudaSiteCall(SetupNeighbourListDM(dm));
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
	initiated = false;
}

void CudaHamiltonianCalculations::heisge(cudaMatrix<real,3,3> &beff, 
		const cudaMatrix<real,3,3> &emomM,
		const cudaMatrix<real,3,3> &external_field) {
	// Kernel call
	parallel.cudaAtomSiteEnsembleCall(HeisgeJij(beff, emomM, external_field, ex, dm));

	//parallel.cudaElementAxisSiteEnsembleCall(HeisgeJijElement(beff, emomM, external_field, ex));
}

