#pragma once

#include <cmath>
#include <cstdio>
#if defined(CUDA_V)
#include <curand.h>
#elif defined(HIP_V)
#include <hiprand/hiprand.h>
#endif
#include "c_headers.hpp"
#include "gpuInterpolation.hpp"
#include "gpuStructures.hpp"
#include "gpu_wrappers.h"
#include "multiscale_kernels.hpp"
#include "real_type.h"
#include "tensor.hpp"

/*
gpuLattice:

   GpuTensor<real, 3> emomM;
   GpuTensor<real, 3> emom;
   GpuTensor<real, 3>  emom2;
   GpuTensor<real, 2>  mmom;
   GpuTensor<real, 2>  mmom0;
   GpuTensor<real, 2>  mmom2;
   GpuTensor<real, 2>  mmomi;

deviceInterpolationInfo:
    unsigned int nrInterpAtoms; // Number of atoms affected
    GpuTensor<int, 1> indices;   //index on firstNeighbour corresponding to each atom. indices(i)
contains 0 if the atom is not affected by the interpolation. GpuTensor<int, 1>  firstNeighbour; //
index of the first neighbour in weights and neighbours GpuTensor<real, 1>  weights; // Per-neighbour
coefficient GpuTensor<int, 1>  neighbours; // Atom indices for neighbours participating in the
interpolation

*/

GpuInterpolation::GpuInterpolation(const unsigned int p_N, const unsigned int p_M,
                                   deviceInterpolationInfo& p_gpuInterpolationInfo,
                                   deviceLattice&           p_gpuLattice)
    : N(p_N), M(p_M), gpuInterpolationInfo(p_gpuInterpolationInfo), gpuLattice(p_gpuLattice)
{
    blocks  = {1, 1, 1};
    threads = {1, 1, 1};
}

GpuInterpolation::~GpuInterpolation() {}

#define ERRCHK(x)                                                                                  \
    do {                                                                                           \
        const auto retval{(x)};                                                                    \
        if (!retval) {                                                                             \
            fprintf(stderr, "%s:%d %s:%d %s\n", __FILE__, __LINE__, #x, retval, "ERRCHK failed");  \
        }                                                                                          \
    } while (0)

/*static void
interpolate_atom(const int index, const int ensemble, const int N, const int* first_neighbour,
                 const int* neighbours, const real* emom, const real* weights, real* v)
{
    const int lo = first_neighbour[index] - 1; // Inclusive
    const int hi = first_neighbour[index + 1]; // Exclusive

    v[0] = v[1] = v[2] = 0;
    for (int i = lo; i < hi; ++i) {
        const int base = 3 * (neighbours[i] - 1 + ensemble * N);
        v[0] += emom[base + 0] * weights[i];
        v[1] += emom[base + 1] * weights[i];
        v[2] += emom[base + 2] * weights[i];
    }

    const real norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;

    // NOTE: the else branch never entered in the model solution
}

static void
interpolate_atoms(const int interp_natoms, const int N, const int M, const int* indices,
                  const int* first_neighbour, const int* neighbours, const real* weights,
                  const real* mmom, real* emomM, real* emom, real* emom2)
{

    for (int atom = 0; atom < interp_natoms; ++atom) {
        const int index = indices[atom] - 1;

        if (index == -1)
            continue;

        for (int ensemble = 0; ensemble < M; ++ensemble) {
            const auto base = 3 * (atom + ensemble * N);

            real v[3];

            interpolate_atom(index, ensemble, N, first_neighbour, neighbours, emom, weights, v);
            emom[base + 0] = v[0];
            emom[base + 1] = v[1];
            emom[base + 2] = v[2];

            interpolate_atom(index, ensemble, N, first_neighbour, neighbours, emom2, weights, v);
            emom2[base + 0] = v[0];
            emom2[base + 1] = v[1];
            emom2[base + 2] = v[2];

            emomM[base + 0] = emom[base + 0] * mmom[atom + ensemble * N];
            emomM[base + 1] = emom[base + 1] * mmom[atom + ensemble * N];
            emomM[base + 2] = emom[base + 2] * mmom[atom + ensemble * N];
        }
    }
}*/

void
GpuInterpolation::interpolate()
{
    // Assert interp_natoms equal to ubound(interfaceInterpolation%indices, 1)
    const int interp_natoms = gpuInterpolationInfo.indices.extent(0);
    ERRCHK(M == gpuLattice.emom.extent(2));

    const dim3 threads = {256, 1, 1};
    const dim3 blocks  = {
        (interp_natoms + threads.x - 1) / threads.y,
        (M + threads.y - 1) / threads.y,
        1,
    };

    real*       emomM           = gpuLattice.emomM.data();
    real*       emom            = gpuLattice.emom.data();
    real*       emom2           = gpuLattice.emom2.data();
    const real* mmom            = gpuLattice.mmom.data();
    const int*  indices         = gpuInterpolationInfo.indices.data();
    const int*  first_neighbour = gpuInterpolationInfo.firstNeighbour.data();
    const real* weights         = gpuInterpolationInfo.weights.data();
    const int*  neighbours      = gpuInterpolationInfo.neighbours.data();
    ERRCHK(indices != nullptr); // Check associated(interp%indices). TODO confirm correct intent.

  /*  interpolate_atoms(interp_natoms,
                      N,
                      M,
                      indices,
                      first_neighbour,
                      neighbours,
                      weights,
                      mmom,
                      emomM,
                      emom,
                      emom2);*/

    interpolate_atoms<<<blocks, threads>>>(interp_natoms,
                                           N,
                                           M,
                                           indices,
                                           first_neighbour,
                                           neighbours,
                                           weights,
                                           mmom,
                                           emomM,
                                           emom,
                                           emom2);
}
