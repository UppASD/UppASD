 #pragma once

#include "c_headers.hpp"
#include "gpu_wrappers.h"
#include "real_type.h"
#include "tensor.hpp"
#include <numeric>

namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// atomInterpolation
static __device__ void interpolate_atom(const int index, const int ensemble, const int N, const int* first_neighbour,
                 const int* neighbours, const GpuTensor<real, 3>  emom, const real* weights, real* v)
{
    const int lo = first_neighbour[index] - 1; // Inclusive
    const int hi = first_neighbour[index + 1] - 1; // Exclusive

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

// NOTE: not tested:
// - Whether correct indexing from Fortran to C
// - interp_atoms and N used correctly
//
// TODO
//
// 1) !!!Fix the race condition with emom/emom2 (also present in the model solution)
// 2) const __restrict__
// 3) autotuning thread block dimensions
// 4) Looking at memory ordering (one thread per vertex? Contiguous memory accesses)
// 5) looking at the profiler for next steps
//
// NOTE: initial draft. Not tested to run nor compile.
//
// multiscaleInterpolateInterfaces
__global__ void interpolate_atoms1(const int interp_natoms, const int N, const int M, const int* indices,
                  const int* first_neighbour, const int* neighbours, const real* weights,
                  const real* mmom, real* emomM, real* emom, const GpuTensor<real, 3> emom_buffer)
{
    const int atom = threadIdx.x + blockIdx.x * blockDim.x;
    if (atom >= interp_natoms)
        return;

    const int ensemble = threadIdx.y + blockIdx.y * blockDim.y;
    if (ensemble >= M)
        return;

    const int index = indices[atom] - 1;
    if (index == -1)
        return;

    const auto base = 3 * (atom + ensemble * N);

    real v[3];

    interpolate_atom(index, ensemble, N, first_neighbour, neighbours, emom_buffer, weights, v);
    emom[base + 0] = v[0];
    emom[base + 1] = v[1];
    emom[base + 2] = v[2];


    emomM[base + 0] = emom[base + 0] * mmom[atom + ensemble * interp_natoms];
    emomM[base + 1] = emom[base + 1] * mmom[atom + ensemble * interp_natoms];
    emomM[base + 2] = emom[base + 2] * mmom[atom + ensemble * interp_natoms];
}

__global__ void interpolate_atoms2(const int interp_natoms, const int N, const int M, const int* indices,
                  const int* first_neighbour, const int* neighbours, const real* weights,
                  real* emom2, const GpuTensor<real, 3> emom_buffer)
{
    const int atom = threadIdx.x + blockIdx.x * blockDim.x;
    if (atom >= interp_natoms)
        return;

    const int ensemble = threadIdx.y + blockIdx.y * blockDim.y;
    if (ensemble >= M)
        return;

    const int index = indices[atom] - 1;
    if (index == -1)
        return;

    const auto base = 3 * (atom + ensemble * N);

    real v[3];

    interpolate_atom(index, ensemble, N, first_neighbour, neighbours, emom_buffer, weights, v);
    emom2[base + 0] = v[0];
    emom2[base + 1] = v[1];
    emom2[base + 2] = v[2];

}

