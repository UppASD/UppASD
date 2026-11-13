// SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
// SPDX-FileCopyrightText: 2026 UppASD contributors
//
// SPDX-License-Identifier: GPL-3.0-or-later

#include "device-implementation.h"

#include "common/errhandler.h"
#include "device-toolkit/errhandler.h"

static __device__ void
interpolate_atom(const int row, const size_t ensemble, const size_t atom_count,
                 const int* __restrict__ first_neighbour, const int* __restrict__ neighbours,
                 const double* __restrict__ weights, const double* __restrict__ emom, double v[3])
{
    // Fortran bounds firstNeighbour(row) .. firstNeighbour(row+1)-1 (1-based, inclusive)
    const int first = first_neighbour[row - 1] - 1;
    const int last  = first_neighbour[row] - 1;

    v[0] = v[1] = v[2] = 0;
    for (int i = first; i < last; ++i) {
        const size_t base = 3 * ((size_t)(neighbours[i] - 1) + ensemble * atom_count);
        v[0] += emom[base + 0] * weights[i];
        v[1] += emom[base + 1] * weights[i];
        v[2] += emom[base + 2] * weights[i];
    }

    // Pitfall: a zero sum yields NaN, as in UppASD atomInterpolation
    const double norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= norm;
    v[1] /= norm;
    v[2] /= norm;
}

static __global__ void
interpolate_interfaces(const size_t atom_count, const size_t ensemble_count,
                       const int* __restrict__ indices, const int* __restrict__ first_neighbour,
                       const int* __restrict__ neighbours, const double* __restrict__ weights,
                       const double* __restrict__ mmom, const double* __restrict__ emom_in,
                       const double* __restrict__ emom2_in, double* __restrict__ emom_out,
                       double* __restrict__ emom2_out, double* __restrict__ emomM_out)
{
    const size_t atom     = threadIdx.x + (size_t)blockIdx.x * blockDim.x;
    const size_t ensemble = threadIdx.y + (size_t)blockIdx.y * blockDim.y;
    if (atom >= atom_count || ensemble >= ensemble_count)
        return;

    const size_t moment = atom + ensemble * atom_count;
    const size_t base   = 3 * moment;
    const int    row    = indices[atom];
    if (row == 0) {
        for (size_t k = 0; k < 3; ++k) {
            emom_out[base + k]  = emom_in[base + k];
            emom2_out[base + k] = emom2_in[base + k];
        }
        return;
    }

    double v[3];
    double v2[3];
    interpolate_atom(row, ensemble, atom_count, first_neighbour, neighbours, weights, emom_in, v);
    interpolate_atom(row, ensemble, atom_count, first_neighbour, neighbours, weights, emom2_in, v2);

    for (size_t k = 0; k < 3; ++k) {
        emom_out[base + k]  = v[k];
        emom2_out[base + k] = v2[k];
        emomM_out[base + k] = v[k] * mmom[moment];
    }
}

void
device_multiscale_interpolate_interfaces(const size_t atom_count, const size_t ensemble_count,
                                         const int* indices, const int* first_neighbour,
                                         const int* neighbours, const double* weights,
                                         const double* mmom, const double* emom_in,
                                         const double* emom2_in, double* emom_out,
                                         double* emom2_out, double* emomM_out)
{
    ERRCHK(atom_count > 0);
    ERRCHK(ensemble_count > 0);
    ERRCHK(indices != nullptr);
    ERRCHK(first_neighbour != nullptr);
    ERRCHK(neighbours != nullptr);
    ERRCHK(weights != nullptr);
    ERRCHK(mmom != nullptr);
    ERRCHK(emom_in != nullptr);
    ERRCHK(emom2_in != nullptr);
    ERRCHK(emom_out != nullptr);
    ERRCHK(emom2_out != nullptr);
    ERRCHK(emomM_out != nullptr);
    ERRCHK(emom_in != emom_out);
    ERRCHK(emom2_in != emom2_out);

    // TODO: autotune the thread block dimensions
    const dim3 threads = {256, 1, 1};
    const dim3 blocks  = {
        (unsigned int)((atom_count + threads.x - 1) / threads.x),
        (unsigned int)((ensemble_count + threads.y - 1) / threads.y),
        1,
    };
    interpolate_interfaces<<<blocks, threads>>>(atom_count,
                                                ensemble_count,
                                                indices,
                                                first_neighbour,
                                                neighbours,
                                                weights,
                                                mmom,
                                                emom_in,
                                                emom2_in,
                                                emom_out,
                                                emom2_out,
                                                emomM_out);
    ERRCHK_CUDA_KERNEL();
}