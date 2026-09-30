#pragma once

#include "c_headers.hpp"
#include "gpu_wrappers.h"
#include "real_type.h"
#include "tensor.hpp"
#include <numeric>

namespace cg = cooperative_groups;
#ifndef M_PI
#define M_PI (3.14159265358979323846)å

#endif

__global__ void interpolate_atoms(const int interp_natoms, const int N, const int M,
                                  const int* indices, const int* first_neighbour,
                                  const int* neighbours, const real* weights, const real* mmom,
                                  real* emomM, real* emom, real* emom2);
