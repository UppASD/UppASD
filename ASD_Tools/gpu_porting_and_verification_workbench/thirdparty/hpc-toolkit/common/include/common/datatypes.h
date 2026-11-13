/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <float.h>

#include "common/macros.h"

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(HPC_REAL_TYPE)
#define HPC_REAL_TYPE double
#endif

typedef HPC_REAL_TYPE real_t;

#define HPC_REAL_EPSILON_float  FLT_EPSILON
#define HPC_REAL_EPSILON_double DBL_EPSILON

/** Machine epsilon of real_t. */
#define REAL_EPSILON HPC_CONCAT2(HPC_REAL_EPSILON_, HPC_REAL_TYPE)

typedef enum [[nodiscard]] {
    HPC_SUCCESS = 0,
    HPC_FAILURE,
} hpc_status_t;

#ifdef __cplusplus
}
#endif
