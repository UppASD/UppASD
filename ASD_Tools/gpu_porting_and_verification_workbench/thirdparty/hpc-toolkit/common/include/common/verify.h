/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stddef.h>

#include "common/datatypes.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Maximum absolute and relative errors and their linear indices, computed in long double. */
typedef struct {
    long double max_abs_error;
    size_t      max_abs_error_index;
    long double max_rel_error;
    size_t      max_rel_error_index;
} verify_result_t;

/**
 * Compares the count elements of actual against expected.
 * The absolute error is |actual - expected|. The relative error is the absolute error divided by
 * max(|expected|, REAL_EPSILON). A NaN in either array is the maximum error, reported at the
 * last index where it occurs.
 */
verify_result_t verify_real_t(const size_t count, const real_t* expected, const real_t* actual);

/** Prints the maximum errors and their linear indices. */
void verify_print(const verify_result_t result);

#ifdef __cplusplus
}
#endif
