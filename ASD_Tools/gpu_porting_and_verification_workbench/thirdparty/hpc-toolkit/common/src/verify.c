/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "common/verify.h"

#include <math.h>
#include <stdio.h>

#include "common/errhandler.h"

#define HPC_FABS(x) _Generic((x), float: fabsf, double: fabs, long double: fabsl)(x)

verify_result_t
verify_real_t(const size_t count, const real_t* expected, const real_t* actual)
{
    ERRCHK(count > 0);
    ERRCHK(expected != nullptr);
    ERRCHK(actual != nullptr);

    verify_result_t result = {
        .max_abs_error       = 0,
        .max_abs_error_index = 0,
        .max_rel_error       = 0,
        .max_rel_error_index = 0,
    };

    const long double epsilon = (long double)REAL_EPSILON;

    for (size_t i = 0; i < count; ++i) {
        const long double e         = (long double)expected[i];
        const long double a         = (long double)actual[i];
        const long double magnitude = HPC_FABS(e);

        const long double abs_error = HPC_FABS(a - e);
        const long double rel_error = abs_error / (magnitude > epsilon ? magnitude : epsilon);

        if (isnan(abs_error) || abs_error > result.max_abs_error) {
            result.max_abs_error       = abs_error;
            result.max_abs_error_index = i;
        }
        if (isnan(rel_error) || rel_error > result.max_rel_error) {
            result.max_rel_error       = rel_error;
            result.max_rel_error_index = i;
        }
    }

    return result;
}

void
verify_print(const verify_result_t result)
{
    printf("Max absolute error: %Lg at index %zu\n", result.max_abs_error,
           result.max_abs_error_index);
    printf("Max relative error: %Lg at index %zu\n", result.max_rel_error,
           result.max_rel_error_index);
}
