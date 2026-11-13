/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <math.h>
#include <stdlib.h>

#include "common/datatypes.h"
#include "common/errhandler.h"
#include "common/math-toolkit.h"
#include "common/verify.h"

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

int
main(void)
{
    const real_t rtol = REAL_EPSILON;
    const real_t atol = REAL_EPSILON;

    {
        const real_t expected[] = {1, 2, 4, 8};
        const real_t actual[]   = {1, 2, 4, 8};

        const verify_result_t result = verify_real_t(ARRAY_SIZE(expected), expected, actual);
        ERRCHK(result.max_abs_error <= 0);
        ERRCHK(result.max_rel_error <= 0);
    }

    {
        const real_t expected[] = {1, 2, 4, 8};
        const real_t actual[]   = {1, (real_t)2.5, 4, 9};

        const verify_result_t result = verify_real_t(ARRAY_SIZE(expected), expected, actual);
        verify_print(result);

        ERRCHK(result.max_abs_error_index == 3);
        ERRCHK(isclose_real_t((real_t)result.max_abs_error, 1, rtol, atol));
        ERRCHK(result.max_rel_error_index == 1);
        ERRCHK(isclose_real_t((real_t)result.max_rel_error, (real_t)0.25, rtol, atol));
    }

    {
        const real_t expected[] = {0, 1};
        const real_t actual[]   = {(real_t)0.5, 1};

        const verify_result_t result = verify_real_t(ARRAY_SIZE(expected), expected, actual);

        ERRCHK(result.max_rel_error_index == 0);
        ERRCHK(isfinite(result.max_rel_error));
        ERRCHK(isclose_real_t((real_t)result.max_rel_error, (real_t)0.5 / REAL_EPSILON, rtol,
                              atol));
    }

    {
        const real_t expected[] = {1, 2, 3};
        const real_t actual[]   = {1, (real_t)NAN, 4};

        const verify_result_t result = verify_real_t(ARRAY_SIZE(expected), expected, actual);

        ERRCHK(isnan(result.max_abs_error));
        ERRCHK(result.max_abs_error_index == 1);
        ERRCHK(isnan(result.max_rel_error));
        ERRCHK(result.max_rel_error_index == 1);
    }

    return EXIT_SUCCESS;
}
