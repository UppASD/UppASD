/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>

#include "common/datatypes.h"
#include "common/errhandler.h"
#include "common/math-toolkit.h"

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

int
main(void)
{
    const real_t rtol = (real_t)1e-5;
    const real_t atol = (real_t)1e-8;

    ERRCHK(isclose_real_t((real_t)1.0, (real_t)1.0, rtol, atol));
    ERRCHK(!isclose_real_t((real_t)1.0, (real_t)1.1, rtol, atol));

    /* Relative tolerance scales with the magnitude of b. */
    ERRCHK(isclose_real_t((real_t)1000.0, (real_t)1000.005, rtol, atol));
    ERRCHK(!isclose_real_t((real_t)1000.0, (real_t)1000.1, rtol, atol));

    /* Relative tolerance alone is not enough near zero; atol dominates. */
    ERRCHK(!isclose_real_t((real_t)1e-9, (real_t)0.0, (real_t)0.0, (real_t)0.0));
    ERRCHK(isclose_real_t((real_t)1e-9, (real_t)0.0, rtol, atol));

    const real_t a[]       = {(real_t)1.0, (real_t)1000.0, (real_t)1e-9};
    const real_t b[]       = {(real_t)1.0, (real_t)1000.005, (real_t)0.0};
    const real_t c[]       = {(real_t)1.0, (real_t)1000.1, (real_t)0.0};
    const size_t a_count   = ARRAY_SIZE(a);

    ERRCHK(allclose_real_t(a_count, a, b, rtol, atol));
    ERRCHK(!allclose_real_t(a_count, a, c, rtol, atol));

    return EXIT_SUCCESS;
}
