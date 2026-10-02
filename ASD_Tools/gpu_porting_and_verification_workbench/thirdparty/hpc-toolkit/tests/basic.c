/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>

#include "common/datatypes.h"
#include "common/errhandler.h"
#include "common/math-toolkit.h"
#include "host-toolkit/memory.h"

static void
fn(const size_t count, double* data)
{
    for (size_t i = 0; i < count; ++i)
        data[i] = (double)i;
}

static void
verify(const size_t count, const double* data)
{
    const real_t rtol = REAL_EPSILON;
    const real_t atol = REAL_EPSILON;

    for (size_t i = 0; i < count; ++i)
        ERRCHK(isclose_real_t(data[i], (double)i, rtol, atol));
}

int
main(void)
{
    const size_t count = 1024;

    buffer_host_double hbuf = buffer_create_host_double(count);

    fn(hbuf.count, hbuf.data);
    verify(hbuf.count, hbuf.data);

    buffer_destroy_host_double(&hbuf);
    return EXIT_SUCCESS;
}
