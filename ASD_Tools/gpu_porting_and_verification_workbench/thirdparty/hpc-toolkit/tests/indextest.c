/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdlib.h>

#include "common/datatypes.h"
#include "common/errhandler.h"
#include "common/index.h"
#include "common/math-toolkit.h"

int
main(void)
{
    {
        const size_t ndims          = 3;
        const size_t extents[]      = {1, 2, 3};
        size_t       coordinates[3];

        const size_t count = prod(ndims, extents);

        for (size_t a = 0; a < count; ++a) {
            to_spatial(a, ndims, extents, coordinates);
            const size_t c = to_linear(ndims, extents, coordinates);
            ERRCHK(a == c);
        }
    }

    {
        const size_t ndims          = 5;
        const size_t extents[]      = {5, 2, 7, 2, 4};
        size_t       coordinates[5];

        const size_t count = prod(ndims, extents);

        for (size_t a = 0; a < count; ++a) {
            to_spatial(a, ndims, extents, coordinates);
            const size_t c = to_linear(ndims, extents, coordinates);
            ERRCHK(a == c);
        }
    }

    return EXIT_SUCCESS;
}
