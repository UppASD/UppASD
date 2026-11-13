/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "common/index.h"

#include "common/errhandler.h"
#include "common/math-toolkit.h"

size_t
to_linear(const size_t ndims, const size_t* extents, const size_t* coordinates)
{
    ERRCHK(ndims > 0);
    ERRCHK(extents != nullptr);
    ERRCHK(coordinates != nullptr);

    size_t index  = 0;
    size_t stride = 1;
    for (size_t d = ndims - 1; d < ndims; --d) {
        ERRCHK(coordinates[d] < extents[d]);
        index += coordinates[d] * stride;
        stride *= extents[d];
    }

    return index;
}

void
to_spatial(const size_t index, const size_t ndims, const size_t* extents, size_t* coordinates)
{
    ERRCHK(ndims > 0);
    ERRCHK(extents != nullptr);
    ERRCHK(coordinates != nullptr);

    ERRCHK(index < prod(ndims, extents));

    size_t remainder = index;
    for (size_t d = ndims - 1; d < ndims; --d) {
        coordinates[d] = remainder % extents[d];
        remainder /= extents[d];
    }
}
