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
#include "host-toolkit/memory.h"
#include "host-toolkit/xcorr.h"

int
main(void)
{
    const size_t ndims           = 2;
    const size_t field_extents[]  = {5, 6};
    const size_t kernel_extents[] = {3, 3};
    const size_t kernel_radius[]  = {kernel_extents[0] / 2, kernel_extents[1] / 2};
    const size_t domain_extents[] = {field_extents[0] - 2 * kernel_radius[0],
                                     field_extents[1] - 2 * kernel_radius[1]};

    const size_t field_count  = field_extents[0] * field_extents[1];
    const size_t kernel_count = kernel_extents[0] * kernel_extents[1];

    buffer_host_double in_field  = buffer_create_host_double(field_count);
    buffer_host_double out_field = buffer_create_host_double(field_count);
    buffer_host_double kernel    = buffer_create_host_double(kernel_count);

    for (size_t i = 0; i < field_count; ++i)
        in_field.data[i] = (double)i;

    for (size_t i = 0; i < kernel_count; ++i)
        kernel.data[i] = 1;

    xcorr(ndims, kernel_extents, kernel.data, field_extents, in_field.data, out_field.data);

    const real_t rtol = (real_t)kernel_count * REAL_EPSILON;
    const real_t atol = REAL_EPSILON;

    for (size_t dr = 0; dr < domain_extents[0]; ++dr) {
        for (size_t dc = 0; dc < domain_extents[1]; ++dc) {
            const size_t row = dr + kernel_radius[0];
            const size_t col = dc + kernel_radius[1];

            double expected = 0;
            for (size_t rr = row - 1; rr <= row + 1; ++rr)
                for (size_t cc = col - 1; cc <= col + 1; ++cc)
                    expected += (double)(rr * field_extents[1] + cc);

            const size_t out_coordinates[] = {row, col};
            const size_t out_index = to_linear(ndims, field_extents, out_coordinates);
            ERRCHK(isclose_real_t(out_field.data[out_index], expected, rtol, atol));
        }
    }

    buffer_destroy_host_double(&kernel);
    buffer_destroy_host_double(&out_field);
    buffer_destroy_host_double(&in_field);

    return EXIT_SUCCESS;
}
