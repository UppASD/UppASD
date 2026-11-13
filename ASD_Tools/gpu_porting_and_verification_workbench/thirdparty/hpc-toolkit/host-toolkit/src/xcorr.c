/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "host-toolkit/xcorr.h"

#include "common/errhandler.h"
#include "common/index.h"
#include "common/math-toolkit.h"
#include "host-toolkit/memory.h"

void
xcorr(const size_t ndims, const size_t* kernel_extents, const real_t* kernel,
      const size_t* field_extents, const real_t* in_field, real_t* out_field)
{
    ERRCHK(ndims > 0);
    ERRCHK(kernel_extents != nullptr);
    ERRCHK(kernel != nullptr);
    ERRCHK(field_extents != nullptr);
    ERRCHK(in_field != nullptr);
    ERRCHK(out_field != nullptr);
    ERRCHK(in_field != out_field);

    for (size_t d = 0; d < ndims; ++d)
        ERRCHK(kernel_extents[d] % 2 == 1);

    const size_t kernel_count = prod(ndims, kernel_extents);

    buffer_host_size_t kernel_radius = buffer_create_host_size_t(ndims);
    div_scalar_size_t(ndims, kernel_extents, 2, kernel_radius.data);

    for (size_t d = 0; d < ndims; ++d)
        ERRCHK(field_extents[d] > 2 * kernel_radius.data[d]);

    buffer_host_size_t domain_extents = buffer_create_host_size_t(ndims);
    mul_scalar_size_t(ndims, kernel_radius.data, 2, domain_extents.data);
    sub_size_t(ndims, field_extents, domain_extents.data, domain_extents.data);

    const size_t domain_count = prod(ndims, domain_extents.data);

    buffer_host_size_t domain_coordinates = buffer_create_host_size_t(ndims);
    buffer_host_size_t out_coordinates    = buffer_create_host_size_t(ndims);
    buffer_host_size_t kernel_coordinates = buffer_create_host_size_t(ndims);
    buffer_host_size_t field_coordinates  = buffer_create_host_size_t(ndims);

    for (size_t i = 0; i < domain_count; ++i) {
        to_spatial(i, ndims, domain_extents.data, domain_coordinates.data);
        add_size_t(ndims, kernel_radius.data, domain_coordinates.data, out_coordinates.data);

        real_t accumulator = 0;
        for (size_t j = 0; j < kernel_count; ++j) {
            to_spatial(j, ndims, kernel_extents, kernel_coordinates.data);

            add_size_t(ndims, out_coordinates.data, kernel_coordinates.data,
                       field_coordinates.data);
            sub_size_t(ndims, field_coordinates.data, kernel_radius.data, field_coordinates.data);

            const size_t in_index = to_linear(ndims, field_extents, field_coordinates.data);
            accumulator += kernel[j] * in_field[in_index];
        }

        const size_t out_index = to_linear(ndims, field_extents, out_coordinates.data);
        out_field[out_index] = accumulator;
    }

    buffer_destroy_host_size_t(&field_coordinates);
    buffer_destroy_host_size_t(&kernel_coordinates);
    buffer_destroy_host_size_t(&out_coordinates);
    buffer_destroy_host_size_t(&domain_coordinates);
    buffer_destroy_host_size_t(&domain_extents);
    buffer_destroy_host_size_t(&kernel_radius);
}
