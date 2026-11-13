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

/**
 * Applies an n-d cross-correlation kernel to the interior of a field.
 *
 * \f[
 *     \mathrm{out\_field}(x) = \sum_{s} \mathrm{kernel}(s) \,
 *     \mathrm{in\_field}(x + s)
 * \f]
 *
 * `in_field` and `out_field` must each hold `prod(field_extents)` elements,
 * matching the halo already present in `in_field`. Only the interior domain
 * (`field_extents[d] - 2 * (kernel_extents[d] / 2)` per dimension) is written;
 * the halo of `out_field` is left untouched. `in_field` and `out_field` must
 * not alias, since each output point reads a neighborhood of the input.
 */
void xcorr(const size_t ndims, const size_t* kernel_extents, const real_t* kernel,
           const size_t* field_extents, const real_t* in_field, real_t* out_field);

#ifdef __cplusplus
}
#endif
