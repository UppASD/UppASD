/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Row-major linear index of an n-d coordinate, last dimension fastest-moving. */
size_t to_linear(const size_t ndims, const size_t* extents, const size_t* coordinates);

/** Row-major n-d coordinate of a linear index, last dimension fastest-moving. */
void to_spatial(const size_t index, const size_t ndims, const size_t* extents,
                size_t* coordinates);

#ifdef __cplusplus
}
#endif
