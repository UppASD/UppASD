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

/** Returns the product of the count elements of data. */
size_t prod(const size_t count, const size_t* data);

/** Returns the sum of the count elements of data. */
size_t sum(const size_t count, const size_t* data);

void add_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out);
void add_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out);

void sub_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out);
void sub_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out);

void mul_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out);
void mul_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out);

void div_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out);
void div_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out);

void add_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out);
void add_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out);

void sub_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out);
void sub_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out);

void mul_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out);
void mul_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out);

void div_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out);
void div_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out);

/** Returns the dot product of the count elements of a and b. */
real_t dot_real_t(const size_t count, const real_t* a, const real_t* b);

/** Returns the Euclidean norm of the count elements of a, computed as sqrt(dot(a, a)). */
real_t norm_real_t(const size_t count, const real_t* a);

/** Writes a divided by its Euclidean norm to out. a must be nonzero. out may alias a. */
void normalize_real_t(const size_t count, const real_t* a, real_t* out);

/** Compares a and b for approximate equality, as with numpy.isclose: |a - b| <= atol + rtol * |b|. */
bool isclose_real_t(const real_t a, const real_t b, const real_t rtol, const real_t atol);

/** Compares a and b element-wise with isclose_real_t, as with numpy.allclose. */
bool allclose_real_t(const size_t count, const real_t* a, const real_t* b, const real_t rtol,
                     const real_t atol);

#ifdef __cplusplus
}
#endif
