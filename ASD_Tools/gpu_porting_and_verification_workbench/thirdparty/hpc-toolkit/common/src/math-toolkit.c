/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "common/math-toolkit.h"

#include <math.h>

#include "common/errhandler.h"

#define HPC_SQRT(x) _Generic((x), float: sqrtf, double: sqrt)(x)
#define HPC_FABS(x) _Generic((x), float: fabsf, double: fabs)(x)

size_t
prod(const size_t count, const size_t* data)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);

    size_t product = 1;
    for (size_t i = 0; i < count; ++i)
        product *= data[i];

    return product;
}

size_t
sum(const size_t count, const size_t* data)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);

    size_t total = 0;
    for (size_t i = 0; i < count; ++i)
        total += data[i];

    return total;
}

void
add_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] + b[i];
}

void
add_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] + b[i];
}

void
sub_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] - b[i];
}

void
sub_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i) {
        ERRCHK(a[i] >= b[i]);
        out[i] = a[i] - b[i];
    }
}

void
mul_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] * b[i];
}

void
mul_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] * b[i];
}

void
div_real_t(const size_t count, const real_t* a, const real_t* b, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] / b[i];
}

void
div_size_t(const size_t count, const size_t* a, const size_t* b, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i) {
        ERRCHK(b[i] != 0);
        out[i] = a[i] / b[i];
    }
}

void
add_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] + scalar;
}

void
add_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] + scalar;
}

void
sub_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] - scalar;
}

void
sub_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i) {
        ERRCHK(a[i] >= scalar);
        out[i] = a[i] - scalar;
    }
}

void
mul_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] * scalar;
}

void
mul_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] * scalar;
}

void
div_scalar_real_t(const size_t count, const real_t* a, const real_t scalar, real_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] / scalar;
}

void
div_scalar_size_t(const size_t count, const size_t* a, const size_t scalar, size_t* out)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(out != nullptr);
    ERRCHK(scalar != 0);

    for (size_t i = 0; i < count; ++i)
        out[i] = a[i] / scalar;
}

real_t
dot_real_t(const size_t count, const real_t* a, const real_t* b)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);

    real_t result = 0;
    for (size_t i = 0; i < count; ++i)
        result += a[i] * b[i];

    return result;
}

real_t
norm_real_t(const size_t count, const real_t* a)
{
    return HPC_SQRT(dot_real_t(count, a, a));
}

void
normalize_real_t(const size_t count, const real_t* a, real_t* out)
{
    ERRCHK(out != nullptr);

    const real_t norm = norm_real_t(count, a);
    ERRCHK(norm > 0);

    div_scalar_real_t(count, a, norm, out);
}

bool
isclose_real_t(const real_t a, const real_t b, const real_t rtol, const real_t atol)
{
    return HPC_FABS(a - b) <= atol + rtol * HPC_FABS(b);
}

bool
allclose_real_t(const size_t count, const real_t* a, const real_t* b, const real_t rtol,
                const real_t atol)
{
    ERRCHK(count > 0);
    ERRCHK(a != nullptr);
    ERRCHK(b != nullptr);

    for (size_t i = 0; i < count; ++i)
        if (!isclose_real_t(a[i], b[i], rtol, atol))
            return false;

    return true;
}
