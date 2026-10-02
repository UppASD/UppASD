/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "common/debug.h"

#include <stdio.h>

void
print_int(const char* label, const int value)
{
    printf("%s: %d\n", label, value);
    fflush(stdout);
}

void
print_size_t(const char* label, const size_t value)
{
    printf("%s: %zu\n", label, value);
    fflush(stdout);
}

void
print_float(const char* label, const float value)
{
    printf("%s: %g\n", label, (double)value);
    fflush(stdout);
}

void
print_double(const char* label, const double value)
{
    printf("%s: %g\n", label, value);
    fflush(stdout);
}

void
print_double_array(const char* label, const size_t count, const double* arr)
{
    printf("%s = {\n", label);
    for (size_t i = 0; i < count; ++i)
        printf("\t%g,\n", arr[i]);
    printf("};\n");
    fflush(stdout);
}

void
print_float_array(const char* label, const size_t count, const float* arr)
{
    printf("%s = {\n", label);
    for (size_t i = 0; i < count; ++i)
        printf("\t%g,\n", (double)arr[i]);
    printf("};\n");
    fflush(stdout);
}
