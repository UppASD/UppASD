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

void print_int(const char* label, const int value);
void print_size_t(const char* label, const size_t value);
void print_float(const char* label, const float value);
void print_double(const char* label, const double value);
void print_double_array(const char* label, const size_t count, const double* arr);
void print_float_array(const char* label, const size_t count, const float* arr);

#ifdef __cplusplus
}
#endif

#define PRINT(x) _Generic((x), \
        int: print_int, \
        size_t: print_size_t, \
        float: print_float, \
        double: print_double \
    )(#x, (x))

#define PRINT_ARRAY(count, arr) _Generic((arr), \
        double*: print_double_array, \
        const double*: print_double_array, \
        float*: print_float_array, \
        const float*: print_float_array \
    )(#arr, (count), (arr))
