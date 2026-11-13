/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stdint.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    struct timespec start;
} stopwatch_t;

void stopwatch_reset(stopwatch_t* stopwatch);

int64_t stopwatch_ns_elapsed(const stopwatch_t stopwatch);

void stopwatch_print_elapsed(const stopwatch_t stopwatch);

#ifdef __cplusplus
}
#endif
