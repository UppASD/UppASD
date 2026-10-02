/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "common/stopwatch.h"

#include <stdio.h>

void
stopwatch_reset(stopwatch_t* stopwatch)
{
    clock_gettime(CLOCK_MONOTONIC, &stopwatch->start);
}

int64_t
stopwatch_ns_elapsed(const stopwatch_t stopwatch)
{
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC, &now);
    return (int64_t)(now.tv_sec - stopwatch.start.tv_sec) * 1000000000LL +
           (now.tv_nsec - stopwatch.start.tv_nsec);
}

void
stopwatch_print_elapsed(const stopwatch_t stopwatch)
{
    const int64_t ns_elapsed = stopwatch_ns_elapsed(stopwatch);
    printf("Time elapsed: %g ms\n", (double)ns_elapsed / 1e6);
}
