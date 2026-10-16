/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <float.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/debug.h"
#include "common/errhandler.h"
#include "common/stopwatch.h"
#include "device-toolkit/memory.h"
#include "device-toolkit/rand.h"
#include "host-toolkit/memory.h"

int
main(int argc, char* argv[])
{
    // Parse parameters
    size_t   count = 1024;
    uint64_t seed  = 1763249876ULL;
    for (int i = 1; i < argc - 1; ++i) {
        if (!strcmp(argv[i], "--count")) {
            count = strtoumax(argv[i++ + 1], NULL, 10);
        }
        else if (!strcmp(argv[i], "--seed")) {
            seed = strtoull(argv[i++ + 1], NULL, 10);
        }
        else {
            fprintf(stderr, "Usage: %s --count <count> --seed <seed>\n", argv[0]);
            ERRCHK(0);
        }
    }
    PRINT(count);
    PRINT(seed);

    stopwatch_t stopwatch;
    stopwatch_reset(&stopwatch);

    buffer_device_double dbuf = buffer_create_device_double(count);
    device_rng_t         drng = device_rng_create(seed, count);
    buffer_host_double   hbuf = buffer_create_host_double(count);

    buffer_destroy_host_double(&hbuf);
    device_rng_destroy(&drng);
    buffer_destroy_device_double(&dbuf);

    stopwatch_print_elapsed(stopwatch);

    printf("Complete\n");
    return EXIT_SUCCESS;
}
