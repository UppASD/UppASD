/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "host-toolkit/rand.h"

#include <stdint.h>
#include <stdio.h>

#include "common/errhandler.h"
#include "host-toolkit/memory.h"

host_rng_t
host_rng_create(const uint64_t seed, const size_t count)
{
    host_rng_t rng = {
        .count  = count,
        .states = alloc_host_uint64_t(count),
    };

    for (size_t idx = 0; idx < rng.count; ++idx)
        rng.states[idx] = seed + UINT64_C(123456) * (uint64_t)idx;

    return rng;
}

void
host_rng_destroy(host_rng_t* rng)
{
    ERRCHK(rng != nullptr);

    free_host_uint64_t(&rng->states);
    rng->count = 0;
}
