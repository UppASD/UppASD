/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "memory.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t    count;
    uint64_t* states;
} host_rng_t;

host_rng_t host_rng_create(const uint64_t seed, const size_t count);
void       host_rng_destroy(host_rng_t* rng);

#ifdef __cplusplus
}
#endif
