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

void reduce_sum(const size_t count, const real_t* data, real_t* result);

#ifdef __cplusplus
}
#endif
