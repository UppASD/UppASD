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

typedef struct {
    size_t bytes;
    void*  data;
} device_reduce_workspace_t;

device_reduce_workspace_t device_reduce_workspace_create(const size_t bytes);
void                      device_reduce_workspace_destroy(device_reduce_workspace_t* workspace);

/**
 * Returns the workspace bytes device_reduce_sum needs for count elements of data.
 * result should not be modified. It is non-const due to the CUB API.
 */
size_t device_reduce_workspace_required_bytes(const size_t count, const real_t* data,
                                              real_t* result);

void device_reduce_sum(const size_t count, const real_t* data, real_t* result,
                       device_reduce_workspace_t* workspace);

#ifdef __cplusplus
}
#endif
