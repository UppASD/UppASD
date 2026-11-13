/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

// Template: no include guard, included once per type
#define HPC_T_NAME(name) HPC_CONCAT3(name, _, HPC_T)

/** Allocates a host array of count elements. */
HPC_T* HPC_T_NAME(alloc_host)(const size_t count);

/** Reallocates the host array at *ptr to count elements, nulling *ptr on success. */
HPC_T* HPC_T_NAME(realloc_host)(const size_t count, HPC_T** ptr);

/** Frees the host array at *ptr and nulls it. */
void HPC_T_NAME(free_host)(HPC_T** ptr);

/** Copies count elements from host src to host dst. */
void HPC_T_NAME(copy_h2h)(const size_t count, const HPC_T* src, HPC_T* dst);

/** Host array of count elements. */
typedef struct {
    size_t count;
    HPC_T* data;
} HPC_T_NAME(buffer_host);

/** Creates a host buffer of count elements. */
HPC_T_NAME(buffer_host) HPC_T_NAME(buffer_create_host)(const size_t count);

/** Frees the data of the host buffer and sets its count to 0. */
void HPC_T_NAME(buffer_destroy_host)(HPC_T_NAME(buffer_host)* buffer);

/** Copies host buffer src to host buffer dst. The counts must match. */
void HPC_T_NAME(buffer_copy_h2h)(const HPC_T_NAME(buffer_host)* src, HPC_T_NAME(buffer_host)* dst);

#undef HPC_T_NAME
