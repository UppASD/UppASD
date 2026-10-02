/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

// Template: no include guard, included once per type
#define HPC_T_NAME(name) HPC_CONCAT3(name, _, HPC_T)

/** Allocates a pinned host array of count elements. */
HPC_T* HPC_T_NAME(alloc_pinned)(const size_t count);

/** Frees the pinned host array at *ptr and nulls it. */
void HPC_T_NAME(free_pinned)(HPC_T** ptr);

/** Allocates a device array of count elements. */
HPC_T* HPC_T_NAME(alloc_device)(const size_t count);

/** Frees the device array at *ptr and nulls it. */
void HPC_T_NAME(free_device)(HPC_T** ptr);

/** Copies count elements from pinned host src to pinned host dst. */
void HPC_T_NAME(copy_pinned_h2h)(const size_t count, const HPC_T* src, HPC_T* dst);

/** Copies count elements from host src to device dst. */
void HPC_T_NAME(copy_h2d)(const size_t count, const HPC_T* src, HPC_T* dst);

/** Copies count elements from device src to host dst. */
void HPC_T_NAME(copy_d2h)(const size_t count, const HPC_T* src, HPC_T* dst);

/** Copies count elements from device src to device dst. */
void HPC_T_NAME(copy_d2d)(const size_t count, const HPC_T* src, HPC_T* dst);

/** Pinned host array of count elements. */
typedef struct {
    size_t count;
    HPC_T* data;
} HPC_T_NAME(buffer_pinned);

/** Device array of count elements. */
typedef struct {
    size_t count;
    HPC_T* data;
} HPC_T_NAME(buffer_device);

/** Creates a pinned host buffer of count elements. */
HPC_T_NAME(buffer_pinned) HPC_T_NAME(buffer_create_pinned)(const size_t count);

/** Frees the data of the pinned host buffer and sets its count to 0. */
void HPC_T_NAME(buffer_destroy_pinned)(HPC_T_NAME(buffer_pinned)* buffer);

/** Creates a device buffer of count elements. */
HPC_T_NAME(buffer_device) HPC_T_NAME(buffer_create_device)(const size_t count);

/** Frees the data of the device buffer and sets its count to 0. */
void HPC_T_NAME(buffer_destroy_device)(HPC_T_NAME(buffer_device)* buffer);

/** Copies pinned host buffer src to pinned host buffer dst. The counts must match. */
void HPC_T_NAME(buffer_copy_pinned_h2h)(const HPC_T_NAME(buffer_pinned)* src,
                                        HPC_T_NAME(buffer_pinned)* dst);

/** Copies pinned host buffer src to device buffer dst. The counts must match. */
void HPC_T_NAME(buffer_copy_h2d)(const HPC_T_NAME(buffer_pinned)* src,
                                 HPC_T_NAME(buffer_device)* dst);

/** Copies device buffer src to pinned host buffer dst. The counts must match. */
void HPC_T_NAME(buffer_copy_d2h)(const HPC_T_NAME(buffer_device)* src,
                                 HPC_T_NAME(buffer_pinned)* dst);

/** Copies device buffer src to device buffer dst. The counts must match. */
void HPC_T_NAME(buffer_copy_d2d)(const HPC_T_NAME(buffer_device)* src,
                                 HPC_T_NAME(buffer_device)* dst);

#undef HPC_T_NAME
