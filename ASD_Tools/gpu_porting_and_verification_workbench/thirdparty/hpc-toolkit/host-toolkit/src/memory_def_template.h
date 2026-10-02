/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

// Template: no include guard, included once per type
#define HPC_T_NAME(name) HPC_CONCAT3(name, _, HPC_T)

HPC_T*
HPC_T_NAME(alloc_host)(const size_t count)
{
    ERRCHK(count > 0);
    ERRCHK(count <= SIZE_MAX / sizeof(HPC_T));

    HPC_T* ptr = (HPC_T*)malloc(count * sizeof(ptr[0]));
    ERRCHK(ptr != nullptr);
    return ptr;
}

HPC_T*
HPC_T_NAME(realloc_host)(const size_t count, HPC_T** ptr)
{
    ERRCHK(count > 0);
    ERRCHK(count <= SIZE_MAX / sizeof(HPC_T));
    ERRCHK(ptr != nullptr);

    HPC_T* new_ptr = (HPC_T*)realloc(*ptr, count * sizeof(new_ptr[0]));
    ERRCHK(new_ptr != nullptr);

    if (new_ptr != nullptr)
        *ptr = nullptr;

    return new_ptr;
}

void
HPC_T_NAME(free_host)(HPC_T** ptr)
{
    ERRCHK(*ptr != nullptr);
    free(*ptr);
    *ptr = nullptr;
}

void
HPC_T_NAME(copy_h2h)(const size_t count, const HPC_T* src, HPC_T* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(count <= SIZE_MAX / sizeof(src[0]));

    const size_t    bytes = count * sizeof(src[0]);
    const uintptr_t s     = (uintptr_t)src;
    const uintptr_t d     = (uintptr_t)dst;
    ERRCHK(s + bytes <= d || d + bytes <= s);

    memcpy(dst, src, bytes);
}

HPC_T_NAME(buffer_host)
HPC_T_NAME(buffer_create_host)(const size_t count)
{
    return (HPC_T_NAME(buffer_host)){
        .count = count,
        .data  = HPC_T_NAME(alloc_host)(count),
    };
}

void
HPC_T_NAME(buffer_destroy_host)(HPC_T_NAME(buffer_host)* buffer)
{
    ERRCHK(buffer != nullptr);

    HPC_T_NAME(free_host)(&buffer->data);
    buffer->count = 0;
}

void
HPC_T_NAME(buffer_copy_h2h)(const HPC_T_NAME(buffer_host)* src, HPC_T_NAME(buffer_host)* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(src->count == dst->count);

    HPC_T_NAME(copy_h2h)(src->count, src->data, dst->data);
}

#undef HPC_T_NAME
