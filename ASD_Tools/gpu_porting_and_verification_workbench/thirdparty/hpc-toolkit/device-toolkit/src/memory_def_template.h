/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

// Template: no include guard, included once per type
#define HPC_T_NAME(name) HPC_CONCAT3(name, _, HPC_T)

HPC_T*
HPC_T_NAME(alloc_pinned)(const size_t count)
{
    ERRCHK(count > 0);
    ERRCHK(count <= SIZE_MAX / sizeof(HPC_T));

    HPC_T* ptr = nullptr;
    ERRCHK_CUDA(cudaMallocHost((void**)&ptr, count * sizeof(ptr[0])));
    ERRCHK(ptr != nullptr);
    return ptr;
}

void
HPC_T_NAME(free_pinned)(HPC_T** ptr)
{
    ERRCHK(*ptr != nullptr);
    ERRCHK_CUDA(cudaFreeHost(*ptr));
    *ptr = nullptr;
}

HPC_T*
HPC_T_NAME(alloc_device)(const size_t count)
{
    ERRCHK(count > 0);
    ERRCHK(count <= SIZE_MAX / sizeof(HPC_T));

    HPC_T* ptr = nullptr;
    ERRCHK_CUDA(cudaMalloc((void**)&ptr, count * sizeof(ptr[0])));
    ERRCHK(ptr != nullptr);
    return ptr;
}

void
HPC_T_NAME(free_device)(HPC_T** ptr)
{
    ERRCHK(*ptr != nullptr);
    ERRCHK_CUDA(cudaFree(*ptr));
    *ptr = nullptr;
}

void
HPC_T_NAME(copy_pinned_h2h)(const size_t count, const HPC_T* src, HPC_T* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(count <= SIZE_MAX / sizeof(src[0]));

    const size_t    bytes = count * sizeof(src[0]);
    const uintptr_t s     = (uintptr_t)src;
    const uintptr_t d     = (uintptr_t)dst;
    ERRCHK(s + bytes <= d || d + bytes <= s);

    ERRCHK_CUDA(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToHost));
}

void
HPC_T_NAME(copy_h2d)(const size_t count, const HPC_T* src, HPC_T* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(count <= SIZE_MAX / sizeof(src[0]));

    ERRCHK_CUDA(cudaMemcpy(dst, src, count * sizeof(src[0]), cudaMemcpyHostToDevice));
}

void
HPC_T_NAME(copy_d2h)(const size_t count, const HPC_T* src, HPC_T* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(count <= SIZE_MAX / sizeof(src[0]));

    ERRCHK_CUDA(cudaMemcpy(dst, src, count * sizeof(src[0]), cudaMemcpyDeviceToHost));
}

void
HPC_T_NAME(copy_d2d)(const size_t count, const HPC_T* src, HPC_T* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(count <= SIZE_MAX / sizeof(src[0]));

    const size_t    bytes = count * sizeof(src[0]);
    const uintptr_t s     = (uintptr_t)src;
    const uintptr_t d     = (uintptr_t)dst;
    ERRCHK(s + bytes <= d || d + bytes <= s);

    ERRCHK_CUDA(cudaMemcpy(dst, src, bytes, cudaMemcpyDeviceToDevice));
}

HPC_T_NAME(buffer_pinned)
HPC_T_NAME(buffer_create_pinned)(const size_t count)
{
    return (HPC_T_NAME(buffer_pinned)){
        .count = count,
        .data  = HPC_T_NAME(alloc_pinned)(count),
    };
}

void
HPC_T_NAME(buffer_destroy_pinned)(HPC_T_NAME(buffer_pinned)* buffer)
{
    ERRCHK(buffer != nullptr);

    HPC_T_NAME(free_pinned)(&buffer->data);
    buffer->count = 0;
}

HPC_T_NAME(buffer_device)
HPC_T_NAME(buffer_create_device)(const size_t count)
{
    return (HPC_T_NAME(buffer_device)){
        .count = count,
        .data  = HPC_T_NAME(alloc_device)(count),
    };
}

void
HPC_T_NAME(buffer_destroy_device)(HPC_T_NAME(buffer_device)* buffer)
{
    ERRCHK(buffer != nullptr);

    HPC_T_NAME(free_device)(&buffer->data);
    buffer->count = 0;
}

void
HPC_T_NAME(buffer_copy_pinned_h2h)(const HPC_T_NAME(buffer_pinned)* src,
                                   HPC_T_NAME(buffer_pinned)* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(src->count == dst->count);

    HPC_T_NAME(copy_pinned_h2h)(src->count, src->data, dst->data);
}

void
HPC_T_NAME(buffer_copy_h2d)(const HPC_T_NAME(buffer_pinned)* src, HPC_T_NAME(buffer_device)* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(src->count == dst->count);

    HPC_T_NAME(copy_h2d)(src->count, src->data, dst->data);
}

void
HPC_T_NAME(buffer_copy_d2h)(const HPC_T_NAME(buffer_device)* src, HPC_T_NAME(buffer_pinned)* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(src->count == dst->count);

    HPC_T_NAME(copy_d2h)(src->count, src->data, dst->data);
}

void
HPC_T_NAME(buffer_copy_d2d)(const HPC_T_NAME(buffer_device)* src, HPC_T_NAME(buffer_device)* dst)
{
    ERRCHK(src != nullptr);
    ERRCHK(dst != nullptr);
    ERRCHK(src->count == dst->count);

    HPC_T_NAME(copy_d2d)(src->count, src->data, dst->data);
}

#undef HPC_T_NAME
