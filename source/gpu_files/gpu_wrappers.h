#pragma once

#if defined(HIP_V)
  #include <hip/hip_runtime.h>
  #include <hiprand/hiprand.h>
  #define GPU_DEVICE_RESET() hipDeviceReset()
  #define GPU_DEVICE_SYNCHRONIZE() hipDeviceSynchronize()
  #define GPU_ERROR_T hipError_t
  #define GPU_EVENT_T hipEvent_t
  #define GPU_EVENT_CREATE(event) hipEventCreate(event)
  #define GPU_EVENT_RECORD(event, stream) hipEventRecord(event, stream)
  #define GPU_FREE(ptr) hipFree(ptr)
  #define GPU_FREE_HOST(ptr) hipHostFree(ptr)
  #define GPU_GET_ERROR_STRING(err) hipGetErrorString(err)
  #define GPU_GET_LAST_ERROR() hipGetLastError()
  #define GPU_MALLOC(ptr, size) hipMalloc(ptr, size)
  #define GPU_MALLOC_HOST(ptr, size) hipHostMalloc(ptr, size)
  #define GPU_MEMCPY(dst, src, size, kind) hipMemcpy(dst, src, size, kind)
  #define GPU_MEMCPY_ASYNC(dst, src, size, kind, stream) hipMemcpyAsync(dst, src, size, kind, stream)
  #define GPU_MEMCPY_HOST_TO_DEVICE hipMemcpyHostToDevice
  #define GPU_MEMCPY_DEVICE_TO_HOST hipMemcpyDeviceToHost
  #define GPU_MEMCPY_DEVICE_TO_DEVICE hipMemcpyDeviceToDevice
  #define GPU_MEMCPY_HOST_TO_HOST hipMemcpyHostToHost
  #define GPU_MEMSET(ptr, value, size) hipMemset(ptr, value, size)
  #define GPU_MEMSET_ASYNC(ptr, value, size, stream) hipMemsetAsync(ptr, value, size, stream)
  #define GPU_NORMAL_DOUBLE(d_state) hiprand_normal_double(d_state)
  #define GPU_RAND_INIT(seed, subsequence, ffset,state) hiprand_init(seed, subsequence, ffset,state)
  #define GPU_SUCCESS hipSuccess
  #define GPU_RAND_STATE hiprandState_t
  #define GPU_STREAM_T hipStream_t
  #define GPU_STREAM_ADD_CALLBACK(stream, callback, userData, flags)	 hipStreamAddCallback(stream, callback, userData, flags)	
  #define GPU_STREAM_CREATE(stream) hipStreamCreate(stream)
  #define GPU_STREAM_DESTROY(stream) hipStreamDestroy(stream)
  #define GPU_STREAM_SYNC(stream) hipStreamSynchronize(stream)
  #define GPU_STREAM_WAIT_EVENT(stream, event, flags) hipStreamWaitEvent(stream, event, flags)
  #define GPU_RAND_RNGTYPE_T hiprandRngType_t
  #define GPU_RAND_GENERATOR_T hiprandGenerator_t
  #define GPU_RAND_RNG_PSEUDO_DEFAULT HIPRAND_RNG_PSEUDO_DEFAULT
  #define GPU_RAND_DESTROY_GEN(gen) hiprandDestroyGenerator(gen)
  #define GPU_RAND_CREATE_GEN(gen, rngType) hiprandCreateGenerator(gen, rngType)//TODO: &
  #define GPU_RAND_STATUS_SUCCESS HIPRAND_STATUS_SUCCESS
  #define GPU_RAND_SET_PSEUDO_RANDOM_GENERATOR_SEED(gen, seed) hiprandSetPseudoRandomGeneratorSeed(gen, seed)
  #define GPU_RAND_SET_STREAM(gen, stream) hiprandSetStream(gen, stream)
  #define GPU_RAND_GENERATE_NORMAL(generator, outputPtr, n, mean, stddev) hiprandGenerateNormal(generator, outputPtr, n, mean, stddev)
  #define GPU_RAND_GENERATE_NORMAL_DOUBLE(generator, outputPtr, n, mean, stddev) hiprandGenerateNormalDouble(generator, outputPtr, n, mean, stddev)

#elif defined(CUDA_V)
  #include <cuda_runtime.h>
  #define GPU_DEVICE_RESET() cudaDeviceReset()
  #define GPU_DEVICE_SYNCHRONIZE() cudaDeviceSynchronize()
  #define GPU_ERROR_T cudaError_t
  #define GPU_EVENT_T cudaEvent_t
  #define GPU_EVENT_CREATE(event) cudaEventCreate(event)//TODO &
  #define GPU_EVENT_RECORD(event, stream) cudaEventRecord(event, stream)
  #define GPU_FREE(ptr) cudaFree(ptr)
  #define GPU_FREE_HOST(ptr) cudaFreeHost(ptr)
  #define GPU_GET_ERROR_STRING(err) cudaGetErrorString(err)
  #define GPU_GET_LAST_ERROR() cudaGetLastError()
  #define GPU_MALLOC(ptr, size) cudaMalloc(ptr, size)
  #define GPU_MALLOC_HOST(ptr, size) cudaMallocHost(ptr, size)
  #define GPU_MEMCPY(dst, src, size, kind) cudaMemcpy(dst, src, size, kind)
  #define GPU_MEMCPY_ASYNC(dst, src, size, kind, stream) cudaMemcpyAsync(dst, src, size, kind, stream)
  #define GPU_MEMCPY_HOST_TO_DEVICE cudaMemcpyHostToDevice
  #define GPU_MEMCPY_DEVICE_TO_HOST cudaMemcpyDeviceToHost
  #define GPU_MEMCPY_DEVICE_TO_DEVICE cudaMemcpyDeviceToDevice
  #define GPU_MEMCPY_HOST_TO_HOST cudaMemcpyHostToHost
  #define GPU_MEMSET(ptr, value, size) cudaMemset(ptr, value, size)
  #define GPU_MEMSET_ASYNC(ptr, value, size, stream) cudaMemsetAsync(ptr, value, size, stream)
  #define GPU_NORMAL_DOUBLE(d_state) curand_normal_double(d_state)
  #define GPU_RAND_INIT(seed, subsequence, ffset,state) curand_init(seed, subsequence, ffset,state)
  #define GPU_SUCCESS cudaSuccess
  #define GPU_RAND_STATE curandState
  #define GPU_STREAM_T cudaStream_t
  #define GPU_STREAM_ADD_CALLBACK(stream, callback, userData, flags)	 cudaStreamAddCallback(stream, callback, userData, flags)	
  #define GPU_STREAM_CREATE(stream) cudaStreamCreate(stream)
  #define GPU_STREAM_DESTROY(stream) cudaStreamDestroy(stream)
  #define GPU_STREAM_SYNC(stream) cudaStreamSynchronize(stream)
  #define GPU_STREAM_WAIT_EVENT(stream, event, flags) cudaStreamWaitEvent(stream, event, flags)
  #define GPU_RAND_RNGTYPE_T curandRngType_t
  #define GPU_RAND_GENERATOR_T curandGenerator_t
  #define GPU_RAND_RNG_PSEUDO_DEFAULT CUPRAND_RNG_PSEUDO_DEFAULT
  #define GPU_RAND_DESTROY_GEN(gen, rngType) curandDestroyGenerator(gen, rngType)/
  #define GPU_RAND_CREATE_GEN(gen, rngType) curandCreateGenerator(gen, rngType)//TODO: &
  #define GPU_RAND_STATUS_SUCCESS CURAND_STATUS_SUCCESS
  #define GPU_RAND_SET_PSEUDO_RANDOM_GENERATOR_SEED(gen, seed) curandSetPseudoRandomGeneratorSeed(gen, seed)
  #define GPU_RAND_SET_STREAM(gen, stream) curandSetStream(gen, stream)
  #define GPU_RAND_GENERATE_NORMAL(generator, outputPtr, n, mean, stddev) curandGenerateNormal(generator, outputPtr, n, mean, stddev)
  #define GPU_RAND_GENERATE_NORMAL_DOUBLE(generator, outputPtr, n, mean, stddev) curandGenerateNormalDouble(generator, outputPtr, n, mean, stddev)


#endif

/*#define ASSERT_GPU(call) \
  do { \
    GPU_ERROR_T err = call; \
    if (err != GPU_SUCCESS) { \
      fprintf(stderr, "GPU error at %s:%d: %s\n", __FILE__, __LINE__, GPU_GET_ERROR_STRING(err)); \
      exit(1); \
    } \
  } while (0)*/