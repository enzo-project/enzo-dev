#ifndef GPU_UTIL_H_
#define GPU_UTIL_H_

#include "cuda_runtime.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define CUDA_SAFE_CALL(call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
      fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
              __FILE__, __LINE__, cudaGetErrorString( err) );              \
      fflush(stderr); \
      exit(1);                                                  \
    } }

template <class T>
class gpu_array
{
public:
  gpu_array() {
    size_ = 0;
    ptr_ = NULL;
  }
  gpu_array(size_t size) {
    size_ = size;
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
  }
  gpu_array(T* data, size_t size) {
    size_ = size;
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyHostToDevice) );
  }
  ~gpu_array() {
    if (ptr_ != NULL) {
      CUDA_SAFE_CALL( cudaFree(ptr_) );
    }
  }
  void malloc(size_t size) {
    CUDA_SAFE_CALL( cudaMalloc((void**)&ptr_, size*sizeof(T)) );
  }
  void free() {
    if (ptr_ != NULL) {
      CUDA_SAFE_CALL( cudaFree(ptr_) );
      ptr_ = NULL;
    }
  }
  void SetPtr(T *ptr) {
    ptr_ = ptr;
  }
  void D2H(T* data) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(data, ptr_, size_*sizeof(T), cudaMemcpyDeviceToHost) );
  }
  void D2H(T* data, size_t size) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(data, ptr_, size*sizeof(T), cudaMemcpyDeviceToHost) );
  }
  void H2D(T *data) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size_*sizeof(T), cudaMemcpyHostToDevice) );
  }
  void H2D(T *data, size_t size) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyHostToDevice) );
  }
  void D2D(T *data) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size_*sizeof(T), cudaMemcpyDeviceToDevice) );
  }
  void D2D(T *data, size_t size) {
    assert( ptr_ != NULL);
    CUDA_SAFE_CALL( cudaMemcpy(ptr_, data, size*sizeof(T), cudaMemcpyDeviceToDevice) );
  }
  T* ptr() { return ptr_; }
  size_t size() { return size_; }
private:
  T *ptr_;
  size_t size_;
};

inline void cumalloc(void **p, size_t size)
{
  assert(size > 0);
  CUDA_SAFE_CALL( cudaMalloc(p, size) );
}

inline void cufree(void *p)
{
  if (p != NULL)
    CUDA_SAFE_CALL( cudaFree(p) );
}

inline void cuh2df(float *dst, float *src, size_t size)
{
  assert(dst != NULL);
  assert(src != NULL);
  assert(size > 0);
  CUDA_SAFE_CALL( cudaMemcpy(dst, src, size, cudaMemcpyHostToDevice) );
}

inline void cud2hf(float *dst, float *src, size_t size)
{
  assert(dst != NULL);
  assert(src != NULL);
  assert(size > 0);
  CUDA_SAFE_CALL( cudaMemcpy(dst, src, size, cudaMemcpyDeviceToHost) );
}

class CUDATimer
{
protected:
  cudaEvent_t start, stop;
public:
  CUDATimer() {
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
  }
  ~CUDATimer() {
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
  }
  void Start() {
    cudaEventRecord(start, 0);
  }
  double GetET() {
    float gpuTime;
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&gpuTime, start, stop);
    return (double)1e-3*gpuTime;
  }
};

inline void gpu_print_sum(float *a, size_t size)
{
  float *tmp = (float*)malloc(size*sizeof(float));
  cudaMemcpy(tmp, a, size*sizeof(float), cudaMemcpyDeviceToHost);
  float sum = 0;
  for (int i = 0; i < size; i++)
    sum += tmp[i];
  printf("%f ", sum);
  free(tmp);
}

extern "C"
int InitGPU(int ProcessRank);

#endif
