#ifndef MEMORYMANAGER
#define MEMORYMANAGER

#if defined(RAJA_ENABLE_CUDA)
#include "RAJA/policy/cuda/raja_cudaerrchk.hpp"
#endif

namespace memoryManager
{

template <typename T>
T *allocate(RAJA::Index_type size)
{
  T *ptr;
#if defined(RAJA_ENABLE_CUDA)
  cudaErrchk(
      cudaMallocManaged((void **)&ptr, sizeof(T) * size, cudaMemAttachGlobal));
#else
  ptr = (T*) malloc(sizeof(T)*size);
#endif
  return ptr;
}

template <typename T>
void deallocate(T *&ptr)
{
  if (ptr) {
#if defined(RAJA_ENABLE_CUDA)
    cudaErrchk(cudaFree(ptr));
#else
    free(ptr);
#endif
    ptr = nullptr;
  }
}

#if defined(RAJA_ENABLE_CUDA)
  template <typename T>
  T *allocate_gpu(RAJA::Index_type size)
  {
    T *ptr;
    cudaErrchk(cudaMalloc((void **)&ptr, sizeof(T) * size));
    return ptr;
  }

  template <typename T>
  void deallocate_gpu(T *&ptr)
  {
    if (ptr) {
      cudaErrchk(cudaFree(ptr));
      ptr = nullptr;
    }
  }
#endif

};  // namespace memoryManager



#endif //  MEMORYMANAGER
