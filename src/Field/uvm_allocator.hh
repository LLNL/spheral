#ifndef __UVM_ALLOCATOR_H
#define __UVM_ALLOCATOR_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <cassert>

namespace uvm_allocator
{
template <class T>
struct UVMAllocator{
  typedef T value_type;

  // Constructor

  UVMAllocator() noexcept {}

  // Copy constructor, necessary to create new value types based on old value types, e.g. T(a) = b

  template <class U> UVMAllocator(const UVMAllocator<U>& other) noexcept {}

  // allocate()

  T* allocate(const std::size_t n){
    void *pt;
    // allocate CUDA managed memory
    cudaMallocManaged(&pt, n*sizeof(T));
    // The following, new(pt), initializes the allocated data using the type's
    // default constructor rather than creating another allocation.
    return new(pt) T[n];
  }

  // deallocate()

  void deallocate(T* ptr, const std::size_t n){
    assert( ptr != nullptr );
    for (std::size_t i=0; i < n; ++i)
       ptr[i].~T();    //calling default destructor
    cudaFree(ptr);     //freeing memory

  }

};

template <class T, class U>
constexpr bool operator== (const UVMAllocator<T>&, const UVMAllocator<U>&) noexcept
{return true;}

template <class T, class U>
constexpr bool operator!= (const UVMAllocator<T>&, const UVMAllocator<U>&) noexcept
{return false;}

}; // end of namespace uvm_allocator

#endif