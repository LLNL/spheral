#ifndef __UVM_ALLOCATOR_H
#define __UVM_ALLOCATOR_H

//#include <cuda.h>
//#include <cuda_runtime.h>
#include <cstdlib>
#include <memory> 

namespace uvm_allocator
{
template <class T>
struct UVMAllocator{

   using value_type=T;
   using pointer = value_type*;
   
   
  // Constructor

  UVMAllocator() = default; 

  // Copy constructor, necessary to create new value types based on old value types, e.g. T(a) = b

  template <class U> UVMAllocator(const UVMAllocator<U>& other) noexcept; 

  // allocate()

  T* allocate(std::size_t n){
    void *pt;
    // allocate CUDA managed memory
//    cudaMallocManaged(&pt, n*sizeof(T));
     pt = std::malloc( n*sizeof(T) );
    // The following, new(pt), initializes the allocated data using the type's
    // default constructor rather than creating another allocation.
    return new(pt) T[n];
  }

  // deallocate()

  void deallocate(T* ptr, std::size_t n){
    for (std::size_t i=0; i < n; ++i)
       ptr[i].~T();    //calling default destructor
//    cudaFree(ptr);     //freeing memory
	std:free((void*)ptr);
  }

};

template <class T, class U>
constexpr bool operator== (const UVMAllocator<T>&, const UVMAllocator<U>&) noexcept {return true;};

template <class T, class U>
constexpr bool operator!= (const UVMAllocator<T>&, const UVMAllocator<U>&) noexcept {return false;};

}; // end of namespace uvm_allocator

#endif
