//------------------------------------------------------------------------------
// Add some safe interfaces to OpenMP from Python.
// Safe here means we can use these regardless of whether OpenMP is on or not.
//------------------------------------------------------------------------------
#ifndef __Spheral_OpenMPhelpers__
#define __Spheral_OpenMPhelpers__

#ifdef _OPENMP
#include "omp.h"
#endif

//------------------------------------------------------------------------------
// omp_get_thread_num
//------------------------------------------------------------------------------
inline
int wrap_omp_get_thread_num() {
  int result = 0;
#ifdef _OPENMP
#pragma omp parallel shared(result)
  {
#pragma omp critical
    result = omp_get_thread_num();
  }
#endif
  return result;
}

//------------------------------------------------------------------------------
// omp_get_num_threads
//------------------------------------------------------------------------------
inline
int wrap_omp_get_num_threads() {
  int result = 1;
#ifdef _OPENMP
#pragma omp parallel shared(result)
  {
#pragma omp critical
    result = omp_get_num_threads();
  }
#endif
  return result;
}

#endif
