//------------------------------------------------------------------------------
// Add some safe interfaces to OpenMP from Python.
// Safe here means we can use these regardless of whether OpenMP is on or not.
//------------------------------------------------------------------------------
#ifdef _OPENMP

#include <omp.h>

#else

int omp_get_thread_num()  { return 0; }
int omp_get_num_threads() { return 1; }

#endif

