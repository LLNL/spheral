//------------------------------------------------------------------------------
// Provide dummy OpenMP methods for when we compile without OpenMP support.
//------------------------------------------------------------------------------
#ifndef __Spheral_OpenMP_wrapper__
#define __Spheral_OpenMP_wrapper__

#ifdef _OPENMP

#include "omp.h"

#else

inline int  omp_get_num_threads() { return 1; }
inline int  omp_get_thread_num()  { return 0; }
inline void omp_set_num_threads() {}

#endif

//------------------------------------------------------------------------------
// An enum type to help with thread reductions
//------------------------------------------------------------------------------
namespace Spheral {
  enum class ThreadReduction {
    MIN = 0,
    MAX = 1,
    SUM = 2
  };
}

#endif
