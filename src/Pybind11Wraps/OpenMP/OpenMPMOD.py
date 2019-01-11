"""
Spheral OpenMP module.

This module provide thin front-end wrappers for the OpenMP methods.
"""

from PYB11Generator import *

PYB11preamble = """
//------------------------------------------------------------------------------
// Provide dummy OpenMP methods for when we compile without OpenMP support.
//------------------------------------------------------------------------------
#ifdef _OPENMP
#include "omp.h"
inline int wrap_omp_get_num_threads() {
  int result;
  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0) result = omp_get_num_threads();
  }
  return result;
}
#else
inline int wrap_omp_get_num_threads() { return 1; }
inline void     omp_set_num_threads() {}
#endif
"""

@PYB11call_guard("py::gil_scoped_release")
@PYB11cppname("wrap_omp_get_num_threads")
def omp_get_num_threads():
    "Get the number of OpenMP threads."
    return

@PYB11call_guard("py::gil_scoped_release")
def omp_set_num_threads():
    "Set the number of OpenMP threads."
    return
