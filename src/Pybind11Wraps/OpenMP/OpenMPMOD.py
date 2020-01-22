"""
Spheral OpenMP module.

This module provide thin front-end wrappers for the OpenMP methods.
"""

from PYB11Generator import *

PYB11includes = ['"Utilities/OpenMP_wrapper.hh"']

PYB11preamble = """
#ifdef _OPENMP
//------------------------------------------------------------------------------
// In order to report the number of threads, we have to check in an OpenMP 
// section.
//------------------------------------------------------------------------------
inline int wrap_omp_get_num_threads() {
  int result;
  #pragma omp parallel
  {
    if (omp_get_thread_num() == 0) result = omp_get_num_threads();
  }
  return result;
}
#endif
"""

@PYB11cppname("wrap_omp_get_num_threads")
def omp_get_num_threads():
    "Get the number of OpenMP threads."
    return

def omp_set_num_threads():
    "Set the number of OpenMP threads."
    return
