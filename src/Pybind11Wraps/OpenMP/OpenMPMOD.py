"""
Spheral OpenMP module.

This module provide thin front-end wrappers for the OpenMP methods.
"""

preamble = """
//------------------------------------------------------------------------------
// Provide dummy OpenMP methods for when we compile without OpenMP support.
//------------------------------------------------------------------------------
#ifdef _OPENMP
#include "omp.h"
#else
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif
"""

def omp_get_thread_num():
    "Get the OpenMP thread ID."
    return

def omp_get_num_threads():
    "Get the number of OpenMP threads."
    return
