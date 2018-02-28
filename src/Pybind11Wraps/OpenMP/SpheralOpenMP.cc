// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"

#ifdef _OPENMP
#include "omp.h"
#else
inline int omp_get_thread_num() { return 0; }
inline int omp_get_num_threads() { return 1; }
#endif

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralOpenMP, m) {
  namespace py = pybind11;
  using namespace pybind11::literals;

  m.doc() = "Spheral OpenMP wrappers.";

  m.def("omp_get_thread_num", &omp_get_thread_num, "Get the OpenMP thread ID.");
  m.def("omp_get_num_threads", &omp_get_num_threads, "Get the number of OpenMP threads.");
}
