#ifndef SPHERAL_BASIC_EXEC_POL_HH
#define SPHERAL_BASIC_EXEC_POL_HH

#include "RAJA/RAJA.hpp"

// The list of execution types we want to possibly run these tests over.
using EXEC_TYPES = ::testing::Types<
  RAJA::loop_exec
#ifdef SPHERAL_ENABLE_CUDA
  ,RAJA::cuda_exec<512>
#endif
>;

#endif // SPHERAL_BASIC_EXEC_POL_HH
