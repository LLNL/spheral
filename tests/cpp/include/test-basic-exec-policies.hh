#ifndef SPHERAL_BASIC_EXEC_POL_HH
#define SPHERAL_BASIC_EXEC_POL_HH

#include "gtest/gtest.h"
#include <RAJA/policy/sequential/policy.hpp>

using SEQ_EXEC_POLICY = RAJA::seq_exec;
// The list of execution types we want to possibly run these tests over.
using EXEC_TYPES = ::testing::Types<SEQ_EXEC_POLICY
#ifdef SPHERAL_ENABLE_CUDA
                                    ,
                                    RAJA::cuda_exec<512>
#endif
#ifdef SPHERAL_ENABLE_HIP
                                    ,
                                    RAJA::hip_exec<512>
#endif
                                    >;

#endif // SPHERAL_BASIC_EXEC_POL_HH
