#include "DeviceTest.hh"


#include "chai/ManagedArray.hpp"

#include<stdio.h>

namespace Spheral
{

int launchCaller(int a, int b)
{
  chai::ManagedArray<int> c(1);

#ifdef RAJA_ENABLE_CUDA
  using EXEC_POL=RAJA::cuda_exec<256>;
  c.move(chai::ExecutionSpace::GPU);
#else
  using EXEC_POL=RAJA::seq_exec;
  c.move(chai::ExecutionSpace::CPU);
#endif
  
  RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0,1),
    [=] RAJA_HOST_DEVICE (int i) {
      add(a,b,&c[0]);
    }
  );
  c.move(chai::ExecutionSpace::CPU);

  return c[0];
}

} // namespace Spheral
