#ifndef __Spheral_DeviceTest_hh__
#define __Spheral_DeviceTest_hh__

#include "RAJA/RAJA.hpp"

namespace Spheral
{

RAJA_HOST_DEVICE void add(int a, int b, int *c)
{
  *c = a + b;
}

int launchCaller(int a, int b);

} // namespace Spehral

#endif // __Spheral_DeviceTest_hh__ 
