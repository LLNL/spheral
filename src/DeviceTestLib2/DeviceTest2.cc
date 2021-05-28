#include "DeviceTest2.hh"


namespace Spheral
{

__device__ void add(int a, int b, int *c)
{
  *c = a + b;
}


} // namespace Spehral
