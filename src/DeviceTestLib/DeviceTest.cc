#include "DeviceTest.hh"


namespace Spheral
{

__device__ void add(int a, int b, int *c)
{
  *c = a + b;
}

#ifdef BUILD_CUDA_TEST_SHARED
__global__ void launch(int a, int b, int *c)
{
  add(a,b,c);
}

__host__ void launchCaller(int a, int b, int *c)
{
  launch<<<1,1>>>(a,b,c);
}
#endif

} // namespace Spehral
