#include "DeviceTest.hh"

namespace Spheral
{

#if defined(RAJA_ENABLE_CUDA)
__device__ void add(int a, int b, int *c)
{
  *c = a + b;
}

__global__ void launch(int a, int b, int *c)
{
  add(a,b,c);
}

int launchCaller(int a, int b)
{
  int c;
  int *d_c;
  cudaMalloc((void**) &d_c, sizeof(int));
  launch<<<1,1>>>(a,b,d_c);
  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
  return c;
}
#else

int launchCaller(int a, int b) {return a + b;}
#endif

} // namespace Spehral
