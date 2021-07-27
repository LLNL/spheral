#include "DeviceTest.hh"

namespace Spheral
{

__device__ void add(int a, int b, int *c)
{
  *c = a + b;
}

__global__ void launch(int a, int b, int *c)
{
  add(a,b,c);
}

__host__ int launchCaller(int a, int b)
{
  int c;
  int *d_c;
  cudaMalloc((void**) &d_c, sizeof(int));
  launch<<<1,1>>>(a,b,d_c);
  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
  return c;
}

} // namespace Spehral
