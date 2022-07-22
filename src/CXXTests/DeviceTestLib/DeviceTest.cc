#include "DeviceTest.hh"

#include<stdio.h>

namespace Spheral
{

#ifdef __CUDACC__
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
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) 
      printf("Error: %s\n", cudaGetErrorString(err));

  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
  return c;
}

#else
int launchCaller(int a, int b)
{
  return a + b;
}

#endif

} // namespace Spehral
