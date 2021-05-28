#include "DeviceTest.hh"
#include "DeviceTestLib2/DeviceTest2.hh"


namespace Spheral
{

//__device__ void add(int a, int b, int *c)
//{
//  *c = a + b;
//}

#ifdef BUILD_CUDA_TEST_SHARED
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
#endif

} // namespace Spehral
