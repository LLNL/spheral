#ifndef __Spheral_DeviceTest_hh__
#define __Spheral_DeviceTest_hh__

#define BUILD_CUDA_TEST_SHARED

namespace Spheral
{

//__device__ void add(int a, int b, int *c);

#ifdef BUILD_CUDA_TEST_SHARED
__global__ void launch(int a, int b, int *c);

__host__ int launchCaller(int a, int b);
#endif

} // namespace Spehral

#endif // __Spheral_DeviceTest_hh__ 
