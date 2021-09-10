#ifndef __Spheral_DeviceTest_hh__
#define __Spheral_DeviceTest_hh__

namespace Spheral
{

#if defined(RAJA_ENABLE_CUDA)
__device__ void add(int a, int b, int *c);

__global__ void launch(int a, int b, int *c);
#endif

int launchCaller(int a, int b);

} // namespace Spehral

#endif // __Spheral_DeviceTest_hh__ 
