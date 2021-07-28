#ifndef __Spheral_DeviceTest_hh__
#define __Spheral_DeviceTest_hh__

namespace Spheral
{

#ifdef __CUDACC__
__device__ void add(int a, int b, int *c);

__global__ void launch(int a, int b, int *c);

#else

int launchCaller(int a, int b);

#endif

} // namespace Spehral

#endif // __Spheral_DeviceTest_hh__ 
