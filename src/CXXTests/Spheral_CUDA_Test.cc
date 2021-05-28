#include <iostream>
#include "DeviceTestLib/DeviceTest.hh"

#ifndef BUILD_CUDA_TEST_SHARED

__global__ void launch(int a, int b, int *c)
{
  Spheral::add(a,b,c);
}

#endif


int main() {
  std::cout << "Hello World\n";

  int a,b,c;

  a = 3; b = 4;

#ifdef BUILD_CUDA_TEST_SHARED
  Spheral::launchCaller(a,b,c);
#else
  int *d_c;
  cudaMalloc((void**) &d_c, sizeof(int));
  launch<<<1,1>>>(a,b,d_c);
  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  cudaFree(d_c);
#endif

  std::cout << "C : " << c << "\n";
  
  return EXIT_SUCCESS;
}
