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
  int *d_c;

  a = 3; b = 4;
  cudaMalloc((void**) &d_c, sizeof(int));

#ifdef BUILD_CUDA_TEST_SHARED
  Spheral::launchCaller(a,b,d_c);
#else
  launch<<<1,1>>>(a,b,d_c);
#endif

  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);
  std::cout << "C : " << c << "\n";
  cudaFree(d_c);
  
  return EXIT_SUCCESS;
}
