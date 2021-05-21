#include <iostream>


__global__ void add(int a, int b, int *c)
{
  *c = a + b;
}

int main() {
  std::cout << "Hello World\n";

  int a,b,c;
  int *d_c;

  a = 3; b = 4;
  cudaMalloc((void**) &d_c, sizeof(int));
  add<<<1,1>>>(a,b,d_c);
  cudaMemcpy(&c, d_c, sizeof(int), cudaMemcpyDeviceToHost);

  std::cout << "C : " << c << "\n";
  cudaFree(d_c);
  
  return EXIT_SUCCESS;
}
