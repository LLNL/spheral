#include <iostream>
#include "DeviceTestLib/DeviceTest.hh"

#include "RAJA/RAJA.hpp"

int main() {
  int a,b,c;
  a = 3; b = 4;

  c = Spheral::launchCaller(a,b);

  std::cout << "C : " << c << "\n";
  
  return EXIT_SUCCESS;
}
