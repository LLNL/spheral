#ifndef SPHERAL_CXX_TESTS
#define SPHERAL_CXX_TESTS

#define ATOMIC_ADD(P, V) RAJA::atomicAdd<RAJA::auto_atomic>(P, V)

#include "DeviceTestLib/DeviceTest.hh"

#include <iostream>

void basicLaunchCallerTest()
{
  int a,b,c;
  a = 3; b = 4;

  std::cout << "Testing LaunchCaller()\n";
  c = Spheral::launchCaller(a,b);
  std::cout << "C : " << c << "\n";
}



#endif //  SPHERAL_CXX_TESTS
