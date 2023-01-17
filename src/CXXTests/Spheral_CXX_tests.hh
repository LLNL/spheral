#ifndef SPHERAL_CXX_TESTS
#define SPHERAL_CXX_TESTS

#include "DeviceTestLib/DeviceTest.hh"
#include "RAJA/RAJA.hpp"
#include "desul/atomics.hpp"

#define ORDER desul::MemoryOrderRelaxed;
#define SCOPE desul::MemoryScopeDevice;

//*****************************************************************************
//
//*****************************************************************************
#include <iostream>

void basicLaunchCallerTest()
{
  int a,b,c;
  a = 3; b = 4;

  std::cout << "Testing LaunchCaller()\n";
  c = Spheral::launchCaller(a,b);
  std::cout << "C : " << c << "\n";
}




//*****************************************************************************
//
//*****************************************************************************
#define ATOMIC_ADD(P, V) RAJA::atomicAdd<RAJA::auto_atomic>(P, V)
//#define ATOMIC_ADD(P, V) desul::atomic_add(P, V, desul::MemoryOrderRelease{}, desul::MemoryScopeDevice{})
#include "Geometry/GeomVector.hh"

void basicAtomicTets()
{

#ifdef RAJA_ENABLE_CUDA
  using EXEC_POL=RAJA::cuda_exec<256>;
#else
  using EXEC_POL=RAJA::loop_exec;
#endif

  std::cout << "Testing RAJA Atomic w/ GeomVector\n";

  using SPH_TYPE = Spheral::GeomVector<1>;
  SPH_TYPE vec1(0);

  RAJA::AtomicRef<SPH_TYPE, RAJA::auto_atomic> vec1_ref(&vec1);

  RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0,7), [=] RAJA_HOST_DEVICE (int i) {
      int d_a, d_b, d_c;
      d_a = 1; d_b = 3;
      Spheral::add(d_a,d_b,&d_c);


      SPH_TYPE inc(d_c);
      vec1_ref += inc;
  });

  std::cout << vec1 << "\n";
}


#include "tests/SPHEvalDerivTest.hh"

#endif //  SPHERAL_CXX_TESTS
