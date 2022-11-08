#include <iostream>
#include "DeviceTestLib/DeviceTest.hh"

#define RAJA_ENABLE_DESUL_ATOMICS
#include "RAJA/RAJA.hpp"

#include "Geometry/GeomVector.hh"

#define ATOMIC_ADD RAJA::atomicAdd<RAJA::omp_atomic>

int main() {

#ifdef RAJA_ENABLE_CUDA
  using EXEC_POL=RAJA::cuda_exec<256>;
#else
  using EXEC_POL=RAJA::loop_exec;
#endif

  int a,b,c;
  a = 3; b = 4;

  std::cout << "Testing LaunchCaller()\n";
  c = Spheral::launchCaller(a,b);
  std::cout << "C : " << c << "\n";

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

  
  return EXIT_SUCCESS;
}
