#include <iostream>
#include "DeviceTestLib/DeviceTest.hh"

#define RAJA_ENABLE_DESUL_ATOMICS
#include "RAJA/RAJA.hpp"

#include "Geometry/GeomVector.hh"

#define ATOMIC_ADD RAJA::atomicAdd<RAJA::omp_atomic>

class MyTestData {
  double data;
public:
  MyTestData() : data(0) {};
  MyTestData(double d) : data(d) {};
  MyTestData(const MyTestData&) = default;

  //MyTestData& operator+=(const MyTestData& rhs) {
  //  this->data += rhs.data;
  //  return *this;
  //}

  friend std::ostream& operator<<(std::ostream& os, const MyTestData& rhs);
  friend MyTestData operator+(const MyTestData& lhs, const MyTestData& rhs);

  //---------------------------------------------------------------------------
  //MyTestData(const volatile MyTestData& rhs) {
  //  data = rhs.data; 
  //};

  //MyTestData& operator=(const volatile MyTestData& rhs) {
  //  this->data = rhs.data;
  //  return *this;
  //}

  //volatile MyTestData& operator+=(const MyTestData& rhs) volatile {
  //  this->data += rhs.data;
  //  return *this;
  //}
  //---------------------------------------------------------------------------

};

MyTestData operator+(const MyTestData& lhs, const MyTestData& rhs) {
  return MyTestData(lhs.data + rhs.data);
}
std::ostream& operator<<(std::ostream& os, const MyTestData& rhs) {
  os << rhs.data;
  return os;
}

int main() {
  //int a,b,c;
  //a = 3; b = 4;

  //c = Spheral::launchCaller(a,b);

  //using SPH_TYPE = double;
  using SPH_TYPE = Spheral::GeomVector<1>;
  //using SPH_TYPE = MyTestData;


  SPH_TYPE vec1(0);
  RAJA::AtomicRef<SPH_TYPE, RAJA::auto_atomic> vec1_ref(&vec1);

  RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0,100), [=](int i) {
      
      SPH_TYPE inc(1);

      vec1_ref += inc;
      //ATOMIC_ADD(vec1, inc);
      //vec1 += inc;

      });

  std::cout << vec1 << "\n";

  //std::cout << "C : " << c << "\n";
  
  return EXIT_SUCCESS;
}
