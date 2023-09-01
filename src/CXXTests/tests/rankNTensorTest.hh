#include <iostream>
#include "Geometry/GeomThirdRankTensor.hh"
#include "Geometry/GeomFourthRankTensor.hh"
#include "NodeList/FluidNodeList.hh"

void rankNTensorTest() {

  using DIM = Spheral::Dim<3>;
  constexpr size_t data_sz = 20;

  Spheral::NodeList<DIM> data_node_list("DataNodeList", data_sz, 0);

  Spheral::Field<DIM, Spheral::GeomFourthRankTensor<2>> testf("TRTField", data_node_list);

  std::cout << "sizeof GeomTRT<1> : " << sizeof(Spheral::GeomThirdRankTensor<1>) << std::endl;
  std::cout << "sizeof GeomTRT<2> : " << sizeof(Spheral::GeomThirdRankTensor<2>) << std::endl;
  std::cout << "sizeof GeomTRT<3> : " << sizeof(Spheral::GeomThirdRankTensor<3>) << std::endl;

  testf[0] = Spheral::GeomFourthRankTensor<2>(5.0);
  testf[2] = Spheral::GeomFourthRankTensor<2>(5.0);

  for (std::ptrdiff_t i = 0; i < data_sz; i++) {
    std::cout << testf[i] << "\n";
  }

}
