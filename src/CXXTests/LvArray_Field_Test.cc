#include "Field/Field.hh"

#include "LvArray/Array.hpp"
#include "LvArray/ChaiBuffer.hpp"

using EXEC_POL = RAJA::cuda_exec<256>;
using HOST_POL = RAJA::seq_exec;

int main() {

  constexpr int N = 10;

  using Dim = Spheral::Dim<3>;
  Spheral::NodeList<Dim> node_list("example node list", N, 0);
  auto n_pos = node_list.positions();

  auto field_view = n_pos.getAccessorView();

  RAJA::RangeSegment range(0, field_view.size());
  RAJA::forall<HOST_POL>(range,
    [=](unsigned int kk) {
      field_view[kk][0]++;
      field_view[kk][1]++;
      field_view[kk][2]++;
  });

  RAJA::forall<EXEC_POL>(range, 
    [=] RAJA_HOST_DEVICE (int kk) {
      field_view[kk][0]++;
      field_view[kk][1]++;
      field_view[kk][2]++;
  });

  bool correctness = true;
  RAJA::forall<HOST_POL>(range,
    [=, &correctness] (int kk) {
      if (field_view[kk] != Spheral::GeomVector<3>(2,2,2)) correctness = false;
  });

  if (correctness)
    std::cout << "PASSED\n";
  else
    std::cout << "FAILED\n";

  return EXIT_SUCCESS;
}
