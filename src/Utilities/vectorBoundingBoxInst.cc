//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "boundingBox.cc"

namespace Spheral {
  template void boundingBox(const vector<Dim<1>::Vector>& positions,
                            Dim<1>::Vector& xmin,
                            Dim<1>::Vector& xmax);
  template void boundingBox(const vector<Dim<2>::Vector>& positions,
                            Dim<2>::Vector& xmin,
                            Dim<2>::Vector& xmax);
  template void boundingBox(const vector<Dim<3>::Vector>& positions,
                            Dim<3>::Vector& xmin,
                            Dim<3>::Vector& xmax);
}
