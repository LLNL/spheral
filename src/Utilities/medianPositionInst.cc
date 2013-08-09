//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "medianPosition.cc"

namespace Spheral {
  template Dim<1>::Vector medianPosition(vector<Dim<1>::Vector>& positions);
  template Dim<2>::Vector medianPosition(vector<Dim<2>::Vector>& positions);
  template Dim<3>::Vector medianPosition(vector<Dim<3>::Vector>& positions);
}
