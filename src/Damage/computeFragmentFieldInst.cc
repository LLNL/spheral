//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "computeFragmentField.cc"

namespace Spheral {
  template Field<Dim<1>, int> computeFragmentField(const NodeList<Dim<1> >&, const double, const Field<Dim<1>, Dim<1>::SymTensor>&, const double, const bool);
  template Field<Dim<2>, int> computeFragmentField(const NodeList<Dim<2> >&, const double, const Field<Dim<2>, Dim<2>::SymTensor>&, const double, const bool);
  template Field<Dim<3>, int> computeFragmentField(const NodeList<Dim<3> >&, const double, const Field<Dim<3>, Dim<3>::SymTensor>&, const double, const bool);
}
