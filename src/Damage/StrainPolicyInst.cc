//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StrainPolicy.cc"

namespace Spheral {
  template class Spheral::StrainPolicy<Dim<1> >;
  template class Spheral::StrainPolicy<Dim<2> >;
  template class Spheral::StrainPolicy<Dim<3> >;
}
