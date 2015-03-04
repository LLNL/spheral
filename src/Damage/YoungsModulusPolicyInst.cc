//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "YoungsModulusPolicy.cc"

namespace Spheral {
  template class YoungsModulusPolicy<Dim<1> >;
  template class YoungsModulusPolicy<Dim<2> >;
  template class YoungsModulusPolicy<Dim<3> >;
}

