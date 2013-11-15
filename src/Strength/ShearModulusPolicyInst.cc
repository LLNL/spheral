//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ShearModulusPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ShearModulusPolicy<Dim<1> >;
  template class ShearModulusPolicy<Dim<2> >;
  template class ShearModulusPolicy<Dim<3> >;
}

