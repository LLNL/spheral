//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SteinbergGuinanLundStrength.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SteinbergGuinanLundStrength<Dim<1> >;
  template class SteinbergGuinanLundStrength<Dim<2> >;
  template class SteinbergGuinanLundStrength<Dim<3> >;
}
