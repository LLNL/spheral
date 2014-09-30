//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPHMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHMassDensityPolicy<Dim<1> >;
  template class SVPHMassDensityPolicy<Dim<2> >;
  template class SVPHMassDensityPolicy<Dim<3> >;
}

