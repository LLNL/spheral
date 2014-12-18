//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "VoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class VoronoiMassDensityPolicy<Dim<1> >;
  template class VoronoiMassDensityPolicy<Dim<2> >;
  template class VoronoiMassDensityPolicy<Dim<3> >;
}

