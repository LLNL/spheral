//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SumVoronoiMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SumVoronoiMassDensityPolicy<Dim<1> >;
  template class SumVoronoiMassDensityPolicy<Dim<2> >;
  template class SumVoronoiMassDensityPolicy<Dim<3> >;
}

