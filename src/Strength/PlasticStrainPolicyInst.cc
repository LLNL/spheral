//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PlasticStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PlasticStrainPolicy<Dim<1> >;
  template class PlasticStrainPolicy<Dim<2> >;
  template class PlasticStrainPolicy<Dim<3> >;
}

