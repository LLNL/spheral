//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPHCorrectionsPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class SVPHCorrectionsPolicy<Dim<1> >;
  template class SVPHCorrectionsPolicy<Dim<2> >;
  template class SVPHCorrectionsPolicy<Dim<3> >;
}

