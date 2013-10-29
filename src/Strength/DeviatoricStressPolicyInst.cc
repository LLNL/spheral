//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "DeviatoricStressPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class DeviatoricStressPolicy<Dim<1> >;
  template class DeviatoricStressPolicy<Dim<2> >;
  template class DeviatoricStressPolicy<Dim<3> >;
}
