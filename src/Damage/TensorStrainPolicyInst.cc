//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "TensorStrainPolicy.cc"

namespace Spheral {
  template class Spheral::TensorStrainPolicy<Dim<1> >;
  template class Spheral::TensorStrainPolicy<Dim<2> >;
  template class Spheral::TensorStrainPolicy<Dim<3> >;
}
