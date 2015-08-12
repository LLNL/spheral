//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MeshPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeshPolicy<Dim<1> >;
  template class MeshPolicy<Dim<2> >;
  template class MeshPolicy<Dim<3> >;
}

