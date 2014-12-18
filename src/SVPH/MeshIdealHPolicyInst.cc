//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MeshIdealHPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MeshIdealHPolicy<Dim<1> >;
  template class MeshIdealHPolicy<Dim<2> >;
  template class MeshIdealHPolicy<Dim<3> >;
}
