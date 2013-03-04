//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GhostNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class GhostNodeIterator< Dim<1> >;
  template class GhostNodeIterator< Dim<2> >;
  template class GhostNodeIterator< Dim<3> >;
}
