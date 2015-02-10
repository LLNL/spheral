//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "InternalNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InternalNodeIterator< Dim<1> >;
  template class InternalNodeIterator< Dim<2> >;
  template class InternalNodeIterator< Dim<3> >;
}
