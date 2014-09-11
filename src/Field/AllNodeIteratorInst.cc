//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "AllNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class AllNodeIterator< Dim<1> >;
  template class AllNodeIterator< Dim<2> >;
  template class AllNodeIterator< Dim<3> >;
}
