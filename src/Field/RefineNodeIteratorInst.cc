//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "RefineNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class RefineNodeIterator< Dim<1> >;
  template class RefineNodeIterator< Dim<2> >;
  template class RefineNodeIterator< Dim<3> >;
}
