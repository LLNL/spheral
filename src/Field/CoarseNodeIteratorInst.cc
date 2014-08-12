//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CoarseNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CoarseNodeIterator< Dim<1> >;
  template class CoarseNodeIterator< Dim<2> >;
  template class CoarseNodeIterator< Dim<3> >;
}
