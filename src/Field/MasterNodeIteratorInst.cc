//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MasterNodeIterator.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class MasterNodeIterator< Dim<1> >;
  template class MasterNodeIterator< Dim<2> >;
  template class MasterNodeIterator< Dim<3> >;
}
