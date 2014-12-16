//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeIteratorBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NodeIteratorBase< Dim<1> >;
  template class NodeIteratorBase< Dim<2> >;
  template class NodeIteratorBase< Dim<3> >;
}
