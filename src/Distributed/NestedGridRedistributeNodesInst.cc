//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "NestedGridRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class NestedGridRedistributeNodes< Dim<1> >;
    template class NestedGridRedistributeNodes< Dim<2> >;
    template class NestedGridRedistributeNodes< Dim<3> >;
  }
}

