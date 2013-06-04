//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "RedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class RedistributeNodes< Dim<1> >;
    template class RedistributeNodes< Dim<2> >;
    template class RedistributeNodes< Dim<3> >;
  }
}

