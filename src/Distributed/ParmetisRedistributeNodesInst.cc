//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ParmetisRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class ParmetisRedistributeNodes< Dim<2> >;
    template class ParmetisRedistributeNodes< Dim<3> >;
  }
}

