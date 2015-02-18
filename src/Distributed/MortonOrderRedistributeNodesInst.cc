//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "MortonOrderRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class MortonOrderRedistributeNodes< Dim<1> >;
    template class MortonOrderRedistributeNodes< Dim<2> >;
    template class MortonOrderRedistributeNodes< Dim<3> >;

  }
}

