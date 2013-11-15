//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "PeanoHilbertOrderRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class PeanoHilbertOrderRedistributeNodes< Dim<1> >;
    template class PeanoHilbertOrderRedistributeNodes< Dim<2> >;
    template class PeanoHilbertOrderRedistributeNodes< Dim<3> >;

  }
}

