//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DistributeByXPosition.cc"

namespace Spheral {
  namespace PartitionSpace {
    template class DistributeByXPosition< Dim<1> >;
    template class DistributeByXPosition< Dim<2> >;
    template class DistributeByXPosition< Dim<3> >;
  }
}

