//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class SpaceFillingCurveRedistributeNodes< Dim<1> >;
    template class SpaceFillingCurveRedistributeNodes< Dim<2> >;
    template class SpaceFillingCurveRedistributeNodes< Dim<3> >;

  }
}

