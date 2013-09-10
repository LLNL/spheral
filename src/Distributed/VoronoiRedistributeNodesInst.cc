//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "VoronoiRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class VoronoiRedistributeNodes< Dim<1> >;
    template class VoronoiRedistributeNodes< Dim<2> >;
    template class VoronoiRedistributeNodes< Dim<3> >;

  }
}

