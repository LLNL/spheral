text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class SpaceFillingCurveRedistributeNodes< Dim< %(ndim)s > >;

  }
}

"""
