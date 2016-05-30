text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {
  namespace PartitionSpace {

    template class SpaceFillingCurveRedistributeNodes< Dim< %(ndim)s > >;

  }
}

"""
