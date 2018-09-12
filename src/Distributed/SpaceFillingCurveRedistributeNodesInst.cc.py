text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {
  template class SpaceFillingCurveRedistributeNodes< Dim< %(ndim)s > >;
}

"""
