text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "SpaceFillingCurveRedistributeNodes.cc"

namespace Spheral {
  template class SpaceFillingCurveRedistributeNodes< Dim< %(ndim)s > >;
}

"""
