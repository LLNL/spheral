text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Boundary/InflowOutflowBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InflowOutflowBoundary< Dim< %(ndim)s > >;
}
"""
