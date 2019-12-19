text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "InflowOutflowBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InflowOutflowBoundary< Dim< %(ndim)s > >;
}
"""
