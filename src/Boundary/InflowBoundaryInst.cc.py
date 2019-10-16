text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "InflowBoundary.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InflowBoundary< Dim< %(ndim)s > >;
}
"""
