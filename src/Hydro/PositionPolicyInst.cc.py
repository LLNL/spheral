text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PositionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PositionPolicy<Dim< %(ndim)s > >;
}
"""
