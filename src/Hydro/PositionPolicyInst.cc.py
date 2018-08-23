text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/PositionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PositionPolicy<Dim< %(ndim)s > >;
}
"""
