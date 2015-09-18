text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CellPressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CellPressurePolicy<Dim< %(ndim)s > >;
}
"""
