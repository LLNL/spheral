text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PressurePolicy<Dim< %(ndim)s > >;
}
"""
