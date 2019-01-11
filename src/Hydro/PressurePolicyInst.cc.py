text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Hydro/PressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PressurePolicy<Dim< %(ndim)s > >;
}
"""
