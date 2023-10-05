text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/YieldStrengthPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class YieldStrengthPolicy<Dim< %(ndim)s > >;
}
"""
