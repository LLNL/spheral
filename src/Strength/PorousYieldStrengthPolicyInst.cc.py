text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousYieldStrengthPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousYieldStrengthPolicy<Dim< %(ndim)s > >;
}
"""
