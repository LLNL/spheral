text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousYieldStrengthPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorousYieldStrengthPolicy<Dim< %(ndim)s > >;
}
"""
