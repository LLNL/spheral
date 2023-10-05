text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousGammaPolicy.cc"

namespace Spheral {
  template class PorousGammaPolicy<Dim< %(ndim)s > >;
}
"""
