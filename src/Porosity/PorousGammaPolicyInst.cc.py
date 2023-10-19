text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousGammaPolicy.cc"

namespace Spheral {
  template class PorousGammaPolicy<Dim< %(ndim)s > >;
}
"""
