text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PalphaPorosity.cc"

namespace Spheral {
  template class PalphaPorosity<Dim< %(ndim)s > >;
}
"""
