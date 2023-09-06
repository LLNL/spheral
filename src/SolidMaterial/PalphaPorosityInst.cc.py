text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PalphaPorosity.cc"

namespace Spheral {
  template class PalphaPorosity<Dim< %(ndim)s > >;
}
"""
