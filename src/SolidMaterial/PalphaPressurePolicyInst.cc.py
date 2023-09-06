text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PalphaPressurePolicy.cc"

namespace Spheral {
  template class PalphaPressurePolicy<Dim< %(ndim)s > >;
}
"""
