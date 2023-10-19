text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorousPressurePolicy.cc"

namespace Spheral {
  template class PorousPressurePolicy<Dim< %(ndim)s > >;
}
"""
