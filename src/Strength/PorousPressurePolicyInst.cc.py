text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Strength/PorousPressurePolicy.cc"

namespace Spheral {
  template class PorousPressurePolicy<Dim< %(ndim)s > >;
}
"""
