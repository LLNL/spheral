text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/StrainPorosity.cc"

namespace Spheral {
  template class StrainPorosity<Dim< %(ndim)s > >;
}
"""
