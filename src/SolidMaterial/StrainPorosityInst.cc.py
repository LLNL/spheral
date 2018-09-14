text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/StrainPorosity.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class StrainPorosity<Dim< %(ndim)s > >;
  }
}
"""
