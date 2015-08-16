text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrainPorosity.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class StrainPorosity<Dim< %(ndim)s > >;
  }
}
"""
