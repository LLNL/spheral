text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/StrainPorosity.cc"

namespace Spheral {
  template class StrainPorosity<Dim< %(ndim)s > >;
}
"""
