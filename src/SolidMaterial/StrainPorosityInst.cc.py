text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrainPorosity.cc"

namespace Spheral {
  template class StrainPorosity<Dim< %(ndim)s > >;
}
"""
