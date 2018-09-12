text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/GruneisenEquationOfState.cc"

namespace Spheral {
  template class GruneisenEquationOfState<Dim< %(ndim)s > >;
}
"""
