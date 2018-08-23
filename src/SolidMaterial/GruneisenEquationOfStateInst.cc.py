text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/GruneisenEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class GruneisenEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
