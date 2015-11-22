text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GruneisenEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class GruneisenEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
