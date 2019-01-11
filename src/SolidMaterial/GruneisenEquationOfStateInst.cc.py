text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GruneisenEquationOfState.cc"

namespace Spheral {
  template class GruneisenEquationOfState<Dim< %(ndim)s > >;
}
"""
