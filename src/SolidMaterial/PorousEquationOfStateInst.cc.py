text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PorousEquationOfState.cc"

namespace Spheral {
  template class PorousEquationOfState<Dim< %(ndim)s > >;
}
"""
