text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
