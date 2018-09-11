text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/PorousEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
