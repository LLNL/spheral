text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/OsborneEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class OsborneEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
