text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "OsborneEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class OsborneEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
