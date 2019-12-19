text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/OsborneEquationOfState.cc"

namespace Spheral {
  template class OsborneEquationOfState<Dim< %(ndim)s > >;
}
"""
