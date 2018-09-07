text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "OsborneEquationOfState.cc"

namespace Spheral {
  template class OsborneEquationOfState<Dim< %(ndim)s > >;
}
"""
