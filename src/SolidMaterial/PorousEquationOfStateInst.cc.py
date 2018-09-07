text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousEquationOfState.cc"

namespace Spheral {
  template class PorousEquationOfState<Dim< %(ndim)s > >;
}
"""
