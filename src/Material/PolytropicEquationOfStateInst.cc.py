text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PolytropicEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PolytropicEquationOfState< Dim< %(ndim)s > >;
}
"""
