text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/PolytropicEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PolytropicEquationOfState< Dim< %(ndim)s > >;
}
"""
