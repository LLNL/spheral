text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "HelmholtzEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class HelmholtzEquationOfState< Dim< %(ndim)s > >;
}
"""
