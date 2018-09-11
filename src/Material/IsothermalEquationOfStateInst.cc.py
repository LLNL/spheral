text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Material/IsothermalEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class IsothermalEquationOfState< Dim< %(ndim)s > >;
  }
}
"""
