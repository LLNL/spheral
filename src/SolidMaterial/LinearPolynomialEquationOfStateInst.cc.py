text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearPolynomialEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class LinearPolynomialEquationOfState<Dim< %(ndim)s > >;
  }
}
"""
