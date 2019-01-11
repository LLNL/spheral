text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidMaterial/LinearPolynomialEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class LinearPolynomialEquationOfState<Dim< %(ndim)s > >;
}
"""
