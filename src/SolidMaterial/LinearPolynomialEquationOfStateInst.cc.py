text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearPolynomialEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class LinearPolynomialEquationOfState<Dim< %(ndim)s > >;
}
"""
