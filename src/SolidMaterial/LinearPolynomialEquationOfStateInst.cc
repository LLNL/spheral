//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "LinearPolynomialEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class LinearPolynomialEquationOfState<Dim<1> >;
    template class LinearPolynomialEquationOfState<Dim<2> >;
    template class LinearPolynomialEquationOfState<Dim<3> >;
  }
}
