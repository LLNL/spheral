//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "EquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class EquationOfState< Dim<1> >;
    template class EquationOfState< Dim<2> >;
    template class EquationOfState< Dim<3> >;
  }
}

