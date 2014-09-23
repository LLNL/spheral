//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "IsothermalEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class IsothermalEquationOfState<Dim<1> >;
    template class IsothermalEquationOfState<Dim<2> >;
    template class IsothermalEquationOfState<Dim<3> >;
  }
}
