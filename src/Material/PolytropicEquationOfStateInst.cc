//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PolytropicEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class PolytropicEquationOfState<Dim<1> >;
    template class PolytropicEquationOfState<Dim<2> >;
    template class PolytropicEquationOfState<Dim<3> >;
  }
}
