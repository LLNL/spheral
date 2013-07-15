//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GruneisenEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class GruneisenEquationOfState<Dim<1> >;
    template class GruneisenEquationOfState<Dim<2> >;
    template class GruneisenEquationOfState<Dim<3> >;
  }
}
