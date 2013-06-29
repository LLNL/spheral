//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousEquationOfState.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousEquationOfState<Dim<1> >;
    template class PorousEquationOfState<Dim<2> >;
    template class PorousEquationOfState<Dim<3> >;
  }
}
