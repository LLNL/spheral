//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "MurnahanEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class MurnahanEquationOfState<Dim<1> >;
    template class MurnahanEquationOfState<Dim<2> >;
    template class MurnahanEquationOfState<Dim<3> >;
  }
}
