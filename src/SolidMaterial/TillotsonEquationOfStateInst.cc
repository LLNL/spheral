//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TillotsonEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class TillotsonEquationOfState<Dim<1> >;
    template class TillotsonEquationOfState<Dim<2> >;
    template class TillotsonEquationOfState<Dim<3> >;
  }
}
