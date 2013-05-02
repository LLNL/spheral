//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ANEOSEquationOfState.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class ANEOSEquationOfState<Dim<1> >;
    template class ANEOSEquationOfState<Dim<2> >;
    template class ANEOSEquationOfState<Dim<3> >;
  }
}
