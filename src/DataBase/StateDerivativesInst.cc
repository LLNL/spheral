//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StateDerivatives.cc"

namespace Spheral {
  template class StateDerivatives<Dim<1> >;
  template class StateDerivatives<Dim<2> >;
  template class StateDerivatives<Dim<3> >;
}
