//------------------------------------------------------------------------------
// Explicit instation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "StateBase.cc"

namespace Spheral {
  template class StateBase<Dim<1> >;
  template class StateBase<Dim<2> >;
  template class StateBase<Dim<3> >;
}
