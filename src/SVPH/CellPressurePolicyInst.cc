//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CellPressurePolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CellPressurePolicy<Dim<1> >;
  template class CellPressurePolicy<Dim<2> >;
  template class CellPressurePolicy<Dim<3> >;
}
