//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrengthModel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class StrengthModel<Dim<1> >;
    template class StrengthModel<Dim<2> >;
    template class StrengthModel<Dim<3> >;
  }
}
