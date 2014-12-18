//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PorousStrengthModel.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class PorousStrengthModel<Dim<1> >;
    template class PorousStrengthModel<Dim<2> >;
    template class PorousStrengthModel<Dim<3> >;
  }
}
