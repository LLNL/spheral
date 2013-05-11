//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrainPorosity.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class StrainPorosity<Dim<1> >;
    template class StrainPorosity<Dim<2> >;
    template class StrainPorosity<Dim<3> >;
  }
}
