//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GammaLawGas.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace Material {
    template class GammaLawGas<Dim<1> >;
    template class GammaLawGas<Dim<2> >;
    template class GammaLawGas<Dim<3> >;
  }
}
