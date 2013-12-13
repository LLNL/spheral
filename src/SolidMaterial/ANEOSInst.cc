//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ANEOS.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class ANEOS<Dim<1> >;
    template class ANEOS<Dim<2> >;
    template class ANEOS<Dim<3> >;
  }
}
