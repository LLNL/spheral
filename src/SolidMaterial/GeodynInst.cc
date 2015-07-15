//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geodyn.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace SolidMaterial {
    template class Geodyn<Dim<1> >;
    template class Geodyn<Dim<2> >;
    template class Geodyn<Dim<3> >;
  }
}
