//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Verlet.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class Verlet< Dim<1> >;
    template class Verlet< Dim<2> >;
    template class Verlet< Dim<3> >;
  }
}
