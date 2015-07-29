//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CheapSynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class CheapSynchronousRK2< Dim<1> >;
    template class CheapSynchronousRK2< Dim<2> >;
    template class CheapSynchronousRK2< Dim<3> >;
  }
}
