//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SynchronousRK2.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class SynchronousRK2< Dim<1> >;
    template class SynchronousRK2< Dim<2> >;
    template class SynchronousRK2< Dim<3> >;
  }
}
