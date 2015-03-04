//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SynchronousRK4.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace IntegratorSpace {
    template class SynchronousRK4< Dim<1> >;
    template class SynchronousRK4< Dim<2> >;
    template class SynchronousRK4< Dim<3> >;
  }
}
